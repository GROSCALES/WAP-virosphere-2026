#!/bin/bash
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=80G
#SBATCH --job-name=defensas_votus
#SBATCH --output=defensas_votus_%j.log

hostname
date
set -uo pipefail

# ==========================================================================
#  Cribado de sistemas de defensa sobre el CATALOGO DE vOTUs
#  PADLOC + DefenseFinder
#
#  Entrada: UN FASTA multi-contig (no 132 genomas como en la version MAGs).
#  El paralelismo sale de trocear el FASTA, no de un bucle de genomas.
#
#  PREDICCION DE GENES
#   prodigal -p meta se corre UNA vez por chunk y alimenta a las dos
#   herramientas, de modo que ven exactamente los mismos ORFs.
#
#   No se usa 'padloc --fna': su prodigal interno aborta con contigs
#   cortos ("Prodigal works best with sequences > 100kbp"), y los vOTUs
#   son mayoritariamente cortos.
#
#   Con --faa/--gff hay que reescribir el GFF: prodigal escribe
#   ID=<num_contig>_<num_gen> (p.ej. ID=1_2) mientras que las cabeceras
#   del .faa son <nombre_contig>_<num_gen>. PADLOC exige que coincidan,
#   y si no, aborta con "errexit on line 425". El bloque awk de abajo
#   reescribe ID= usando la columna 1 del GFF (nombre real del contig).
#
#  CRUCE POSTERIOR CON VIBRANT
#   Hacerlo por CONTIG + COORDENADAS, nunca por el string del ID de
#   proteina: los profagos escindidos de VIBRANT llevan otro marco de
#   coordenadas y otra numeracion.
#
#  INTERPRETAR LA SALIDA DE DefenseFinder
#   'systems.tsv' solo lista SISTEMAS COMPLETOS. Los fagos llevan
#   componentes sueltos, asi que suele venir vacio. El fichero util es
#   'hmmer.tsv', con todos los hits de perfil HMM.
#
#  PILOTO:  PILOT=1 sbatch 08_defensas_votus.sh
#           Comprueba al final la linea "replicones": si se acerca al
#           numero de contigs con hit, gembase parsea bien y puedes
#           lanzar el run completo quitando PILOT=1.
# ==========================================================================

# ---------- 1. Rutas ----------
WAP_BASE="${WAP_BASE:-$HOME/gabi/vOTUs_total}"
IN_FASTA="${IN_FASTA:-$WAP_BASE/vOTUs_catalogue_rep_seq.fasta}"
OUT_DIR="${OUT_DIR:-$WAP_BASE/analisis_defensas_votus}"

# entornos conda
ENV_DF="defensefinder_env"
ENV_PADLOC="padloc"

# prodigal vive dentro del entorno de padloc; ruta absoluta para evitar
# el overhead y la captura de salida de 'conda run'
PRODIGAL_BIN="${PRODIGAL_BIN:-$HOME/miniconda3/envs/padloc/bin/prodigal}"
export PRODIGAL_BIN

PILOT="${PILOT:-0}"
PILOT_N=200

mkdir -p "$OUT_DIR"
WORK="$OUT_DIR/work"
LOGS="$OUT_DIR/logs"
mkdir -p "$WORK" "$LOGS"
ERROR_LOG="$OUT_DIR/errores.log"
: > "$ERROR_LOG"

if [[ ! -s "$IN_FASTA" ]]; then
    echo "ERROR: no encuentro $IN_FASTA" | tee -a "$ERROR_LOG"
    exit 1
fi
if [[ ! -x "$PRODIGAL_BIN" ]]; then
    echo "ERROR: no encuentro prodigal en $PRODIGAL_BIN" | tee -a "$ERROR_LOG"
    exit 1
fi

# ---------- 2. Cores ----------
CORES_TOTAL="${SLURM_CPUS_PER_TASK:-30}"
PARALLEL_JOBS=6
THREADS_PER_JOB=$(( CORES_TOTAL / PARALLEL_JOBS ))
[[ "$THREADS_PER_JOB" -lt 1 ]] && THREADS_PER_JOB=1

echo "=========================================================="
echo " CRIBADO DE DEFENSAS EN vOTUs"
echo " Entrada : $IN_FASTA"
echo " Salida  : $OUT_DIR"
echo " Cores   : $CORES_TOTAL | chunks en paralelo: $PARALLEL_JOBS | hilos/chunk: $THREADS_PER_JOB"
[[ "$PILOT" == "1" ]] && echo " MODO PILOTO: solo los primeros $PILOT_N contigs"
echo "=========================================================="

# ---------- 3. Versiones (para Metodos: apuntalo ahora o se pierde) ----------
VERS="$OUT_DIR/VERSIONS.txt"
{
    echo "fecha: $(date -Iseconds)"
    echo "host : $(hostname)"
    echo "input: $IN_FASTA"
    echo "md5  : $(md5sum "$IN_FASTA" | cut -d' ' -f1)"
    echo "--- prodigal ---"
    "$PRODIGAL_BIN" -v 2>&1 | head -3
    echo "--- padloc ---"
    conda run -n "$ENV_PADLOC" padloc --version 2>&1 | head -5
    echo "--- padloc db ---"
    conda run -n "$ENV_PADLOC" padloc --db-version 2>&1 | head -5
    echo "--- defense-finder ---"
    conda run -n "$ENV_DF" pip show mdmparis-defense-finder 2>&1 | head -3
    conda run -n "$ENV_DF" macsyfinder --version 2>&1 | head -3
} > "$VERS" 2>&1
echo "[i] Versiones guardadas en $VERS"

# ---------- 4. Recuento y troceado ----------
N_CONTIGS=$(grep -c '^>' "$IN_FASTA")
echo "[i] Contigs en el catalogo: $N_CONTIGS"

SRC="$IN_FASTA"
if [[ "$PILOT" == "1" ]]; then
    SRC="$WORK/pilot.fna"
    awk -v n="$PILOT_N" '/^>/{c++} c>n{exit} {print}' "$IN_FASTA" > "$SRC"
    echo "[i] Piloto: $(grep -c '^>' "$SRC") contigs"
fi

rm -f "$WORK"/chunk_*.fna
awk -v n="$PARALLEL_JOBS" -v pre="$WORK/chunk" '
    /^>/ { i++; c = ((i-1) % n) + 1 }
    { print > (pre "_" c ".fna") }
' "$SRC"
echo "[i] Troceado en $(ls "$WORK"/chunk_*.fna 2>/dev/null | wc -l) chunks"

# ---------- 5. Funcion que procesa UN chunk ----------
procesar_chunk() {
    local fna="$1" threads="$2" out_dir="$3" error_log="$4" env_df="$5" env_padloc="$6"
    local tag faa gff gffp cdir n_gff n_common

    tag=$(basename "$fna" .fna)
    cdir="$out_dir/chunks/$tag"
    mkdir -p "$cdir/padloc" "$cdir/defense_finder"

    faa="$cdir/$tag.faa"
    gff="$cdir/$tag.gff"          # gff crudo de prodigal (ID=1_2)
    gffp="$cdir/${tag}_padloc.gff" # gff con ID=<contig>_2, para PADLOC

    # ---- 5.1 Prediccion de genes ----
    echo ">> $tag : prediccion de genes (prodigal -p meta)"
    if ! "$PRODIGAL_BIN" -i "$fna" -a "$faa" -f gff -o "$gff" -p meta -q 2> "$cdir/prodigal.err" ; then
        echo "ERROR prodigal en $tag :: $(tail -3 "$cdir/prodigal.err" | tr '\n' ' ')" | tee -a "$error_log"
        return
    fi
    # cabeceras del faa: quitar ' # start # end # strand # attrs'
    # deja <contig>_<n>, que es lo que MacSyFinder necesita para gembase
    sed -i 's/ #.*//' "$faa"

    # ---- 5.2 Reescribir ID= del GFF para que case con el faa ----
    awk 'BEGIN{FS=OFS="\t"}
         /^#/{print; next}
         NF>=9{
             n=$9; sub(/^.*ID=/,"",n); sub(/;.*/,"",n);      # n = "1_2"
             g=n;  sub(/^[0-9]+_/,"",g);                     # g = "2"
             rest=$9; sub(/^ID=[^;]*;?/,"",rest);
             $9 = "ID=" $1 "_" g (rest=="" ? "" : ";" rest);
             print
         }' "$gff" > "$gffp"

    # ---- 5.3 Comprobacion: los ID del gff deben existir en el faa ----
    grep '^>' "$faa" | sed 's/^>//' | sort -u > "$cdir/.faa_ids"
    grep -v '^#' "$gffp" | awk -F'\t' 'NF>=9{print $9}' \
        | sed 's/^ID=//; s/;.*//' | sort -u > "$cdir/.gff_ids"
    n_gff=$(wc -l < "$cdir/.gff_ids")
    n_common=$(comm -12 "$cdir/.faa_ids" "$cdir/.gff_ids" | wc -l)
    echo "   [check] $tag : IDs del gff presentes en el faa = $n_common / $n_gff"
    if [[ "$n_common" -eq 0 ]]; then
        echo "ERROR $tag :: ningun ID del gff casa con el faa, PADLOC va a fallar" | tee -a "$error_log"
    fi

    # ---- 5.4 PADLOC ----
    echo ">> $tag : PADLOC"
    if ! conda run -n "$env_padloc" padloc \
            --faa "$faa" --gff "$gffp" \
            --outdir "$cdir/padloc" \
            --cpu "$threads" 2> "$cdir/padloc.err" ; then
        echo "ERROR PADLOC en $tag :: $(tail -3 "$cdir/padloc.err" | tr '\n' ' ')" | tee -a "$error_log"
    fi

    # ---- 5.5 DefenseFinder ----
    # --db-type gembase : cada contig se trata como un replicon independiente.
    #   IMPRESCINDIBLE con un ensamblaje multi-contig. Con ordered_replicon
    #   MacSyFinder creeria que todo es UN genoma y montaria sistemas que
    #   cruzan fronteras de contig (falsos positivos). Con unordered pierdes
    #   la sintenia, que es justo lo que aporta DefenseFinder.
    echo ">> $tag : DefenseFinder"
    if ! conda run -n "$env_df" defense-finder run \
            --out-dir "$cdir/defense_finder" \
            --db-type gembase \
            -w "$threads" \
            "$faa" 2> "$cdir/defensefinder.err" ; then
        echo "ERROR DefenseFinder en $tag :: $(tail -3 "$cdir/defensefinder.err" | tr '\n' ' ')" | tee -a "$error_log"
    fi

    echo "OK $tag"
}
export -f procesar_chunk

# ---------- 6. Lanzar ----------
ls "$WORK"/chunk_*.fna | xargs -P "$PARALLEL_JOBS" -I {} \
    bash -c 'procesar_chunk "$1" "$2" "$3" "$4" "$5" "$6"' _ \
    {} "$THREADS_PER_JOB" "$OUT_DIR" "$ERROR_LOG" "$ENV_DF" "$ENV_PADLOC"

# ---------- 7. Concatenar resultados ----------
echo "[i] Concatenando salidas..."
MERGED="$OUT_DIR/merged"
mkdir -p "$MERGED"

# --- PADLOC: un *_padloc.csv por chunk (dentro de chunks/*/padloc/) ---
first=1
: > "$MERGED/padloc_all.csv"
while IFS= read -r f; do
    if [[ "$first" == "1" ]]; then cat "$f" >> "$MERGED/padloc_all.csv"; first=0
    else tail -n +2 "$f" >> "$MERGED/padloc_all.csv"; fi
done < <(find "$OUT_DIR/chunks" -path "*/padloc/*" -name "*.csv" | sort)

# --- DefenseFinder: systems / genes / hmmer ---
for tabla in systems genes hmmer; do
    first=1
    : > "$MERGED/defense_finder_${tabla}_all.tsv"
    while IFS= read -r f; do
        if [[ "$first" == "1" ]]; then cat "$f" >> "$MERGED/defense_finder_${tabla}_all.tsv"; first=0
        else tail -n +2 "$f" >> "$MERGED/defense_finder_${tabla}_all.tsv"; fi
    done < <(find "$OUT_DIR/chunks" -name "*defense_finder_${tabla}.tsv" | sort)
done

# --- PADLOC: salida HMM por dominio ---
# config.py lee F_PADLOC_DOMTBL. PADLOC escribe un domtblout por chunk
# pero no los une, asi que se hace aqui y no a mano.
: > "$MERGED/padloc_domtblout_all.txt"
find "$OUT_DIR/chunks" -path "*/padloc/*" -name "*domtblout*" \
    -exec cat {} + >> "$MERGED/padloc_domtblout_all.txt" 2>/dev/null || true

# --- faa y gff unificados (cruce con VIBRANT por contig + coordenadas) ---
cat "$OUT_DIR"/chunks/*/*.faa            > "$MERGED/vOTUs_prodigal.faa" 2>/dev/null
cat "$OUT_DIR"/chunks/*/*_padloc.gff     > "$MERGED/vOTUs_prodigal.gff" 2>/dev/null

# ---------- 8. Resumen y comprobaciones ----------
n_line() { [[ -s "$1" ]] && echo $(( $(wc -l < "$1") - 1 )) || echo 0; }

echo "=========================================================="
echo " RESUMEN"
echo "=========================================================="
echo " Contigs de entrada          : $(grep -c '^>' "$SRC")"
echo " Proteinas predichas         : $(grep -c '^>' "$MERGED/vOTUs_prodigal.faa" 2>/dev/null || echo 0)"
echo "----------------------------------------------------------"
if [[ $(n_line "$MERGED/padloc_all.csv") -gt 0 ]]; then
    echo " PADLOC - lineas             : $(n_line "$MERGED/padloc_all.csv")"
    echo " PADLOC - contigs con hit    : $(tail -n +2 "$MERGED/padloc_all.csv" | cut -d, -f1 | sort -u | wc -l)"
    echo " PADLOC - sistemas distintos : $(tail -n +2 "$MERGED/padloc_all.csv" | cut -d, -f4 | sort -u | wc -l)"
else
    echo " PADLOC - SIN SALIDA (revisa errores.log y chunks/*/padloc.err)"
fi
echo "----------------------------------------------------------"
echo " DF - sistemas COMPLETOS     : $(n_line "$MERGED/defense_finder_systems_all.tsv")   <- suele ser 0 en fagos"
echo " DF - genes en sistema       : $(n_line "$MERGED/defense_finder_genes_all.tsv")"
echo " DF - hits HMM (util)        : $(n_line "$MERGED/defense_finder_hmmer_all.tsv")   <- este es el fichero a analizar"
if [[ $(n_line "$MERGED/defense_finder_hmmer_all.tsv") -gt 0 ]]; then
    echo " DF - contigs con hit HMM    : $(tail -n +2 "$MERGED/defense_finder_hmmer_all.tsv" | awk -F'\t' '{print $2}' | sort -u | wc -l)"
    echo " DF - familias distintas     : $(tail -n +2 "$MERGED/defense_finder_hmmer_all.tsv" | awk -F'\t' '{print $5}' | sed 's/__.*//' | sort -u | wc -l)"
fi
echo "----------------------------------------------------------"
echo " Errores registrados         : $(wc -l < "$ERROR_LOG")"
echo "----------------------------------------------------------"
echo " COMPROBACION gembase:"
echo "   los replicones vistos por DefenseFinder deben acercarse al"
echo "   numero de contigs con hit. Si sale 1, gembase NO ha parseado"
echo "   las cabeceras y NO lances el run completo."
if [[ $(n_line "$MERGED/defense_finder_hmmer_all.tsv") -gt 0 ]]; then
    echo "   -> replicones: $(tail -n +2 "$MERGED/defense_finder_hmmer_all.tsv" | awk -F'\t' '{print $2}' | sort -u | wc -l)"
fi
echo "=========================================================="
echo " Resultados fusionados: $MERGED"
echo " Versiones            : $VERS"
echo " Errores              : $ERROR_LOG"
echo "=========================================================="
date
