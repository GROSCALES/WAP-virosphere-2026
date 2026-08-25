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
#  vOTU Defence System Screening (PADLOC + DefenseFinder)
#
#  INPUT & PARALLELIZATION
#   Input: A single multi-contig FASTA. Parallelization is achieved by 
#   chunking the FASTA file to optimize cluster resource usage.
#
#  UNIFIED GENE PREDICTION
#   `prodigal -p meta` is executed independently once per chunk. This ensures 
#   that both DefenseFinder and PADLOC operate on the exact same set of ORFs 
#   and correctly handles short viral contigs.
#
#   GFF Harmonization: The awk block rewrites the GFF ID field using the 
#   true contig name (Column 1) to guarantee the strict identifier matching 
#   required by PADLOC's parser.
#
#  DOWNSTREAM INTEGRATION
#   Downstream merging with VIBRANT annotations must be done via 
#   CONTIG + COORDINATES, never by protein ID string, to correctly map 
#   excised prophages which possess an independent coordinate frame.
#
#  DEFENSEFINDER OUTPUT
#   The pipeline prioritizes 'hmmer.tsv' (raw HMM profile hits) over 
#   'systems.tsv', as viral genomes frequently carry isolated or orphan 
#   defence components rather than complete, canonical systems.
#
#  PILOT MODE
#   Usage: PILOT=1 sbatch 08_defensas_votus.sh
#   Allows validation of gembase parsing (ensuring the 'replicons' count 
#   matches the expected contig hit count) before launching the full job.
# ==========================================================================

# ---------- 1. Paths ----------
WAP_BASE="${WAP_BASE:-$HOME/gabi/vOTUs_total}"
IN_FASTA="${IN_FASTA:-$WAP_BASE/vOTUs_catalogue_rep_seq.fasta}"
OUT_DIR="${OUT_DIR:-$WAP_BASE/analisis_defensas_votus}"

# Conda environments
ENV_DF="defensefinder_env"
ENV_PADLOC="padloc"

# Prodigal binary absolute path (bypasses 'conda run' overhead)
PRODIGAL_BIN="${PRODIGAL_BIN:-$HOME/miniconda3/envs/padloc/bin/prodigal}"
export PRODIGAL_BIN

PILOT="${PILOT:-0}"
PILOT_N=200

mkdir -p "$OUT_DIR"
WORK="$OUT_DIR/work"
LOGS="$OUT_DIR/logs"
mkdir -p "$WORK" "$LOGS"
ERROR_LOG="$OUT_DIR/errors.log"
: > "$ERROR_LOG"

if [[ ! -s "$IN_FASTA" ]]; then
    echo "ERROR: Cannot find input FASTA: $IN_FASTA" | tee -a "$ERROR_LOG"
    exit 1
fi
if [[ ! -x "$PRODIGAL_BIN" ]]; then
    echo "ERROR: Cannot find Prodigal binary at $PRODIGAL_BIN" | tee -a "$ERROR_LOG"
    exit 1
fi

# ---------- 2. Core Allocation ----------
CORES_TOTAL="${SLURM_CPUS_PER_TASK:-30}"
PARALLEL_JOBS=6
THREADS_PER_JOB=$(( CORES_TOTAL / PARALLEL_JOBS ))
[[ "$THREADS_PER_JOB" -lt 1 ]] && THREADS_PER_JOB=1

echo "=========================================================="
echo " vOTU DEFENCE SCREENING"
echo " Input      : $IN_FASTA"
echo " Output     : $OUT_DIR"
echo " Cores      : $CORES_TOTAL | Parallel chunks: $PARALLEL_JOBS | Threads/chunk: $THREADS_PER_JOB"
[[ "$PILOT" == "1" ]] && echo " PILOT MODE : Processing first $PILOT_N contigs only"
echo "=========================================================="

# ---------- 3. Version Tracking (Required for Methods) ----------
VERS="$OUT_DIR/VERSIONS.txt"
{
    echo "Date : $(date -Iseconds)"
    echo "Host : $(hostname)"
    echo "Input: $IN_FASTA"
    echo "MD5  : $(md5sum "$IN_FASTA" | cut -d' ' -f1)"
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
echo "[i] Tool versions documented in $VERS"

# ---------- 4. Chunking ----------
N_CONTIGS=$(grep -c '^>' "$IN_FASTA")
echo "[i] Total contigs in catalogue: $N_CONTIGS"

SRC="$IN_FASTA"
if [[ "$PILOT" == "1" ]]; then
    SRC="$WORK/pilot.fna"
    awk -v n="$PILOT_N" '/^>/{c++} c>n{exit} {print}' "$IN_FASTA" > "$SRC"
    echo "[i] Pilot mode active: $(grep -c '^>' "$SRC") contigs extracted"
fi

rm -f "$WORK"/chunk_*.fna
awk -v n="$PARALLEL_JOBS" -v pre="$WORK/chunk" '
    /^>/ { i++; c = ((i-1) % n) + 1 }
    { print > (pre "_" c ".fna") }
' "$SRC"
echo "[i] Catalogue split into $(ls "$WORK"/chunk_*.fna 2>/dev/null | wc -l) chunks"

# ---------- 5. Single Chunk Processing Function ----------
procesar_chunk() {
    local fna="$1" threads="$2" out_dir="$3" error_log="$4" env_df="$5" env_padloc="$6"
    local tag faa gff gffp cdir n_gff n_common

    tag=$(basename "$fna" .fna)
    cdir="$out_dir/chunks/$tag"
    mkdir -p "$cdir/padloc" "$cdir/defense_finder"

    faa="$cdir/$tag.faa"
    gff="$cdir/$tag.gff"          # Raw prodigal GFF (ID=1_2)
    gffp="$cdir/${tag}_padloc.gff" # Harmonized GFF for PADLOC (ID=<contig>_2)

    # ---- 5.1 Gene Prediction ----
    echo ">> $tag : Gene prediction (prodigal -p meta)"
    if ! "$PRODIGAL_BIN" -i "$fna" -a "$faa" -f gff -o "$gff" -p meta -q 2> "$cdir/prodigal.err" ; then
        echo "ERROR prodigal on $tag :: $(tail -3 "$cdir/prodigal.err" | tr '\n' ' ')" | tee -a "$error_log"
        return
    fi
    
    # Clean FAA headers: remove ' # start # end # strand # attrs'
    # Leaves strictly <contig>_<n>, matching gembase/MacSyFinder requirements
    sed -i 's/ #.*//' "$faa"

    # ---- 5.2 Rewrite GFF ID= to match FAA headers ----
    awk 'BEGIN{FS=OFS="\t"}
         /^#/{print; next}
         NF>=9{
             n=$9; sub(/^.*ID=/,"",n); sub(/;.*/,"",n);      # n = "1_2"
             g=n;  sub(/^[0-9]+_/,"",g);                     # g = "2"
             rest=$9; sub(/^ID=[^;]*;?/,"",rest);
             $9 = "ID=" $1 "_" g (rest=="" ? "" : ";" rest);
             print
         }' "$gff" > "$gffp"

    # ---- 5.3 Verification: GFF IDs must exactly match FAA IDs ----
    grep '^>' "$faa" | sed 's/^>//' | sort -u > "$cdir/.faa_ids"
    grep -v '^#' "$gffp" | awk -F'\t' 'NF>=9{print $9}' \
        | sed 's/^ID=//; s/;.*//' | sort -u > "$cdir/.gff_ids"
    n_gff=$(wc -l < "$cdir/.gff_ids")
    n_common=$(comm -12 "$cdir/.faa_ids" "$cdir/.gff_ids" | wc -l)
    echo "   [check] $tag : GFF IDs present in FAA = $n_common / $n_gff"
    if [[ "$n_common" -eq 0 ]]; then
        echo "ERROR $tag :: No GFF IDs match FAA headers; PADLOC will fail" | tee -a "$error_log"
    fi

    # ---- 5.4 PADLOC ----
    echo ">> $tag : PADLOC"
    if ! conda run -n "$env_padloc" padloc \
            --faa "$faa" --gff "$gffp" \
            --outdir "$cdir/padloc" \
            --cpu "$threads" 2> "$cdir/padloc.err" ; then
        echo "ERROR PADLOC on $tag :: $(tail -3 "$cdir/padloc.err" | tr '\n' ' ')" | tee -a "$error_log"
    fi

    # ---- 5.5 DefenseFinder ----
    # --db-type gembase: Treats each contig as an independent replicon.
    # CRITICAL for multi-contig assemblies to prevent cross-contig false 
    # positive systems, while preserving synteny (unlike 'unordered' mode).
    echo ">> $tag : DefenseFinder"
    if ! conda run -n "$env_df" defense-finder run \
            --out-dir "$cdir/defense_finder" \
            --db-type gembase \
            -w "$threads" \
            "$faa" 2> "$cdir/defensefinder.err" ; then
        echo "ERROR DefenseFinder on $tag :: $(tail -3 "$cdir/defensefinder.err" | tr '\n' ' ')" | tee -a "$error_log"
    fi

    echo "OK $tag"
}
export -f procesar_chunk

# ---------- 6. Execution ----------
ls "$WORK"/chunk_*.fna | xargs -P "$PARALLEL_JOBS" -I {} \
    bash -c 'procesar_chunk "$1" "$2" "$3" "$4" "$5" "$6"' _ \
    {} "$THREADS_PER_JOB" "$OUT_DIR" "$ERROR_LOG" "$ENV_DF" "$ENV_PADLOC"

# ---------- 7. Concatenate Results ----------
echo "[i] Merging outputs..."
MERGED="$OUT_DIR/merged"
mkdir -p "$MERGED"

# --- PADLOC: One *_padloc.csv per chunk ---
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

# --- PADLOC: HMM domain output ---
# config.py relies on F_PADLOC_DOMTBL. Concatenating domtblout securely.
: > "$MERGED/padloc_domtblout_all.txt"
find "$OUT_DIR/chunks" -path "*/padloc/*" -name "*domtblout*" \
    -exec cat {} + >> "$MERGED/padloc_domtblout_all.txt" 2>/dev/null || true

# --- Unified FAA and GFF (For VIBRANT mapping by contig+coords) ---
cat "$OUT_DIR"/chunks/*/*.faa            > "$MERGED/vOTUs_prodigal.faa" 2>/dev/null
cat "$OUT_DIR"/chunks/*/*_padloc.gff     > "$MERGED/vOTUs_prodigal.gff" 2>/dev/null

# ---------- 8. Summary and Validation ----------
n_line() { [[ -s "$1" ]] && echo $(( $(wc -l < "$1") - 1 )) || echo 0; }

echo "=========================================================="
echo " SUMMARY"
echo "=========================================================="
echo " Input contigs             : $(grep -c '^>' "$SRC")"
echo " Predicted proteins        : $(grep -c '^>' "$MERGED/vOTUs_prodigal.faa" 2>/dev/null || echo 0)"
echo "----------------------------------------------------------"
if [[ $(n_line "$MERGED/padloc_all.csv") -gt 0 ]]; then
    echo " PADLOC - Total lines      : $(n_line "$MERGED/padloc_all.csv")"
    echo " PADLOC - Contigs w/ hits  : $(tail -n +2 "$MERGED/padloc_all.csv" | cut -d, -f1 | sort -u | wc -l)"
    echo " PADLOC - Distinct systems : $(tail -n +2 "$MERGED/padloc_all.csv" | cut -d, -f4 | sort -u | wc -l)"
else
    echo " PADLOC - NO OUTPUT (Check errors.log and chunks/*/padloc.err)"
fi
echo "----------------------------------------------------------"
echo " DF - COMPLETE systems     : $(n_line "$MERGED/defense_finder_systems_all.tsv")  (Often 0 in viral genomes)"
echo " DF - Genes in systems     : $(n_line "$MERGED/defense_finder_genes_all.tsv")"
echo " DF - HMM hits (core file) : $(n_line "$MERGED/defense_finder_hmmer_all.tsv")"
if [[ $(n_line "$MERGED/defense_finder_hmmer_all.tsv") -gt 0 ]]; then
    echo " DF - Contigs w/ HMM hits  : $(tail -n +2 "$MERGED/defense_finder_hmmer_all.tsv" | awk -F'\t' '{print $2}' | sort -u | wc -l)"
    echo " DF - Distinct families    : $(tail -n +2 "$MERGED/defense_finder_hmmer_all.tsv" | awk -F'\t' '{print $5}' | sed 's/__.*//' | sort -u | wc -l)"
fi
echo "----------------------------------------------------------"
echo " Registered Errors         : $(wc -l < "$ERROR_LOG")"
echo "----------------------------------------------------------"
echo " GEMBASE PARSING VALIDATION:"
echo "   The number of replicons detected by DefenseFinder must align closely"
echo "   with the number of contigs with hits. If it reads 1, gembase failed"
echo "   to parse the headers properly."
if [[ $(n_line "$MERGED/defense_finder_hmmer_all.tsv") -gt 0 ]]; then
    echo "   -> Replicons parsed     : $(tail -n +2 "$MERGED/defense_finder_hmmer_all.tsv" | awk -F'\t' '{print $2}' | sort -u | wc -l)"
fi
echo "=========================================================="
echo " Merged Results Directory : $MERGED"
echo " Version log              : $VERS"
echo " Error log                : $ERROR_LOG"
echo "=========================================================="
date
