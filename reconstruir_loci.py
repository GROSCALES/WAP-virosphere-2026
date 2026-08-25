#!/usr/bin/env python3
"""
Reconstructor de loci de defensa.
Entrada : 07_long_table_FINAL.tsv  +  mapa_v4.tsv
Salida  : 10_loci_defensa.tsv  (una fila por LOCUS, no por gen)
Regla   : agrupa hits de defensa contiguos dentro de la ventana del mapa y
          comprueba si estan los componentes que exige la literatura.
"""
import pandas as pd, numpy as np, re, sys

LONG = sys.argv[1] if len(sys.argv)>1 else "07_long_table_FINAL.tsv"
MAPA = sys.argv[2] if len(sys.argv)>2 else "data/mapa_v5.tsv"

# ============================================================
# 1. ARQUITECTURAS EN FORMATO MAQUINA
#    sintaxis:  A + B + (C|D)      ->  todos los '+' obligatorios
#               '-'                ->  sin requisito (efector/modulo autosuficiente)
# ============================================================
ARQ = {
 'Gabija'      : r'GajA + GajB',
 'Thoeris'     : r'ThsA + ThsB',
 'Hachiman'    : r'HamA + HamB',
 'Septu'       : r'PtuA + PtuB',
 'Lamassu'     : r'LmuA + LmuB',
 'Wadjet'      : r'JetA + JetB + JetC + JetD',
 'Zorya'       : r'ZorA + ZorB + (ZorC|ZorD|ZorE)',
 'Menshen'     : r'NsnA + NsnB',
 'Mza'         : r'MzaA + MzaB',
 'DISARM'      : r'DrmA + DrmB',
 'BREX'        : r'(PglX|PglZ) + (BrxA|BrxB|BrxC|BrxL)',
 'Druantia'    : r'(DruA|DruB|DruC|DruD) + (DruE|DruK|DruM)',
 'CBASS'       : r'(CD-NTase|cyclase|adenylate cyclase|E4\.6\.1\.1) + (Cap\d|Effector|E1E2|SAVED|TIR|phospholipase|transmembrane)',
 'Pycsar'      : r'(PycC|cyclase|Pycsar) + (Pyc_3TM|effector|transmembrane|TIR)',
 'Retron'      : r'(reverse transcript|RNA-directed DNA polymerase|RVT|ltrA)',
 'Ssp'         : r'SspA + SspB + SspC + SspD + (SspE|SspF|SspG|SspH)',
 'Dnd'         : r'(DndA|IscS|iscS) + DndB + DndC + DndD + DndE + (DndF|DndG|DndH|DndI)',
 'Dpd'         : r'DpdA + (folE|FolE) + (queD|QueD) + (queE|QueE) + (queC|QueC)',
 'MTasa huerfana'   : r'(methyltransferas|methylas|MTase|\bDam\b|\bdcm\b|DNMT|MT-A70) + !(endonucleas|REase|restriction enzyme|HsdR)',
 'Nucleasa putativa': r'(endonucleas|REase|restriction enzyme|nucleas)',
 'RM Type II'  : r'(methyltransferas|methylas|MTase|\bDam\b|\bdcm\b) + (endonucleas|REase|restriction enzyme)',
 'RM Type III' : r'(Mod|methyltransferas) + (Res|endonucleas|restriction enzyme)',
 'RM Type I'   : r'(HsdR|hsdR) + (HsdM|hsdM) + (HsdS|hsdS)',
 'RM Type IV'  : r'(McrB|McrC|Mrr|mcrA|mREase|RM_Type_IV)',          # NO exige MTasa
 'RM Type IIG' : r'(RM_Type_IIG|Type IIG)',                          # fusionado
 'CRISPR-Cas'  : r'(CRISPR array|crispr_array) + (cas\d|Cas\d)',     # sin array -> fragmento
 'TA/Abi'      : r'(toxin|interferase|MazF|CcdB|YoeB|RelE|HigB|Zeta|DarT) + (antitoxin|MazE|CcdA|YefM|RelB|HigA|DarG)',
}
# familias normalizadas: como mapear el 'system' del mapa a una clave de ARQ
NORM = [(r'Gabija|Gaj','Gabija'),(r'Thoeris|Ths','Thoeris'),(r'Hachiman|Ham','Hachiman'),
        (r'Septu|Ptu','Septu'),(r'Lamassu|Lmu','Lamassu'),(r'Wadjet|Jet','Wadjet'),
        (r'Zorya|Zor','Zorya'),(r'Menshen|Nsn','Menshen'),(r'Mza','Mza'),
        (r'DISARM|Drm','DISARM'),(r'BREX|Brx|Pgl','BREX'),(r'Druantia|Dru','Druantia'),
        (r'CBASS','CBASS'),(r'Pycsar|Pyc','Pycsar'),(r'Retron|retron|DRT','Retron'),
        (r'DNA methylation|anti-restriction|MTase_II|oMTase','MTasa huerfana'),(r'Ssp','Ssp'),(r'Dnd','Dnd'),(r'Dpd','Dpd'),
        (r'RM Type IIG','RM Type IIG'),(r'RM Type IV','RM Type IV'),
        (r'RM Type III','RM Type III'),(r'RM Type II|Restriction-Modification|^RM$','RM Type II'),
        (r'HsdR','RM Type I'),(r'CRISPR|Cas','CRISPR-Cas'),
        (r'Toxin-Antitoxin|TA|Abi|MazEF|HicAB|DarTG','TA/Abi')]

def familia(system, cat_raw):
    t = f"{system} {cat_raw}"
    for pat,f in NORM:
        if re.search(pat, t): return f
    return system

def evalua(expr, texto, textos=None):
    """devuelve (cumple, n_encontrados, n_requeridos, faltan[])

    REGLA DE GENES DISTINTOS. Cada termino POSITIVO lo tiene que aportar un
    gen distinto del bloque. Sin esto, una etiqueta como
    'Toxin_Antitoxin_and_Abi' satisface a la vez (toxin|...) y
    (antitoxin|...) ella sola, y el locus sale completo aunque no exista
    ninguna pareja toxina-antitoxina. Afectaba a 104 de 112 loci TA/Abi y a
    los 8 de Druantia ('Druantia_II' contiene "DruA").

    Las NEGACIONES se evaluan sobre el texto COMPLETO del bloque: la
    ausencia tiene que darse en todo el vecindario, no en un gen suelto.
    Por eso 'MTasa huerfana' no se ve afectada, y es correcto que no lo
    este.

    textos = lista con el texto de cada gen. Si no se pasa, se comporta
    como antes (compatibilidad).
    """
    if expr in ('-','', None): return True, 1, 1, []
    if textos is None: textos = [texto]
    tokens = [t.strip() for t in expr.split('+')]
    pos, neg = [], []
    for tk in tokens:
        (neg if tk.startswith('!') else pos).append(tk)

    ok, faltan = 0, []

    # --- negaciones: sobre todo el texto del bloque ---
    for tk in neg:
        pat = tk[1:]
        pat = pat[1:-1] if pat.startswith('(') and pat.endswith(')') else pat
        if not re.search(pat, texto, re.I): ok += 1
        else: faltan.append('NO debe estar: ' + tk)

    # --- positivos: cada uno de un gen distinto (emparejamiento de Kuhn) ---
    cands = []
    for tk in pos:
        pat = tk[1:-1] if tk.startswith('(') and tk.endswith(')') else tk
        rx = re.compile(pat, re.I)
        cands.append([i for i, t in enumerate(textos) if rx.search(t)])

    emparejado = {}          # gen -> indice de termino
    def intenta(i, visto):
        for g in cands[i]:
            if g in visto: continue
            visto.add(g)
            if g not in emparejado or intenta(emparejado[g], visto):
                emparejado[g] = i
                return True
        return False

    for i in sorted(range(len(pos)), key=lambda i: len(cands[i])):
        intenta(i, set())
    asignados = set(emparejado.values())
    ok += len(asignados)
    for i, tk in enumerate(pos):
        if i not in asignados: faltan.append(tk)

    return ok == len(tokens), ok, len(tokens), faltan

# ============================================================
# 2. CARGA
# ============================================================
l = pd.read_csv(LONG, sep='\t', low_memory=False, keep_default_na=False)
m = pd.read_csv(MAPA, sep='\t', keep_default_na=False)
f = l[l.Passes_catalog=="Yes"].copy() if 'Passes_catalog' in l else l.copy()
d = f[f.gene_class=="Defence"].copy()

mm = m.set_index('category_raw')
for col in ['eje','subcategoria','system','nivel_entidad','tipo_de_llamada_original',
            'ventana_maxima','confianza_si_aislado','confianza_si_completo','metrica']:
    if col in m.columns:
        d[col] = d.category.map(mm[col]) if col!='system' else d.category.map(mm['system'])
d = d[~d.eje.isna()].copy()
d['familia'] = [familia(s,c) for s,c in zip(d.system.astype(str), d.category.astype(str))]
d['txt'] = (d.ko_name.astype(str)+' '+d.pfam_name.astype(str)+' '+d.df_system.astype(str)+' '+
            d.pl_system.astype(str)+' '+d.pl_profile.astype(str)+' '+d.category.astype(str))
d['start_n'] = pd.to_numeric(d.start, errors='coerce')
d['orf_n']   = pd.to_numeric(d.orf_num, errors='coerce')
d['len_n']   = pd.to_numeric(d.contig_length, errors='coerce')

def ventana(v):
    mg = re.search(r'(\d+)\s*genes', str(v)); mk = re.search(r'(\d+)\s*kb', str(v))
    return (int(mg.group(1)) if mg else 10, (int(mk.group(1)) if mk else 20)*1000)

# ============================================================
# 3. AGRUPAR EN LOCI
# ============================================================
loci, lid = [], 0
for ctg, sub in d.groupby('votu_id'):
    sub = sub.sort_values('start_n').reset_index(drop=True)
    grupo, actual = [], []
    for i, r in sub.iterrows():
        if not actual: actual = [i]; continue
        prev = sub.loc[actual[-1]]
        mg, mk = ventana(r.ventana_maxima)
        dg = abs(r.orf_n - prev.orf_n) if pd.notna(r.orf_n) and pd.notna(prev.orf_n) else 99
        db = abs(r.start_n - prev.start_n) if pd.notna(r.start_n) and pd.notna(prev.start_n) else 9e9
        if dg <= mg and db <= mk: actual.append(i)
        else: grupo.append(actual); actual = [i]
    if actual: grupo.append(actual)

    for g in grupo:
        blk = sub.loc[g]
        texto = ' '.join(blk.txt)
        textos = list(blk.txt)
        for fam, fsub in blk.groupby('familia'):
            if str(fsub.eje.mode().iat[0]).startswith('X'): continue
            lid += 1
            expr = ARQ.get(fam, '-')
            cumple, nok, nreq, faltan = evalua(expr, texto, textos)
            # un perfil de COMPONENTE sin arquitectura definida no se auto-completa
            nivel = str(fsub.nivel_entidad.mode().iat[0]) if 'nivel_entidad' in fsub else ''
            if expr == '-' and nivel == 'componente':
                cumple, nok, nreq = False, 0, 1
                faltan = ['arquitectura de sistema no reconstruida']
            eje  = fsub.eje.mode().iat[0]
            borde = (blk.start_n.min() < 2000) or (pd.notna(blk.len_n.iat[0]) and blk.end.astype(float).max() > blk.len_n.iat[0]-2000)
            # escision RM: un RM incompleto se resuelve segun que mitad tenga
            reasignado = ''
            if not cumple and fam.startswith('RM Type'):
                tiene_mt  = re.search(r'methyltransferas|methylas|MTase|\bDam\b|\bdcm\b|DNMT|MT-A70', texto, re.I)
                tiene_nuc = re.search(r'endonucleas|REase|restriction enzyme|HsdR', texto, re.I)
                if tiene_mt and not tiene_nuc:
                    fam, eje_new, reasignado = 'MTasa huerfana', '1 Contra-defensa', 'RM incompleto -> MTasa huerfana'
                    expr = ARQ['MTasa huerfana']; cumple, nok, nreq, faltan = evalua(expr, texto, textos)
                elif tiene_nuc and not tiene_mt:
                    fam, eje_new, reasignado = 'Nucleasa putativa', '2 Defensa', 'RM incompleto -> nucleasa putativa'
                    expr = ARQ['Nucleasa putativa']; cumple, nok, nreq, faltan = evalua(expr, texto, textos)
                else: eje_new = None
                if reasignado: eje = eje_new
            if cumple:
                comp, conf = 'completo', fsub.confianza_si_completo.mode().iat[0]
                if reasignado: conf = 'candidata' if 'MTasa' in fam else 'media'
            elif borde:
                comp, conf = 'fragmento (borde de contig)', 'media'
            else:
                comp, conf = 'fragmento', fsub.confianza_si_aislado.mode().iat[0]
            loci.append(dict(
                locus_id=f'L{lid:05d}', contig_id=ctg, familia=fam, eje=eje,
                subcategoria=fsub.subcategoria.mode().iat[0],
                start=int(blk.start_n.min()) if pd.notna(blk.start_n.min()) else -1,
                end=int(pd.to_numeric(blk.end,errors='coerce').max()),
                n_genes_locus=len(fsub), n_genes_bloque=len(blk),
                componentes_requeridos=expr, n_componentes_observados=nok,
                n_componentes_requeridos=nreq, faltan='; '.join(faltan),
                completitud_locus=comp, confianza_locus=conf,
                borde_de_contig=bool(borde),
                metrica='confirmada' if (cumple and conf=='alta') else 'candidata',
                reasignacion=reasignado,
                genes='; '.join(sorted(set(fsub.category)))))

L = pd.DataFrame(loci)
L.to_csv('10_loci_defensa.tsv', sep='\t', index=False)

# ============================================================
# 4. INFORME
# ============================================================
print('='*84); print(' RECONSTRUCCION DE LOCI')
print(f'  genes de defensa de entrada : {len(d)}')
print(f'  LOCI reconstruidos          : {len(L)}   en {L.contig_id.nunique()} contigs')
print(f'  genes agrupados en loci     : {L.n_genes_locus.sum()}  (excluidos aparte)')
print(f'  loci con >1 gen (colapso)   : {(L.n_genes_locus>1).sum()}')
print(); print(' completitud:'); print(L.completitud_locus.value_counts().to_string())
print(); print(' LAS DOS CIFRAS QUE NO SE SUMAN:')
print(L.metrica.value_counts().to_string())
print(); print(' por eje:')
print(pd.crosstab(L.eje, L.completitud_locus).to_string())
print(); print(' familias multigenicas: cuantas reconstruyen?')
t = L[L.n_componentes_requeridos>1].groupby('familia').agg(
    loci=('locus_id','size'), completos=('completitud_locus', lambda s:(s=='completo').sum()))
t['%'] = (100*t.completos/t.loci).round(0)
print(t.sort_values('loci',ascending=False).to_string())
print(f"\n -> 10_loci_defensa.tsv")
