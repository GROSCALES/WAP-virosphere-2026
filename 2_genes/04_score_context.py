#!/usr/bin/env python3
"""
04_score_context.py — genomic context, DRAM-v style, applied to EVERY ORF.

The old pipeline required >=1 viral-like neighbour anywhere within +/-10
ORFs. In DRAM-v terms that is auxiliary score 4 (signal on one flank only),
i.e. LOW confidence — and it was applied only to the regex-derived defence
hits, never to the genes reclassified from VIBRANT's AMG list. Two evidence
standards inside the same category.

Here one rule is computed for every ORF in the table, so AVGs and defence
genes are held to exactly the same criterion.

DRAM-v auxiliary score (Shaffer et al. 2020, NAR 48:8883; the default of
`DRAM-v distill` is Max auxiliary score: 3):

    1  hallmark gene on BOTH flanks                      high confidence
    2  hallmark on one flank, viral-like on the other    high confidence
    3  viral-like on BOTH flanks                         medium confidence
    4  signal on ONE flank only                          low
    5  no signal on either flank, or gene at contig end  discard

Definitions used here, mirroring DRAM-v's two tiers:

    hallmark   : a structural/virion keyword in ANY name field
                 (capsid, terminase, portal, tail, baseplate, holin, ...)
    viral-like : the ORF has a VOG assignment at all

Note on which field is searched: for CONTEXT we deliberately do look at
'VOG name'. That field is the description of a viral reference protein,
which is exactly the wrong thing for deciding a gene's own identity (it
caused every false-friend collision in the old pipeline) but exactly the
right thing for asking "does this neighbourhood look viral?".

Adds: n_hall_left, n_viral_left, n_hall_right, n_viral_right,
      aux_score, at_contig_end, passes_context

Input : 02_gene_table_L3_evidence.tsv
Output: 03_gene_table_context.tsv
"""

import os
from collections import defaultdict

from config import (
    F_GENE_TABLE_L3, F_GENE_TABLE_CTX,
    CTX_WINDOW_ORFS, CTX_VIRAL_KEYWORDS, require,
)

KEYWORDS = [k.lower() for k in CTX_VIRAL_KEYWORDS]


def is_hallmark(ko_name, pfam_name, vog_name):
    blob = f"{ko_name} {pfam_name} {vog_name}".lower()
    return any(k in blob for k in KEYWORDS)


def aux_score(hl, vl, hr, vr, at_end):
    """Return the DRAM-v auxiliary score for one ORF."""
    if at_end:
        return 5
    left_h, right_h = hl > 0, hr > 0
    left_v, right_v = (hl + vl) > 0, (hr + vr) > 0
    if left_h and right_h:
        return 1
    if (left_h and right_v) or (right_h and left_v):
        return 2
    if left_v and right_v:
        return 3
    if left_v or right_v:
        return 4
    return 5


def main():
    require(F_GENE_TABLE_L3)
    print("=" * 66)
    print(" GENOMIC CONTEXT — DRAM-v auxiliary score, every ORF")
    print("=" * 66)

    with open(F_GENE_TABLE_L3) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    for c in ("protein_id", "votu_id", "orf_num", "ko_name", "pfam_name",
              "vog_name", "vog"):
        if c not in idx:
            raise SystemExit(f"[SCHEMA ERROR] column '{c}' missing from L3")

    # ---- pass 1: per contig, record each ORF's hallmark / viral-like state
    print("[1] scanning ORFs ...")
    contigs = defaultdict(dict)          # votu -> {orf_num: (hall, viral)}
    n = 0
    with open(F_GENE_TABLE_L3) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            n += 1
            votu = c[idx["votu_id"]]
            orf = int(c[idx["orf_num"]])
            hall = is_hallmark(c[idx["ko_name"]], c[idx["pfam_name"]],
                               c[idx["vog_name"]])
            viral = bool(c[idx["vog"]].strip())
            contigs[votu][orf] = (hall, viral)
    print(f"    ORFs: {n}   contigs: {len(contigs)}")

    # ---- pass 2: score, streaming so memory stays flat --------------------
    print("[2] scoring and writing ...")
    new_cols = ["n_hall_left", "n_viral_left", "n_hall_right", "n_viral_right",
                "aux_score", "at_contig_end", "passes_context"]
    dist = defaultdict(int)
    n_pass = 0

    with open(F_GENE_TABLE_L3) as fin, open(F_GENE_TABLE_CTX, "w") as out:
        next(fin)
        out.write("\t".join(header + new_cols) + "\n")
        for line in fin:
            c = line.rstrip("\n").split("\t")
            votu = c[idx["votu_id"]]
            orf = int(c[idx["orf_num"]])
            neigh = contigs[votu]
            last = max(neigh)

            hl = vl = hr = vr = 0
            for d in range(1, CTX_WINDOW_ORFS + 1):
                p = neigh.get(orf - d)
                if p:
                    hl += p[0]
                    vl += p[1] and not p[0]
                q = neigh.get(orf + d)
                if q:
                    hr += q[0]
                    vr += q[1] and not q[0]

            # "gene at the end of a contig": no ORF on one side at all
            at_end = (orf == 1) or (orf == last)
            score = aux_score(hl, vl, hr, vr, at_end)
            ok = score <= 3
            dist[score] += 1
            n_pass += ok

            c += [str(hl), str(vl), str(hr), str(vr), str(score),
                  "1" if at_end else "0", "1" if ok else "0"]
            out.write("\t".join(c) + "\n")

    # ---- report -----------------------------------------------------------
    print("\n" + "=" * 66)
    print(" AUXILIARY SCORE DISTRIBUTION")
    print("=" * 66)
    labels = {1: "hallmark both flanks      (high)",
              2: "hallmark + viral-like     (high)",
              3: "viral-like both flanks    (medium)",
              4: "signal one flank only     (low)",
              5: "no signal / contig end    (discard)"}
    for s in (1, 2, 3, 4, 5):
        print(f"  {s}  {labels[s]:<34} {dist[s]:>8}  ({100.0*dist[s]/n:.1f}%)")
    print("-" * 66)
    print(f"  passes_context (score <= 3)  : {n_pass}  ({100.0*n_pass/n:.1f}%)")
    print("=" * 66)
    print(" This is a FLAG, not a filter. Nothing was removed. Scripts")
    print(" downstream decide what to do with it, and the sensitivity")
    print(" analysis is a filter() in R, not another run of this pipeline.")
    print(f"\n Output: {F_GENE_TABLE_CTX}")


if __name__ == "__main__":
    main()
