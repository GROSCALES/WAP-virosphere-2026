#!/usr/bin/env python3
"""
genes/02_add_evidence.py — LEVEL 3: Integrates specialist tool evidence into the gene table.

Appends raw output from DefenseFinder and PADLOC as independent columns:
  DefenseFinder : df_system, df_profile, df_ievalue, df_profile_cov,
                  df_seq_cov, df_in_complete_system
  PADLOC        : pl_system, pl_profile, pl_protein, pl_evalue,
                  pl_hmm_cov, pl_target_cov, pl_in_complete_system

No genes are classified or filtered at this stage. Tool outputs are preserved 
independently; final classification is resolved downstream (in step 05) via 
explicit precedence rules.

PADLOC Thresholds:
The raw hmmsearch output is filtered using PADLOC's own per-profile thresholds 
(from data/hmm_meta.txt): e.val.threshold, hmm.coverage.threshold, and 
target.coverage.threshold. 
Coverage is summed over all domains of a (protein, profile) pair. This strictly 
reproduces PADLOC's gene-level filtering prior to any system-level architectural 
evaluation.

Input : 01_gene_table_L2.tsv
Output: 02_gene_table_L3_evidence.tsv
"""
import os
import re
from collections import defaultdict

from config import (
    F_GENE_TABLE_L2, F_GENE_TABLE_L3,
    F_DF_HMMER, F_DF_GENES, F_PADLOC_CSV, F_PADLOC_DOMTBL, F_PADLOC_HMM_META,
    DF_HMMER_COLS, check_header, require,
)

TRAILING_NUM = re.compile(r"_\d+$")


def family_from_profile(name):
    """DefenseFinder: 'RM_Type_II__Type_II_MTases_FAM_22' -> 'RM_Type_II'
       PADLOC      : 'MTase_II_00164' -> 'MTase_II' ; 'DinG' -> 'DinG'"""
    if "__" in name:
        return name.split("__", 1)[0]
    return TRAILING_NUM.sub("", name)


def to_float(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ==========================================================================
# DefenseFinder
# ==========================================================================
def load_defensefinder(f_hmmer, f_genes, valid_ids):
    print(f"[1] DefenseFinder — {os.path.basename(f_hmmer)}")
    check_header(f_hmmer, DF_HMMER_COLS)

    best = {}
    n_rows = 0
    with open(f_hmmer) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            n_rows += 1
            pid, profile = c[0], c[4]
            ev = to_float(c[5])
            if ev is None:
                continue
            prev = best.get(pid)
            if prev is None or ev < prev["df_ievalue_f"]:
                best[pid] = {
                    "df_system": family_from_profile(profile),
                    "df_profile": profile,
                    "df_ievalue": c[5],
                    "df_ievalue_f": ev,
                    "df_profile_cov": c[7],
                    "df_seq_cov": c[8],
                }
    print(f"    hit rows                 : {n_rows}")
    print(f"    proteins with a hit      : {len(best)}")

    # --- proteins belonging to a COMPLETE system -------------------------
    in_system = set()
    if os.path.exists(f_genes) and os.path.getsize(f_genes) > 0:
        with open(f_genes) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows = [l.rstrip("\n").split("\t") for l in fh]
        # find the column that actually holds protein IDs (data-driven, so a
        # change in DefenseFinder's output shape cannot silently break this)
        best_col, best_hits = None, 0
        for i in range(len(header)):
            vals = {r[i] for r in rows if len(r) > i}
            hits = len(vals & valid_ids)
            if hits > best_hits:
                best_col, best_hits = i, hits
        if best_col is not None and best_hits:
            in_system = {r[best_col] for r in rows if len(r) > best_col} & valid_ids
            print(f"    complete-system genes    : {len(in_system)} "
                  f"(column '{header[best_col]}')")
        else:
            print("    [warn] could not locate a protein-ID column in "
                  f"{os.path.basename(f_genes)}")
    else:
        print("    complete-system genes    : 0 (file empty)")

    return best, in_system


# ==========================================================================
# PADLOC
# ==========================================================================
def load_padloc_thresholds(path):
    print(f"[2] PADLOC thresholds — {os.path.basename(path)}")
    thr = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(h.strip() for h in header)}
        need = ["hmm.name", "protein.name", "e.val.threshold",
                "hmm.coverage.threshold", "target.coverage.threshold"]
        missing = [n for n in need if n not in idx]
        if missing:
            raise SystemExit(f"[SCHEMA ERROR] hmm_meta.txt lacks {missing}")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) <= idx["target.coverage.threshold"]:
                continue
            thr[c[idx["hmm.name"]]] = {
                "protein": c[idx["protein.name"]],
                "e": to_float(c[idx["e.val.threshold"]], 1e-5),
                "hcov": to_float(c[idx["hmm.coverage.threshold"]], 0.0),
                "tcov": to_float(c[idx["target.coverage.threshold"]], 0.0),
            }
    print(f"    profiles with thresholds : {len(thr)}")
    return thr


def load_padloc(f_domtbl, thr, f_csv, valid_ids):
    print(f"[3] PADLOC hits — {os.path.basename(f_domtbl)}")

    # accumulate domain coverage per (protein, profile)
    acc = {}
    n_rows = 0
    with open(f_domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split()
            if len(c) < 20:
                continue
            n_rows += 1
            target, tlen = c[0], to_float(c[2])
            query, qlen = c[3], to_float(c[5])
            ev = to_float(c[6])
            hfrom, hto = to_float(c[15]), to_float(c[16])
            afrom, ato = to_float(c[17]), to_float(c[18])
            if None in (tlen, qlen, ev, hfrom, hto, afrom, ato) or not (tlen and qlen):
                continue
            k = (target, query)
            d = acc.get(k)
            if d is None:
                d = acc[k] = {"tlen": tlen, "qlen": qlen, "e": ev,
                              "h": 0.0, "a": 0.0}
            d["h"] += hto - hfrom + 1
            d["a"] += ato - afrom + 1

    print(f"    domain rows              : {n_rows}")
    print(f"    (protein, profile) pairs : {len(acc)}")

    # apply each profile's own thresholds
    best = {}
    n_pass = 0
    no_thr = set()
    for (target, query), d in acc.items():
        t = thr.get(query)
        if t is None:
            no_thr.add(query)
            continue
        hcov = d["h"] / d["qlen"]
        tcov = d["a"] / d["tlen"]
        if d["e"] > t["e"] or hcov < t["hcov"] or tcov < t["tcov"]:
            continue
        n_pass += 1
        prev = best.get(target)
        if prev is None or d["e"] < prev["pl_evalue_f"]:
            best[target] = {
                "pl_system": family_from_profile(query),
                "pl_profile": query,
                "pl_protein": t["protein"],
                "pl_evalue": f"{d['e']:.2e}",
                "pl_evalue_f": d["e"],
                "pl_hmm_cov": f"{hcov:.3f}",
                "pl_target_cov": f"{tcov:.3f}",
            }

    print(f"    pairs passing thresholds : {n_pass}")
    print(f"    proteins retained        : {len(best)}")
    if no_thr:
        print(f"    [warn] {len(no_thr)} profiles absent from hmm_meta.txt "
              f"(e.g. {sorted(no_thr)[:3]})")

    # --- proteins in a COMPLETE system (PADLOC's own CSV) -----------------
    in_system = set()
    if os.path.exists(f_csv) and os.path.getsize(f_csv) > 0:
        with open(f_csv) as fh:
            header = fh.readline().rstrip("\n").split(",")
            idx = {n.strip(): i for i, n in enumerate(header)}
            if "target.name" in idx:
                j = idx["target.name"]
                for line in fh:
                    c = line.rstrip("\n").split(",")
                    if len(c) > j:
                        in_system.add(c[j])
                in_system &= valid_ids
            else:
                print("    [warn] 'target.name' not found in padloc_all.csv")
    print(f"    complete-system genes    : {len(in_system)}")

    return best, in_system


# ==========================================================================
# MAIN — stream L2, append columns, write L3
# ==========================================================================
def main():
    require(F_GENE_TABLE_L2, F_DF_HMMER, F_PADLOC_DOMTBL, F_PADLOC_HMM_META)

    print("=" * 66)
    print(" LEVEL 3 — DefenseFinder + PADLOC as evidence columns")
    print("=" * 66)

    with open(F_GENE_TABLE_L2) as fh:
        l2_header = fh.readline().rstrip("\n").split("\t")
        valid_ids = {l.split("\t", 1)[0] for l in fh}
    print(f"[0] spine: {len(valid_ids)} ORFs\n")

    df_best, df_sys = load_defensefinder(F_DF_HMMER, F_DF_GENES, valid_ids)
    thr = load_padloc_thresholds(F_PADLOC_HMM_META)
    pl_best, pl_sys = load_padloc(F_PADLOC_DOMTBL, thr, F_PADLOC_CSV, valid_ids)

    df_cols = ["df_system", "df_profile", "df_ievalue",
               "df_profile_cov", "df_seq_cov", "df_in_complete_system"]
    pl_cols = ["pl_system", "pl_profile", "pl_protein", "pl_evalue",
               "pl_hmm_cov", "pl_target_cov", "pl_in_complete_system"]

    print(f"\n[4] Writing {os.path.basename(F_GENE_TABLE_L3)}")
    n_df = n_pl = n_both = 0
    with open(F_GENE_TABLE_L2) as fin, open(F_GENE_TABLE_L3, "w") as out:
        next(fin)
        out.write("\t".join(l2_header + df_cols + pl_cols) + "\n")
        for line in fin:
            row = line.rstrip("\n").split("\t")
            pid = row[0]
            d, p = df_best.get(pid), pl_best.get(pid)
            if d:
                n_df += 1
            if p:
                n_pl += 1
            if d and p:
                n_both += 1
            row += [(d or {}).get(c, "") for c in df_cols[:-1]]
            row.append("1" if pid in df_sys else "0")
            row += [(p or {}).get(c, "") for c in pl_cols[:-1]]
            row.append("1" if pid in pl_sys else "0")
            out.write("\t".join(row) + "\n")

    print("\n" + "=" * 66)
    print(" VERIFICATION")
    print("=" * 66)
    print(f" ORFs in table                    : {len(valid_ids)}")
    print(f" with a DefenseFinder hit         : {n_df}")
    print(f" with a PADLOC hit (thresholded)  : {n_pl}")
    print(f" with BOTH                        : {n_both}")
    print(f" union (any specialist evidence)  : {n_df + n_pl - n_both}")
    print("-" * 66)
    print(f" in a complete DefenseFinder system: {len(df_sys)}")
    print(f" in a complete PADLOC system       : {len(pl_sys)}")
    orphan_base = n_df + n_pl - n_both
    if orphan_base:
        in_sys = len(df_sys | pl_sys)
        print(f" ORPHAN COMPONENTS                 : "
              f"{100.0 * (orphan_base - in_sys) / orphan_base:.1f}% "
              f"of defence-associated proteins are NOT part of a complete system")
    print("=" * 66)
    print(f" Output: {F_GENE_TABLE_L3}")


if __name__ == "__main__":
    main()
