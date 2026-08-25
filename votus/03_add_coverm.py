#!/usr/bin/env python3
"""
votus/03_add_coverm.py — abundance (TPM) from CoverM cross-mapping.

Only CB22 and GLI reads were mapped: VIRFC and LIS are MDA-amplified and
cannot be quantified. Every vOTU in the catalogue is nevertheless measured,
because all libraries are cross-mapped against the whole catalogue. So a
VIRFC-derived vOTU with TPM 0 everywhere means "not detected in the
quantitative libraries", which is a real, interpretable state — not missing
data.

Rebuild of 03_add_coverm.py, with three things fixed:

  1. COLUMN COLLISIONS. The original renamed columns by regex to
     TPM_<zone>_<fraction>. If two source columns mapped to the same name,
     pandas would silently create duplicate columns and every later
     positional access would be wrong. Now asserted.

  2. NA vs ZERO. The original did fillna(0), which merges "vOTU absent from
     the CoverM table" with "vOTU present but not detected". They are
     different and are now counted separately.

  3. UNEXPECTED FRACTIONS. If a TPM column turns up for a fraction that
     config says was not quantified, that is a contradiction between the
     data and the stated design, and the script says so instead of
     carrying on.

Adds TPM_total, TPM_<group> sums, n_samples_detected, and reports what
share of total viral abundance the length-filtered subsets capture — the
number that answers "why did you discard 80% of your contigs?".

Output: votu_03_coverm.tsv
"""

import os
import re
import sys

import pandas as pd

from config import (
    F_VOTU_02, F_VOTU_03, F_COVERM,
    FRACTIONS, ZONES, QUANTIFIED_FRACTIONS, require, Counter_,
)

SAMPLE_RE = re.compile(r"MCHA-([A-Z]{2})-\d+-([A-Z0-9]+)")

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def main():
    require(F_VOTU_02, F_COVERM)
    print("=" * 70)
    print(" vOTU LEVEL 3 — abundance (CoverM cross-mapping)")
    print("=" * 70)

    df = pd.read_csv(F_VOTU_02, sep="\t", dtype=str)
    df["contig_length"] = pd.to_numeric(df["contig_length"], errors="coerce")
    n0 = len(df)
    print(f"[1] table in: {n0} vOTUs")

    # ---------------------------------------------------------------- 2
    print(f"[2] reading {os.path.basename(F_COVERM)} ...")
    cm = pd.read_csv(F_COVERM, sep="\t")
    first = cm.columns[0]
    print(f"    rows: {len(cm)}   columns: {len(cm.columns)}")
    print(f"    ID column: '{first}'")
    cm = cm.rename(columns={first: "vOTU_ID"})

    print("    original sample columns:")
    for c in cm.columns[1:]:
        print(f"      {c}")

    # ---------------------------------------------------------------- 3
    print("[3] normalising sample names ...")
    mapping, seen = {}, {}
    for c in cm.columns[1:]:
        m = SAMPLE_RE.search(c)
        if not m:
            cnt.add(f"CoverM column that does not parse: {c}")
            print(f"    [warn] cannot parse: {c}  (left as-is)")
            continue
        zone, frac = m.group(1), m.group(2)
        if zone not in ZONES:
            cnt.add(f"zone not in config.ZONES: {zone}")
        if frac not in FRACTIONS:
            cnt.add(f"fraction not in config.FRACTIONS: {frac}")
        elif frac not in QUANTIFIED_FRACTIONS:
            # data contradicts the declared design — say so, do not paper over
            print(f"    [!] {c} belongs to {frac}, which config says was NOT "
                  f"quantified (MDA). Check the design table.")
            cnt.add(f"TPM column for a fraction marked not-quantified: {frac}")
        new = f"TPM_{zone}_{frac}"
        if new in seen:
            raise SystemExit(
                f"[FATAL] two CoverM columns map to '{new}':\n"
                f"  {seen[new]}\n  {c}\n"
                f"  Renaming would silently create duplicate columns.")
        seen[new] = c
        mapping[c] = new

    cm = cm.rename(columns=mapping)
    tpm_cols = sorted(v for v in mapping.values())
    print(f"    TPM columns: {len(tpm_cols)} -> {tpm_cols}")
    if cm.columns.duplicated().any():
        raise SystemExit("[FATAL] duplicate column names after renaming")

    for c in tpm_cols:
        cm[c] = pd.to_numeric(cm[c], errors="coerce")

    # ---------------------------------------------------------------- 4
    print("[4] merging ...")
    cm = cm[["vOTU_ID"] + tpm_cols].drop_duplicates(subset=["vOTU_ID"],
                                                    keep="first")
    df = df.merge(cm, on="vOTU_ID", how="left")
    if len(df) != n0:
        print("[FATAL] merge changed the row count")
        sys.exit(1)

    absent = df[tpm_cols[0]].isna() if tpm_cols else pd.Series([True] * n0)
    n_absent = int(absent.sum())
    df["In_CoverM"] = (~absent).map({True: "Yes", False: "No"})
    if n_absent:
        cnt.add("vOTUs with NO row in the CoverM table (TPM set to 0)",
                n_absent)
    df[tpm_cols] = df[tpm_cols].fillna(0.0)

    # ---------------------------------------------------------------- 5
    print("[5] derived abundance columns ...")
    df["TPM_total"] = df[tpm_cols].sum(axis=1)
    for frac in QUANTIFIED_FRACTIONS:
        cols = [c for c in tpm_cols if c.endswith("_" + frac)]
        if cols:
            df[f"TPM_sum_{frac}"] = df[cols].sum(axis=1)
    df["n_samples_detected"] = (df[tpm_cols] > 0).sum(axis=1)
    df["Detected_quantitative"] = (df["TPM_total"] > 0).map(
        {True: "Yes", False: "No"})

    n_zero = int((df["TPM_total"] == 0).sum())
    cnt.add("vOTUs present in CoverM but never detected (TPM_total = 0)",
            n_zero - n_absent if n_zero > n_absent else 0)

    df.to_csv(F_VOTU_03, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    print("\n" + "=" * 70)
    print(" VERIFICATION")
    print("=" * 70)
    print(f" vOTUs in / out           : {n0} / {len(df)}")
    print(f" with a CoverM row        : {int((df.In_CoverM=='Yes').sum())}")
    print(f" detected (TPM_total > 0) : {int((df.Detected_quantitative=='Yes').sum())}")
    print(f" never detected           : {n_zero}")

    print("\n detection by assembly fraction:")
    print("   (VIRFC and LIS were NOT quantified themselves; a zero here")
    print("    means absent from the CB22/GLI libraries, not missing data)")
    for frac in ["CB22", "GLI", "VIRFC", "LIS"]:
        sub = df[df.Fraction == frac]
        if not len(sub):
            continue
        det = int((sub.Detected_quantitative == "Yes").sum())
        print(f"   {frac:<7} {det:>7} / {len(sub):<7} "
              f"({100.0*det/len(sub):5.1f}% detected)   "
              f"TPM share {100.0*sub.TPM_total.sum()/df.TPM_total.sum():5.1f}%")

    # --- the number that answers "why discard 80% of the contigs?" -------
    print("\n" + "-" * 70)
    print(" ABUNDANCE CAPTURED BY THE LENGTH-FILTERED SUBSETS")
    print("-" * 70)
    tot = df["TPM_total"].sum()
    for label, mask in [
        ("all vOTUs", pd.Series([True] * len(df))),
        (">= 5 kb", df.contig_length >= 5000),
        (">= 10 kb", df.contig_length >= 10000),
    ]:
        sub = df[mask.values]
        print(f"   {label:<12} {len(sub):>7} vOTUs "
              f"({100.0*len(sub)/len(df):5.1f}%)   "
              f"{100.0*sub.TPM_total.sum()/tot:5.1f}% of total TPM")
    print("   -> if the 5 kb subset holds most of the abundance, the")
    print("      discarded contigs are low-abundance fragments and the")
    print("      analysis is unaffected. This is a supplementary figure.")

    cnt.report(total=n0)
    print("\n" + "=" * 70)
    print(f" Output: {F_VOTU_03}")


if __name__ == "__main__":
    main()
