#!/usr/bin/env python3
"""
votus/04_add_phabox_gtdb.py — PhaBOX taxonomy/host, and GTDB-Tk for MAGs.

Integrates PhaBOX predictions (taxonomy, lifestyle, host) and GTDB-Tk 
classifications for MAGs. Both bacterial and archaeal GTDB summaries are 
read and concatenated to ensure archaeal MAGs retain their taxonomy.

MAG conflict resolution: All MAGs linked to a vOTU are resolved. If a vOTU 
maps to several MAGs with different taxonomies, the conflict is explicitly 
flagged in MAG_tax_conflict. MAG-derived host predictions are prioritized 
over PhaBOX database inferences, as they represent physically assembled 
genomes from the same environment.

Output: votu_04_phabox_gtdb.tsv
"""

import os
import sys

import pandas as pd

from config import (
    F_VOTU_03, F_VOTU_04, F_PHABOX, F_GTDB_BAC, F_GTDB_AR,
    require, Counter_,
)

NULLS = {"-", "", "nan", "none", "unknown", "filtered", "na", "n/a"}

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def norm(x):
    return "Unknown" if str(x).strip().lower() in NULLS else str(x).strip()


def load_gtdb():
    """Read bacterial and archaeal summaries and concatenate them."""
    print("[3] GTDB-Tk ...")
    frames = []
    for label, path in (("bacteria", F_GTDB_BAC), ("archaea", F_GTDB_AR)):
        if not os.path.exists(path):
            print(f"    [warn] missing {label} file: {path}")
            cnt.add(f"GTDB {label} summary missing")
            continue
        t = pd.read_csv(path, sep="\t", dtype=str)
        if "user_genome" not in t.columns or "classification" not in t.columns:
            raise SystemExit(
                f"[SCHEMA ERROR] {os.path.basename(path)} lacks "
                f"user_genome/classification; found {list(t.columns)[:8]}")
        t = t[["user_genome", "classification"]].copy()
        t["domain_file"] = label
        print(f"    {label:<9} {len(t):>5} MAGs")
        frames.append(t)
    if not frames:
        raise SystemExit("[FATAL] no GTDB-Tk summary could be read")
    g = pd.concat(frames, ignore_index=True)
    dup = int(g["user_genome"].duplicated().sum())
    if dup:
        cnt.add("MAGs appearing in BOTH GTDB files (kept first)", dup)
        g = g.drop_duplicates(subset=["user_genome"], keep="first")
    print(f"    total unique MAGs: {len(g)}")
    return dict(zip(g["user_genome"], g["classification"]))


def main():
    require(F_VOTU_03, F_PHABOX)
    print("=" * 70)
    print(" vOTU LEVEL 4 — PhaBOX taxonomy/host + GTDB-Tk for MAG-derived")
    print("=" * 70)

    df = pd.read_csv(F_VOTU_03, sep="\t", dtype=str)
    n0 = len(df)
    print(f"[1] table in: {n0} vOTUs")

    # ---------------------------------------------------------------- 2
    print(f"[2] reading {os.path.basename(F_PHABOX)} ...")
    pb = pd.read_csv(F_PHABOX, sep="\t", dtype=str)
    print(f"    rows: {len(pb)}   columns: {list(pb.columns)}")

    want = {
        "Accession": "vOTU_ID",
        "Pred": "PhaBOX_Pred",
        "PhaMerConfidence": "PhaBOX_Confidence",
        "Lineage": "PhaBOX_Taxonomy",
        "TYPE": "PhaBOX_Lifestyle",
        "Host": "PhaBOX_Host_Species",
        "Host_NCBI_lineage": "PhaBOX_Host_NCBI",
        "Host_GTDB_lineage": "PhaBOX_Host_GTDB",
    }
    missing = [c for c in want if c not in pb.columns]
    if missing:
        raise SystemExit(f"[SCHEMA ERROR] PhaBOX lacks columns: {missing}")

    pb = pb[list(want)].rename(columns=want)
    dup = int(pb["vOTU_ID"].duplicated().sum())
    if dup:
        cnt.add("duplicate vOTU_IDs in the PhaBOX output (kept first)", dup)
        pb = pb.drop_duplicates(subset=["vOTU_ID"], keep="first")

    for c in pb.columns:
        if c != "vOTU_ID":
            n_null = int(pb[c].apply(lambda x: str(x).strip().lower()
                                     in NULLS).sum())
            pb[c] = pb[c].apply(norm)
            if n_null:
                cnt.add(f"PhaBOX '{c}' empty/filtered", n_null)

    df = df.merge(pb, on="vOTU_ID", how="left")
    if len(df) != n0:
        print("[FATAL] PhaBOX merge changed the row count")
        sys.exit(1)
    for c in want.values():
        if c != "vOTU_ID":
            df[c] = df[c].fillna("Unknown")

    n_pb = int((df["PhaBOX_Pred"] != "Unknown").sum())
    print(f"    vOTUs with a PhaBOX prediction: {n_pb} "
          f"({100.0*n_pb/n0:.1f}%)")

    # ---------------------------------------------------------------- 3
    gtdb = load_gtdb()

    print("[4] resolving MAG taxonomy (all MAGs, not just the first) ...")

    def resolve(mag_ids):
        s = str(mag_ids).strip()
        if s in ("None", "nan", ""):
            return pd.Series(["Unknown", "No", 0])
        mags = [m for m in s.split(";") if m]
        taxa = []
        for m in mags:
            t = gtdb.get(m)
            if t is None:
                cnt.add("MAG absent from the GTDB-Tk summaries")
            else:
                taxa.append(t)
        if not taxa:
            return pd.Series(["Unknown_MAG_Tax", "No", len(mags)])
        uniq = sorted(set(taxa))
        if len(uniq) == 1:
            return pd.Series([uniq[0], "No", len(mags)])
        cnt.add("vOTUs whose MAGs disagree on taxonomy")
        return pd.Series([";".join(uniq), "Yes", len(mags)])

    df[["MAG_taxonomy", "MAG_tax_conflict", "MAG_resolved_n"]] = \
        df["MAG_ids"].apply(resolve)

    n_mag = int((df["MAG_taxonomy"] != "Unknown").sum())
    n_ok = int((~df["MAG_taxonomy"].isin(["Unknown", "Unknown_MAG_Tax"])).sum())
    print(f"    vOTUs linked to a MAG        : {n_mag}")
    print(f"    of those, taxonomy resolved  : {n_ok}")

    # ---------------------------------------------------------------- 5
    print("[5] consolidating host prediction ...")
    # MAG-derived host beats PhaBOX: it is a physically assembled genome
    # from the same samples, not a database inference.
    def host(r):
        if r["MAG_taxonomy"] not in ("Unknown", "Unknown_MAG_Tax"):
            return "MAG: " + r["MAG_taxonomy"], "MAG"
        if r["PhaBOX_Host_GTDB"] != "Unknown":
            return "PhaBOX: " + r["PhaBOX_Host_GTDB"], "PhaBOX_GTDB"
        if r["PhaBOX_Host_NCBI"] != "Unknown":
            return "PhaBOX: " + r["PhaBOX_Host_NCBI"], "PhaBOX_NCBI"
        return "Unknown", "None"

    res = df.apply(host, axis=1, result_type="expand")
    df["Host_Prediction"] = res[0]
    df["Host_Source"] = res[1]

    df.to_csv(F_VOTU_04, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    print("\n" + "=" * 70)
    print(" VERIFICATION")
    print("=" * 70)
    print(f" vOTUs in / out : {n0} / {len(df)}")
    print("\n host prediction source:")
    for k, v in df["Host_Source"].value_counts().items():
        print(f"   {v:>8}  {k}   ({100.0*v/n0:.1f}%)")
    print("\n PhaBOX lifestyle:")
    for k, v in df["PhaBOX_Lifestyle"].value_counts().head(8).items():
        print(f"   {v:>8}  {k}")
    print("\n top MAG-derived hosts:")
    sub = df[df["Host_Source"] == "MAG"]
    for k, v in sub["MAG_taxonomy"].value_counts().head(10).items():
        print(f"   {v:>5}  {k[:90]}")

    cnt.report(total=n0)
    print("\n" + "=" * 70)
    print(f" Output: {F_VOTU_04}")


if __name__ == "__main__":
    main()
