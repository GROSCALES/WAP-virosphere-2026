#!/usr/bin/env python3
"""
votus/06_filter_catalog.py — The final catalogue filter and cascade audit.

Applies the structural and evidence criteria to define the final viral catalogue:
  1. Viral evidence: Must possess >=1 independent line of evidence (CheckV provirus, 
     CT3 hallmark, VIBRANT prediction, PhaBOX prediction, or high CheckV quality).
  2. Length/Topology: >= 10 kb (conservative) OR >= 5 kb (main) OR circular OR 
     short but rescued (possesses both a hallmark and provirus linkage).
  3. Domain: Exclusion of eukaryotic viruses and virophages.

Lifestyle assignment follows a strict three-tier evidence hierarchy:
  - Temperate (integrated): Physical evidence (CheckV/cluster provirus).
  - Temperate (predicted): Algorithmic prediction (VIBRANT or PhaBOX).
  - Lytic: No temperate signal detected (cannot strictly exclude pseudolysogeny).

Nothing is dropped from the table. The full table is kept with 'Passes_catalog' 
and 'Conservative_10kb' boolean flags, so any filter can be revisited in R 
without re-running the pipeline.

Output: votu_06_catalog.tsv 
"""

import os
import sys

import pandas as pd

from config import (
    F_VOTU_05, F_VOTU_06, LEN_MAIN_MIN, LEN_CONSERVATIVE_MIN,
    require, Counter_,
)

# Eukaryotic viruses and virophages — not the target of this study.
EUK_KEYWORDS = [
    'nucleocytoviricota', 'megaviricetes', 'mimiviridae', 'pithovirus',
    'pithoviruses', 'medusavirus', 'kaumoebavirus', 'tetraselmis',
    'theiavirus', 'imitervirales', 'cressdnaviricota', 'shotokuvirae',
    'arfiviricetes', 'circoviridae', 'parvoviridae', 'cossaviricota',
    'crab', 'ctenophore', 'pond-associated circular dna', 'maveriviricetes',
    'lavidaviridae', 'lavidavirales', 'virophage', 'virophaviricetes',
    'sputniviroviridae', 'burtonviroviridae', 'mividavirales',
    'priklausovirales', 'sputnikvirus', 'polintoviricetes', 'adintoviridae',
    'aquintoviricetes', 'archintovirales', 'amphintovirales',
    'naldaviricetes', 'lefavirales', 'orthopolintovirales', 'phypoliviridae',
    'pharingeaviricetes', 'rowavirales', 'polisuviricotina',
    'alphadintovirus', 'betadintovirus',
]

cnt = Counter_("COUNTERS — nothing is dropped in silence")


def main():
    require(F_VOTU_05)
    print("=" * 70)
    print(" vOTU LEVEL 6 — catalogue filter (flags, full table preserved)")
    print("=" * 70)

    df = pd.read_csv(F_VOTU_05, sep="\t", dtype=str)
    df["contig_length"] = pd.to_numeric(df["contig_length"], errors="coerce").fillna(0)
    for c in ("CT3_Virion_Hallmarks", "CT3_Rep_Hallmarks", "CT3_RDRP_Hallmarks"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    n0 = len(df)
    print(f"[1] table in: {n0} vOTUs")

    # ---------------------------------------------------------------- 2
    print("[2] independent lines of evidence ...")
    ev = pd.DataFrame({
        "VIBRANT": df["Vibrant_predicted"] == "Yes",
        "PhaBOX": df["PhaBOX_Pred"] != "Unknown",
        "CT3_hallmark": df["CT3_Any_Hallmark"] == "Yes",
        "CheckV_quality": df["checkv_quality"].str.lower().isin(
            ["medium-quality", "high-quality", "complete"]),
        "Provirus": (df["CheckV_Provirus"].str.strip() == "Yes") |
                    (df["Cluster_has_Provirus"] == "Yes"),
    })
    df["N_evidence"] = ev.sum(axis=1)
    df["Evidence_sources"] = ev.apply(
        lambda r: "+".join([c for c in ev.columns if r[c]]) or "None", axis=1)
    df["Has_evidence"] = (df["N_evidence"] > 0).map({True: "Yes", False: "No"})

    # ---------------------------------------------------------------- 3
    print("[3] length / rescue ...")
    long_enough = df["contig_length"] >= LEN_MAIN_MIN
    circular = (df["is_circular_checkv"] == "Yes") | (df["CT3_Circular"] == "Yes")
    hallmark = (df["CT3_Virion_Hallmarks"] > 0) | (df["CT3_Rep_Hallmarks"] > 0)
    # "golden fragment": short, but carries a hallmark AND is provirus-linked.
    # NOTE the change of criterion. The original rescued short contigs that
    # carried an AMG, which is circular reasoning in a paper whose result is
    # about AMGs: you keep the contig because it has one, then report it.
    # Provirus evidence is independent of gene content.
    provirus = ev["Provirus"]
    rescue = (~long_enough) & hallmark & provirus
    df["Passes_length"] = (long_enough | circular | rescue).map(
        {True: "Yes", False: "No"})
    df["Rescued_short"] = rescue.map({True: "Yes", False: "No"})
    print(f"    >= {LEN_MAIN_MIN} bp        : {int(long_enough.sum())}")
    print(f"    circular (rescue)   : {int((circular & ~long_enough).sum())}")
    print(f"    short + hallmark + provirus (rescue): {int(rescue.sum())}")
    cnt.add("short contigs rescued (hallmark + provirus)", int(rescue.sum()))

    # ---------------------------------------------------------------- 4
    print("[4] eukaryotic viruses and virophages ...")
    tax_low = df["Taxonomy"].fillna("").str.lower()
    euk = pd.Series(False, index=df.index)
    for kw in EUK_KEYWORDS:
        hit = tax_low.str.contains(kw, regex=False)
        if hit.any():
            cnt.add(f"excluded as eukaryotic/virophage: {kw}", int(hit.sum()))
        euk |= hit
    df["Is_prokaryotic"] = (~euk).map({True: "Yes", False: "No"})
    print(f"    flagged eukaryotic/virophage: {int(euk.sum())}")

    # ---------------------------------------------------------------- 5
    print("[5] lifestyle ...")

    def lifestyle(r):
        if str(r["CheckV_Provirus"]).strip() == "Yes":
            return pd.Series(["Temperate (integrated)", "CheckV provirus"])
        if r["Cluster_has_Provirus"] == "Yes":
            return pd.Series(["Temperate (integrated)", "cluster provirus"])
        srcs = []
        if str(r["Vibrant_Lifestyle"]).lower() == "lysogenic":
            srcs.append("VIBRANT")
        if str(r["PhaBOX_Lifestyle"]).lower() == "temperate":
            srcs.append("PhaBOX")
        if srcs:
            return pd.Series(["Temperate (predicted)", "+".join(srcs)])
        return pd.Series(["Lytic", "no temperate signal"])

    df[["Lifestyle", "Lifestyle_Source"]] = df.apply(lifestyle, axis=1)

    # ---------------------------------------------------------------- 6
    df["Passes_catalog"] = (
        (df["Has_evidence"] == "Yes") &
        (df["Passes_length"] == "Yes") &
        (df["Is_prokaryotic"] == "Yes")).map({True: "Yes", False: "No"})
    df["Conservative_10kb"] = (
        (df["Passes_catalog"] == "Yes") &
        ((df["contig_length"] >= LEN_CONSERVATIVE_MIN) |
         (df["is_circular_checkv"] == "Yes") |
         (df["CT3_Circular"] == "Yes"))).map({True: "Yes", False: "No"})

    df.to_csv(F_VOTU_06, sep="\t", index=False)

    # ---------------------------------------------------------------- report
    m_ev = df["Has_evidence"] == "Yes"
    m_len = df["Passes_length"] == "Yes"
    m_pro = df["Is_prokaryotic"] == "Yes"
    n1 = int(m_ev.sum())
    n2 = int((m_ev & m_len).sum())
    n3 = int((m_ev & m_len & m_pro).sum())

    print("\n" + "=" * 70)
    print(" FILTERING CASCADE")
    print("=" * 70)
    print(f"   [INPUT]  vOTUs                       {n0:>8}")
    print(f"   [-] no viral evidence               -{n0-n1:>8}")
    print(f"   [-] too short and not rescued       -{n1-n2:>8}")
    print(f"   [-] eukaryotic / virophage          -{n2-n3:>8}")
    print("   " + "-" * 44)
    print(f"   [OUTPUT] FINAL CATALOGUE             {n3:>8}")
    print(f"            conservative (>=10 kb)      "
          f"{int((df.Conservative_10kb=='Yes').sum()):>8}")

    fin = df[df["Passes_catalog"] == "Yes"]

    print("\n" + "-" * 70)
    print(" EVIDENCE WITHIN THE FINAL CATALOGUE")
    print("-" * 70)
    for k in range(1, 6):
        n = int((fin["N_evidence"] == k).sum())
        if n:
            print(f"   {n:>7}  {k} independent source(s)   "
                  f"({100.0*n/len(fin):.1f}%)")
    n_multi = int((fin["N_evidence"] >= 2).sum())
    print(f"\n   >=2 sources: {n_multi} ({100.0*n_multi/len(fin):.1f}%)  "
          f"<- the number to put in Methods")

    print("\n most frequent evidence combinations:")
    for k, v in fin["Evidence_sources"].value_counts().head(8).items():
        print(f"   {v:>7}  {k}")

    print("\n" + "-" * 70)
    print(" LIFESTYLE (final catalogue)")
    print("-" * 70)
    for k, v in fin["Lifestyle"].value_counts().items():
        print(f"   {v:>7}  {k}")
    print("\n by culture group:")
    tab = pd.crosstab(fin["Fraction_group"], fin["Lifestyle"])
    print(tab.to_string())
    print("\n by fraction (the mitomycin control is GLI vs CB22):")
    tab2 = pd.crosstab(fin["Fraction"], fin["Lifestyle"], normalize="index") * 100
    print(tab2.round(1).to_string())

    print("\n" + "-" * 70)
    print(" FINAL CATALOGUE COMPOSITION")
    print("-" * 70)
    for col in ("Fraction", "Fraction_group", "Zone"):
        print(f"\n {col}:")
        for k, v in fin[col].value_counts().items():
            print(f"   {v:>7}  {k}")

    cnt.report(total=n0)
    print("\n" + "=" * 70)
    print(f" Output: {F_VOTU_06}")
    print(" The FULL table is written, with Passes_catalog / Conservative_10kb")
    print(" as flags. Filtering happens in R at reporting time, so changing")
    print(" the criteria never means re-running the pipeline.")


if __name__ == "__main__":
    main()
