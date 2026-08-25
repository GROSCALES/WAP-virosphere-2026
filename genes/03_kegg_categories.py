#!/usr/bin/env python3
"""
03_kegg_categories.py — Maps KOs to KEGG categories and measures ambiguity.

A KO can belong to several KEGG pathways. The primary category is defined 
as the one in which the KO occurs in the GREATEST NUMBER OF DISTINCT PATHWAYS 
(a property of the KO, not the pathway). Ties are broken by a declared list.
"""

import os
from collections import defaultdict

from config import F_VIB_AMG_PATH, OUT_DIR, require

F_OUT = os.path.join(OUT_DIR, "03_ko_kegg_categories.tsv")

# Tie-break only: applied when a KO occurs in the same number of pathways in
# two or more categories. Ordered from most to least specific, so a KO with
# a specific role is not swallowed by a broad category.
TIEBREAK = [
    "Energy metabolism",
    "Nucleotide metabolism",
    "Metabolism of cofactors and vitamins",
    "Amino acid metabolism",
    "Metabolism of other amino acids",
    "Lipid metabolism",
    "Glycan biosynthesis and metabolism",
    "Carbohydrate metabolism",
    "Metabolism of terpenoids and polyketides",
    "Biosynthesis of other secondary metabolites",
    "Xenobiotics biodegradation and metabolism",
    "Folding, sorting and degradation",
]
RANK = {c: i for i, c in enumerate(TIEBREAK)}


def main():
    require(F_VIB_AMG_PATH)
    print("=" * 70)
    print(" KO -> KEGG category, with ambiguity measured")
    print("=" * 70)

    ko_paths = defaultdict(list)          # ko -> [(metabolism, pathway)]
    n_rows = 0
    with open(F_VIB_AMG_PATH) as fh:
        header = fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            n_rows += 1
            metabolism = c[1].strip().strip('"')
            pathway = c[2].strip().strip('"')
            for ko in c[4].split(","):
                ko = ko.strip()
                if ko:
                    ko_paths[ko].append((metabolism, pathway))

    print(f"[1] pathway rows read : {n_rows}")
    print(f"    distinct KOs      : {len(ko_paths)}")

    # ---- assign ---------------------------------------------------------
    rows = []
    multi_cat = 0
    for ko, entries in sorted(ko_paths.items()):
        per_cat = defaultdict(set)
        for met, path in entries:
            per_cat[met].add(path)
        if len(per_cat) > 1:
            multi_cat += 1
        # most pathways wins; ties by declared order
        primary = sorted(
            per_cat.items(),
            key=lambda kv: (-len(kv[1]), RANK.get(kv[0], 99), kv[0])
        )[0][0]
        rows.append({
            "KO": ko,
            "kegg_metabolism_primary": primary,
            "n_categories": len(per_cat),
            "n_pathways": len({p for _, p in entries}),
            "kegg_metabolism_all": ";".join(sorted(per_cat)),
            "kegg_pathway_all": ";".join(sorted({p for _, p in entries})),
        })

    cols = ["KO", "kegg_metabolism_primary", "n_categories", "n_pathways",
            "kegg_metabolism_all", "kegg_pathway_all"]
    with open(F_OUT, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- report ---------------------------------------------------------
    pct = 100.0 * multi_cat / len(rows) if rows else 0
    print("\n" + "=" * 70)
    print(" AMBIGUITY")
    print("=" * 70)
    print(f" KOs in a single category   : {len(rows) - multi_cat}")
    print(f" KOs in MULTIPLE categories : {multi_cat}  ({pct:.1f}%)")
    print("   -> this percentage is what goes in Methods; it is the honest")
    print("      measure of how much the single-category figure simplifies.")

    print("\n Most ambiguous KOs (categories / pathways):")
    for r in sorted(rows, key=lambda r: (-r["n_categories"], -r["n_pathways"]))[:15]:
        print(f"   {r['KO']}  cats={r['n_categories']} paths={r['n_pathways']}"
              f"  -> {r['kegg_metabolism_primary']}")
        print(f"      all: {r['kegg_metabolism_all']}")

    per_primary = defaultdict(int)
    for r in rows:
        per_primary[r["kegg_metabolism_primary"]] += 1
    print("\n KOs assigned per primary category:")
    for cat, n in sorted(per_primary.items(), key=lambda kv: -kv[1]):
        print(f"   {n:5d}  {cat}")

    print("\n" + "=" * 70)
    print(f" Output: {F_OUT}")
    print(" Join this onto the gene table by KO in script 04.")


if __name__ == "__main__":
    main()
