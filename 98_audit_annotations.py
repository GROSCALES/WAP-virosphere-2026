#!/usr/bin/env python3
"""
98_audit_annotations.py — Comprehensive biological audit of gene classifications.

While structural QC ensures data integrity (row counts, non-duplication), this 
script validates the biological consistency of the pipeline: it asserts that 
every gene's final classification is supported by its underlying annotation.

Validation Checks:
  1. CORE-FUNCTION LEAKAGE (Strict): Asserts that no gene classified as Defence 
     carries an identity matching config.CORE_FUNCTION_REJECT. This enforces 
     the absolute integrity of the pre-classification veto mechanism. Exits 
     non-zero (fails the pipeline) if violated.
  2. SUSPECT VOCABULARY: Flags Defence genes whose annotations contain terms 
     associated with cellular replication, repair, translation, or structural 
     phage components for manual review.
  3. FULL DEFENCE INVENTORY: Compiles every distinct (category, KO name, Pfam 
     name, classification rule) tuple with its frequency.
  4. AVG INVENTORY: Compiles the equivalent exhaustive inventory for Auxiliary 
     Metabolic Genes (AVGs).
  5. REFERENCE COVERAGE: Asserts that every raw defence label maps successfully 
     to the explicitly declared semantic vocabulary.

Output: 98_audit_report.txt (and a non-zero exit code if Check 1 fails).
"""

import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

from config import OUT_DIR, CORE_FUNCTION_REJECT, DEFENCE_CATEGORIES, require

F_S3 = os.path.join(OUT_DIR, "Table_S3_AVG_Defence_genes.tsv")
F_MAP = os.path.join(OUT_DIR, "06_category_map.tsv")
F_OUT = os.path.join(OUT_DIR, "98_audit_report.txt")

VETO = re.compile(CORE_FUNCTION_REJECT, re.I)

# Words that should make you look twice at a "defence" gene. Deliberately
# broad: this is a prompt to check, not a filter.
SUSPECT = {
    "replication/repair": r"(DNA polymerase|primase|helicase|replicat|repair|"
                          r"recombinat|excision|topoisomeras|gyrase|ligase|"
                          r"single-strand|SOS response|RecA|LexA)",
    "ribosome/translation": r"(ribosom|tRNA synthetase|elongation factor|"
                            r"initiation factor)",
    "phage structural": r"(capsid|portal|terminase|tail|baseplate|virion|"
                        r"head protein|scaffold|major coat|minor coat|"
                        r"holin|endolysin|lysozyme|assembly)",
    "central metabolism": r"(dehydrogenase|synthetase|synthase|kinase|"
                          r"reductase|transferase family|ATP synthase)",
}
SUSPECT = {k: re.compile(v, re.I) for k, v in SUSPECT.items()}

_out = []
_fail = []


def say(s=""):
    print(s)
    _out.append(s)


def rule(title):
    say("\n" + "=" * 78)
    say(" " + title)
    say("=" * 78)


def ident(r):
    return f"{r.KO_name} {r.Pfam_name}".strip()


def path_of(r):
    """Which rule produced this call.

    Gene_class is tested FIRST. A protein that carries a DefenseFinder or
    PADLOC hit but was vetoed falls through to rule 3 and ends up an AVG
    while KEEPING its specialist columns, so testing those columns first
    would mislabel exactly the genes this audit cares most about.
    """
    if r.Gene_class == "AVG":
        return "rule 3 (VIBRANT AMG)"
    if r.DefenseFinder_raw_system or r.PADLOC_raw_system:
        return "rule 1 (specialist)"
    return "rule 2 (SEARCH_TERMS)"


def matches(series, rx):
    """Regex test that does not go through pandas' str.contains.

    CORE_FUNCTION_REJECT and the SUSPECT patterns are written with capture
    groups, and str.contains warns about that on every single call, which
    would bury the report in warnings. .search() is also plainer to read.
    """
    return series.map(lambda x: bool(rx.search(str(x))))


def main():
    require(F_S3, F_MAP)
    s3 = pd.read_csv(F_S3, sep="\t", low_memory=False, keep_default_na=False)
    cmap = pd.read_csv(F_MAP, sep="\t", keep_default_na=False)

    s3["identity"] = (s3.KO_name.astype(str) + " " +
                      s3.Pfam_name.astype(str)).str.strip()
    s3["path"] = s3.apply(path_of, axis=1)
    dfc = s3[s3.Gene_class == "Defence"]
    avg = s3[s3.Gene_class == "AVG"]

    say("=" * 78)
    say(" ANNOTATION AUDIT — every gene against its own annotation")
    say("=" * 78)
    say(f" Defence genes : {len(dfc)}")
    say(f" AVG genes     : {len(avg)}")
    say("\n by the rule that produced the call:")
    for k, v in s3.path.value_counts().items():
        say(f"   {v:>6}  {k}")

    # ---------------------------------------------------------------- 1
    rule("1. CORE-FUNCTION LEAK  (hard invariant: must be 0)")
    leak = dfc[matches(dfc.identity, VETO)]
    if len(leak):
        say(f"   [FAIL] {len(leak)} Defence genes carry a core-function identity")
        _fail.append("core-function leak")
        for (cat, ko, pf, p), n in Counter(
                zip(leak.Category, leak.KO_name, leak.Pfam_name,
                    leak.path)).most_common():
            say(f"     {n:>5}  [{p}] {cat}")
            say(f"            KO   : {ko[:90]}")
            say(f"            Pfam : {pf[:90]}")
    else:
        say("   [PASS] no Defence gene matches CORE_FUNCTION_REJECT.")
        say("          The veto is evaluated once per protein, before any")
        say("          rule, so this cannot be bypassed by any route.")

    # ---------------------------------------------------------------- 2
    rule("2. SUSPECT VOCABULARY  (look at every line; not automatically wrong)")
    any_flag = False
    for label, rx in SUSPECT.items():
        hit = dfc[matches(dfc.identity, rx)]
        if not len(hit):
            say(f"\n   {label}: none")
            continue
        any_flag = True
        say(f"\n   {label}: {len(hit)} genes")
        for (cat, ko, pf, p), n in Counter(
                zip(hit.Category, hit.KO_name, hit.Pfam_name,
                    hit.path)).most_common():
            say(f"     {n:>5}  [{p}] {cat}")
            say(f"            {(ko or pf)[:96]}")
    if not any_flag:
        say("\n   nothing flagged.")
    say("\n   NOTE: a genuine restriction endonuclease contains 'endonuclease';")
    say("   a genuine Gabija contains 'helicase'. This section is a reading")
    say("   list, not a verdict. What must NOT appear is a protein you cannot")
    say("   justify as part of the system it was assigned to.")

    # ---------------------------------------------------------------- 3
    rule("3. FULL DEFENCE INVENTORY  (complete; nothing truncated)")
    for cat in DEFENCE_CATEGORIES:
        sub = dfc[dfc.Category == cat]
        if not len(sub):
            continue
        say(f"\n --- {cat}  ({len(sub)} genes) " + "-" * max(0, 50 - len(cat)))
        for (ko, pf, p), n in Counter(
                zip(sub.KO_name, sub.Pfam_name, sub.path)).most_common():
            lbl = ko or pf or "(no annotation)"
            extra = f" | {pf[:44]}" if (ko and pf) else ""
            say(f"   {n:>5}  [{p:<20}] {lbl[:66]}{extra}")

    # ---------------------------------------------------------------- 4
    rule("4. AVG INVENTORY  (never audited before)")
    say(" AVGs come from VIBRANT's AMG flag (column 3 of its annotation file)")
    say(" plus a KEGG pathway lookup. The KEGG category is therefore only as")
    say(" good as VIBRANT's AMG call. Read this list the same way as above.")
    for cat, n_cat in avg.Category.value_counts().items():
        sub = avg[avg.Category == cat]
        say(f"\n --- {cat}  ({n_cat} genes) " + "-" * max(0, 46 - len(str(cat))))
        for (ko, pf), n in Counter(zip(sub.KO_name, sub.Pfam_name)).most_common():
            lbl = ko or pf or "(no annotation)"
            say(f"   {n:>5}  {lbl[:88]}")

    # ---------------------------------------------------------------- 5
    rule("5. REFERENCE COVERAGE")
    cats_seen = set(dfc.Category.unique())
    undeclared = cats_seen - set(DEFENCE_CATEGORIES)
    if undeclared:
        say(f"   [FAIL] categories not in DEFENCE_CATEGORIES: {sorted(undeclared)}")
        _fail.append("undeclared category")
    else:
        say(f"   [PASS] all {len(cats_seen)} categories are declared in config.")
    unmapped = cmap[cmap.category == "UNMAPPED"]
    if len(unmapped):
        say(f"   [FAIL] {len(unmapped)} raw labels UNMAPPED in the category map")
        _fail.append("unmapped label")
    else:
        say("   [PASS] every raw label resolves to a category.")

    say("\n   raw labels per category (from 06_category_map.tsv):")
    for cat, g in cmap.groupby("category"):
        say(f"     {int(g.n_genes.sum()):>6} genes / {len(g):>3} labels  {cat}")

    # ---------------------------------------------------------------- end
    with open(F_OUT, "w") as fh:
        fh.write("\n".join(_out) + "\n")
    say("\n" + "=" * 78)
    if _fail:
        say(f" {len(_fail)} CHECK(S) FAILED: {', '.join(_fail)}")
        say("=" * 78)
        sys.exit(1)
    say(" NO HARD FAILURES.")
    say(" Sections 2, 3 and 4 still need reading by a human — they are the")
    say(" part a script cannot decide.")
    say(f" Report written to {F_OUT}")
    say("=" * 78)


if __name__ == "__main__":
    main()
