#!/usr/bin/env python3
"""
genes/05_classify.py — the ONE place where a gene is called Defence, AVG or Other.

Everything before this script only gathered evidence into columns. Here a
single declared rule turns that evidence into a call, and records who said
what, so any row can be traced back.

PRECEDENCE
  1. specialist_defence : DefenseFinder or PADLOC assigned the protein to a
                          defence system, each using its OWN thresholds.
  2. vibrant_defence    : KO name or Pfam name matches a SEARCH_TERMS
                          category, with the e-value of the field that
                          actually matched <= EVALUE_MAX. Covers toxins,
                          virulence, phage communication, superinfection
                          exclusion and SOS modulation, which PADLOC and
                          DefenseFinder do not look for.
  3. vibrant_amg        : VIBRANT flagged the gene as an AMG (column 3).
  4. other

Rule 1 over rule 3 replaces the old hand-picked RECLASSIFY_DEFENSE list of
four KOs with a stated criterion.

Searches for rule 2 look ONLY at 'KO name' and 'Pfam name'. 'VOG name' is
the description of a reference protein from another organism and is never
used for identity — that was the source of every false-friend collision.

--------------------------------------------------------------------------
 THE CORE-FUNCTION VETO  (added after cross-checking the specialist output)
--------------------------------------------------------------------------
Cross-referencing every specialist call against the protein's OWN VIBRANT
annotation showed that ~25% of them were core cellular machinery captured
by HMM profiles of poorly characterised defence families:

    UG9        609 genes   527x Pfam "DNA polymerase family A", 523x KO polA
    UG6        166 genes   Toprim / DnaB-like helicase, KO dnaG, PEO1
    UG39        75 genes   VirE N-domain (primase)
    DEDDh_*    238 genes   75x KO dnaQ (Pol III proofreading epsilon subunit)
    cd09637    166 genes   98x KO K07465 RecB exonuclease; only 1 was cas4
    pfam13640  163 genes   110x KO K07336 / P4HA prolyl hydroxylase (metabolic)
    SspBCDE    168 genes   37x KO cysH (PAPS reductase) vs only 3x dndC
    Dodola     129 genes   AAA+ ATPases: clpA, clpC, cobS
    RM          96 genes   58x radD, 24x IRC3, 83x VOG uvsW

That last row is the SAME false-friend problem the original pipeline had:
radD, IRC3 and uvsW are repair helicases sharing PF04851 "Type III
restriction enzyme, res subunit". The old FALSOS_AMIGOS blacklist excluded
them correctly; PADLOC's generic "RM" profile pulls them back in.

So the veto below applies the same logic the original pipeline already
used, now to the specialists' output as well:

    a specialist defence call is REJECTED when the protein's own KEGG
    orthologue identifies it as a core replication, repair or metabolic
    enzyme (config.CORE_FUNCTION_REJECT).

Deliberately conservative, at a known cost. Cas4 is genuinely a RecB-family
exonuclease, so rejecting KO K07465 discards a handful of real Cas4 calls
along with ~98 generic exonucleases. With only 1 of 166 cd09637 hits
actually annotated cas4, erring towards rejection is the right trade: a
reviewer who opens the supplementary table and finds 609 "defence genes"
that are DNA polymerase I stops trusting the whole table.

Vetoed genes are NOT discarded. They fall through to rules 2 and 3, so a
cysH or cobS rejected here is picked up as an AVG, which is what it is.
Every veto is counted, and the reason is written into class_source.

Nothing is filtered. aux_score and passes_context ride along as columns so
the sensitivity analysis is a filter() in R.

Input : 03_gene_table_context.tsv , 03_ko_kegg_categories.tsv
Output: 04_gene_table_classified.tsv
"""

import os
import re
from collections import Counter

from config import (
    F_GENE_TABLE_CTX, F_GENE_TABLE_CLASS, OUT_DIR,
    SEARCH_TERMS, EVALUE_MAX, CORE_FUNCTION_REJECT, require,
)

F_KO_CATS = os.path.join(OUT_DIR, "03_ko_kegg_categories.tsv")
COMPILED = {cat: re.compile(pat, re.I) for cat, pat in SEARCH_TERMS.items()}
VETO = re.compile(CORE_FUNCTION_REJECT, re.I)


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_ko_categories(path):
    print(f"[1] KEGG categories from {os.path.basename(path)}")
    cats = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i = {c: j for j, c in enumerate(header)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(header):
                continue
            cats[c[i["KO"]]] = (c[i["kegg_metabolism_primary"]],
                                c[i["kegg_metabolism_all"]])
    print(f"    KOs mapped: {len(cats)}")
    return cats


def vibrant_defence_match(ko_name, ko_eval, pfam_name, pfam_eval):
    """Rule 2. Returns (category, matched_field, evalue) or None.

    The e-value checked belongs to the field that matched, not simply the
    first number in the row — which is what the old get_evalue() did,
    always returning the KO e-value even for a Pfam-driven call.
    """
    for cat, rx in COMPILED.items():
        if ko_name and rx.search(ko_name):
            ev = to_float(ko_eval)
            if ev is not None and ev <= EVALUE_MAX:
                return cat, "KO name", ev
        if pfam_name and rx.search(pfam_name):
            ev = to_float(pfam_eval)
            if ev is not None and ev <= EVALUE_MAX:
                return cat, "Pfam name", ev
    return None


def main():
    require(F_GENE_TABLE_CTX, F_KO_CATS)
    print("=" * 74)
    print(" CLASSIFICATION — one declared rule, full traceability")
    print("=" * 74)

    ko_cats = load_ko_categories(F_KO_CATS)

    with open(F_GENE_TABLE_CTX) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    i = {c: j for j, c in enumerate(header)}
    for c in ("ko", "is_amg", "ko_name", "ko_eval", "pfam_name", "pfam_eval",
              "df_system", "pl_system", "aux_score", "passes_context"):
        if c not in i:
            raise SystemExit(f"[SCHEMA ERROR] column '{c}' missing")

    new_cols = ["gene_class", "class_source", "category", "category_all",
                "evidence_sources", "match_field", "match_evalue",
                "specialist_vetoed"]

    n = 0
    by_class = Counter()
    by_source = Counter()
    def_cats = Counter()
    avg_cats = Counter()
    reclassified = Counter()
    ctx_ok = Counter()
    vetoed = Counter()           # KO gene symbol -> n vetoed
    vetoed_fam = Counter()       # specialist family -> n vetoed
    vetoed_to = Counter()        # where the vetoed gene ended up

    print("[2] classifying ...")
    with open(F_GENE_TABLE_CTX) as fin, open(F_GENE_TABLE_CLASS, "w") as out:
        next(fin)
        out.write("\t".join(header + new_cols) + "\n")
        for line in fin:
            c = line.rstrip("\n").split("\t")
            n += 1

            ko = c[i["ko"]].strip()
            ko_name = c[i["ko_name"]]
            is_amg = c[i["is_amg"]] == "1"
            df_sys = c[i["df_system"]].strip()
            pl_sys = c[i["pl_system"]].strip()

            gene_class, source = "Other", "none"
            category = category_all = ""
            evidence, match_field, match_eval = [], "", ""
            veto_flag = "No"

            # ================================================================
            # THE CORE-FUNCTION VETO IS A PROPERTY OF THE PROTEIN, NOT OF THE
            # RULE THAT HAPPENS TO FIRE.  <-- this is the fix.
            # ================================================================
            # It used to be evaluated only INSIDE the `if df_sys or pl_sys`
            # block, i.e. only for specialist calls. A protein with NO
            # DefenseFinder and NO PADLOC hit never had its identity checked
            # at all: it fell straight through to rule 2, where
            # SEARCH_TERMS["Restriction_Modification"] matched the pattern
            # "restriction enzyme" against Pfam PF04851, "Type III restriction
            # enzyme, res subunit".
            #
            # PF04851 is THE documented false friend of this pipeline: it is
            # shared between genuine Type III restriction res subunits and the
            # repair helicases radD, IRC3 and uvsW. So 123 DNA-repair proteins
            # were entering the largest defence category through the front
            # door while the veto guarded only the back one:
            #
            #     73  radD;    DNA repair protein RadD
            #     43  IRC3;    ATP-dependent helicase IRC3
            #      5  ERCC3, XPB; DNA excision repair protein ERCC-3
            #      2  SMARCAL1, HARP; SWI/SNF-related chromatin remodeller
            #
            # Found by auditing the genes that entered via rule 2 alone
            # (DefenseFinder and PADLOC columns both empty) against their own
            # annotations -- the same audit that had already been done for the
            # specialist path.
            #
            # Computing the veto ONCE, here, before any rule is dispatched,
            # makes it structurally impossible for a core-function protein to
            # be called Defence by any route. Genuine restriction enzymes are
            # unaffected: hsdR, hsdM, mcrA and "res; type III restriction
            # enzyme" carry no core-function orthologue, so the veto does not
            # fire for them.
            _id = f"{ko_name} {c[i['pfam_name']]}".strip()
            core_function = bool(_id and VETO.search(_id))
            if core_function:
                veto_flag = "Yes"
                fam = df_sys or pl_sys or "(no specialist call - rule 2 path)"
                vetoed_fam[fam] += 1
                sym = (ko_name.split(";")[0].strip()[:40] or
                       c[i["pfam_name"]].strip()[:40] or ko or "(no KO)")
                vetoed[sym] += 1
                # falls through to rule 3 only (see below)

            # ---- rule 1 : specialists, subject to the core-function veto
            if df_sys or pl_sys:
                if core_function:
                    pass          # vetoed above; do not make a specialist call
                else:
                    gene_class, source = "Defence", "specialist_defence"
                    category = df_sys or pl_sys
                    category_all = ";".join(
                        sorted({x for x in (df_sys, pl_sys) if x}))
                    if df_sys:
                        evidence.append("DefenseFinder")
                    if pl_sys:
                        evidence.append("PADLOC")
                    def_cats[category] += 1
                    if is_amg:
                        reclassified[ko or "(no KO)"] += 1

            # ---- rule 2 : VIBRANT KO / Pfam name -----------------------
            # FIXED — `and veto_flag == "No"` is the whole change.
            #
            # Without it a gene the veto had just rejected arrived here with
            # gene_class still == "Other", and rule 2 then searched THE SAME
            # ko_name / pfam_name that had triggered the veto one moment
            # earlier. SEARCH_TERMS["Restriction_Modification"] contains the
            # pattern "restriction enzyme", which matches the Pfam name
            # "Type III restriction enzyme, res subunit" -- PF04851, the very
            # domain that radD, IRC3 and uvsW share and that the docstring
            # above names as the original false friend.
            #
            # Effect in the last run: 96 genes (73x radD "DNA repair protein
            # RadD", 22x IRC3 helicase, 1x SMARCAL1 chromatin remodeller)
            # were vetoed by rule 1 and immediately re-admitted as Defence by
            # rule 2, every one of them landing in "Nucleic acid degradation".
            # All 96 carried evidence_sources == "VIBRANT", which is the
            # fingerprint of this path.
            #
            # The docstring's stated intent was always that vetoed genes fall
            # through to be picked up as METABOLIC genes ("a vetoed cysH or
            # cobS is picked up as an AVG, which is what it is") -- that is
            # rule 3. Rule 2 sat in between and silently re-made the exact
            # defence call the veto had refused. Vetoed genes now reach rule 3
            # only, as described. The 29 vetoed genes that legitimately become
            # AVGs (cysH, cobS, clpA...) are unaffected: rule 2 never matched
            # them in the first place.
            if gene_class == "Other" and veto_flag == "No":
                hit = vibrant_defence_match(ko_name, c[i["ko_eval"]],
                                            c[i["pfam_name"]],
                                            c[i["pfam_eval"]])
                if hit:
                    cat, field, ev = hit
                    gene_class, source = "Defence", "vibrant_defence"
                    category = category_all = cat
                    match_field, match_eval = field, f"{ev:.2e}"
                    evidence.append("VIBRANT")
                    def_cats[cat] += 1

            # ---- rule 3 : VIBRANT AMG flag -----------------------------
            if gene_class == "Other" and is_amg:
                gene_class, source = "AVG", "vibrant_amg"
                prim, allc = ko_cats.get(ko, ("Unclassified", "Unclassified"))
                category, category_all = prim, allc
                evidence.append("VIBRANT")
                avg_cats[prim] += 1

            if veto_flag == "Yes":
                source = f"{source}|specialist_vetoed_core_function"
                vetoed_to[gene_class] += 1

            by_class[gene_class] += 1
            by_source[source.split("|")[0]] += 1
            if c[i["passes_context"]] == "1":
                ctx_ok[gene_class] += 1

            c += [gene_class, source, category, category_all,
                  "+".join(evidence), match_field, match_eval, veto_flag]
            out.write("\t".join(c) + "\n")

    # ------------------------------------------------------------------
    print("\n" + "=" * 74)
    print(" RESULT")
    print("=" * 74)
    print(f" ORFs classified : {n}\n")
    for k in ("Defence", "AVG", "Other"):
        tot, ok = by_class[k], ctx_ok[k]
        pct = 100.0 * ok / tot if tot else 0
        print(f"   {k:<8} {tot:>8}     of which aux_score<=3: {ok} ({pct:.1f}%)")

    print("\n By rule:")
    for s in ("specialist_defence", "vibrant_defence", "vibrant_amg", "none"):
        print(f"   {s:<20} {by_source[s]:>8}")

    print("\n" + "-" * 74)
    print(" CORE-FUNCTION VETO — specialist calls overruled by the")
    print(" protein's own KEGG identity (see module docstring)")
    print("-" * 74)
    print(f"   total vetoed : {sum(vetoed.values())}")
    print("\n   by specialist family:")
    for fam, k in vetoed_fam.most_common(15):
        print(f"     {k:>6}  {fam}")
    print("\n   by the KEGG identity that overruled it:")
    for sym, k in vetoed.most_common(15):
        print(f"     {k:>6}  {sym}")
    print("\n   where the vetoed genes ended up:")
    for cls, k in vetoed_to.most_common():
        print(f"     {k:>6}  {cls}")
    print("   (nothing is lost: a vetoed cysH or cobS is picked up as an AVG,")
    print("    which is what it is. NOTE: since the rule-2 fix, no vetoed gene")
    print("    can come back as Defence -- if any 'Defence' appears in the")
    print("    list above, something is wrong and it must be investigated.)")
    if vetoed_to.get("Defence"):
        print(f"\n   [!] {vetoed_to['Defence']} vetoed genes are still Defence "
              f"-- this should be 0 after the rule-2 fix.")

    print("\n" + "-" * 74)
    print(" Top DEFENCE categories (raw labels; harmonised in script 06):")
    for cat, k in def_cats.most_common(20):
        print(f"   {k:>7}  {cat}")

    print("\n AVG categories (primary):")
    for cat, k in avg_cats.most_common():
        print(f"   {k:>7}  {cat}")

    print("\n" + "-" * 74)
    print(" AMG -> DEFENCE reclassified automatically by rule 1")
    print(" (replaces the hand-picked RECLASSIFY_DEFENSE list)")
    print("-" * 74)
    for ko, k in reclassified.most_common(25):
        print(f"   {k:>7}  {ko}")
    print(f"   TOTAL: {sum(reclassified.values())} genes, "
          f"{len(reclassified)} distinct KOs")

    print("\n" + "=" * 74)
    print(f" Output: {F_GENE_TABLE_CLASS}")


if __name__ == "__main__":
    main()
