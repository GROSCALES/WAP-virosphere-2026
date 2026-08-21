#!/usr/bin/env python3
"""
genes/06_harmonise_categories.py — one vocabulary for the defence categories.

Resolution chain, in order, for every raw label found in the gene table:

  1. DEFENCE_EXCLUDE          -> dropped from defence, reason recorded
  2. OWN_CATEGORIES           -> our SEARCH_TERMS labels; we know what they are
  3. PADLOC sys/*.yaml        -> profile family -> component protein -> system
                                 (case-insensitive: comparing literally failed
                                 on cd09637 -> "Cas4" vs the yaml's "cas4")
  4. sys_meta.txt name match  -> DefenseFinder names vs PADLOC yaml names
  5. SYSTEM_TO_CATEGORY, matched by LONGEST PREFIX, case-insensitive, over
     three candidates IN THIS ORDER:
        a) the resolved system name
        b) the PROTEIN names behind the profile, from hmm_meta
        c) the raw label itself

Order (b) before (c) matters. Some PADLOC profiles are named after a raw
database accession — cd09637, pfam13640, COG1468 — which carries no meaning,
and matching those on their prefix sent 329 CRISPR-associated proteins
(Cas4, 2OG) into "Uncharacterised". Their protein.name does carry meaning,
so it is consulted first. Accession prefixes remain in the map only as a
last resort.

Order (c) is what catches DefenseFinder's own system names (Gabija, Lamassu,
CBASS, Pycsar, UG9, DS-4...), which never appear in PADLOC's files at all.

The reasoning behind the category vocabulary lives in config.py section 7b:
it is a scientific decision, so it belongs with the data, in one place.

Output: 06_category_map.tsv   (category_raw -> system -> category)
"""

import os
import re
import glob
from collections import defaultdict, Counter

from config import (
    F_PADLOC_SYS_DIR, F_PADLOC_HMM_META, F_GENE_TABLE_CLASS, OUT_DIR,
    DEFENCE_EXCLUDE, SYSTEM_TO_CATEGORY, LABEL_OVERRIDE, DEFENCE_CATEGORIES,
    require,
)

F_SYS_META = os.path.join(os.path.dirname(F_PADLOC_SYS_DIR), "sys_meta.txt")
F_OUT = os.path.join(OUT_DIR, "06_category_map.tsv")

TRAILING_NUM = re.compile(r"_\d+$")

OWN_CATEGORIES = {
    "Restriction_Modification": "Restriction-Modification",
    "DNA_Methylation": "DNA methylation / anti-restriction",
    "Toxin_Antitoxin_and_Abi": "Toxin-Antitoxin / Abi",
    "Toxins_and_Virulence_Factors": "Toxins and virulence factors",
    "Phage_Communication": "Phage communication",
    "Superinfection_Exclusion": "Superinfection exclusion",
    "Phage_Nucleus_Assembly": "Phage nucleus",
    "AMR_and_Metal_Resistance": "AMR / metal resistance",
    "SOS_response_modulation": "SOS response modulation",
}

# longest key first: "DEDDh_I_II" must win over "DE", "RM Type" over "RM"
_KEYS = sorted(SYSTEM_TO_CATEGORY, key=len, reverse=True)
_OVR = sorted(LABEL_OVERRIDE, key=len, reverse=True)


def override_category(raw):
    """LABEL_OVERRIDE, matched by longest prefix on the RAW label only.

    Consulted BEFORE the normal resolution chain, because the labels it
    covers are precisely the ones the chain gets wrong: cd09637 and COG1468
    resolve through hmm_meta's protein.name to "Cas4", and TnsC through
    PADLOC's yaml to "CRISPR-Cas Type I-F3", so neither would ever reach the
    raw-label step where an ordinary SYSTEM_TO_CATEGORY entry could catch
    them. See the reasoning in config.py section 7b-ter.
    """
    low = str(raw).strip().lower()
    for k in _OVR:
        if low.startswith(k.lower()):
            return LABEL_OVERRIDE[k], k, raw
    return None, None, None


# --------------------------------------------------------------------------
# FIXED — choosing between several PADLOC systems for the same component.
#
# A PADLOC component often belongs to more than one system, so comp2sys
# returns a SET. The previous code resolved that with sorted(hits)[0], i.e.
# whichever system came first ALPHABETICALLY, with no rule stated anywhere.
#
# That silently contradicted the policy declared in 07_join_and_report.py:
#
#     "CATEGORIES ARE SYSTEM-LEVEL, NOT GENE-LEVEL. [...] an RM system
#      comprises a nucleolytic R subunit and a methyltransferase M subunit,
#      and both are filed under the system's primary mechanism."
#
# because PADLOC's catch-all bucket is literally named "DNA-modification
# systems (other)", and D sorts before R. So:
#
#     RM_Type_II  -> "RM Type II"                       -> Nucleic acid degradation
#     MTase_II    -> "DNA-modification systems (other)" -> Genome modification
#     REase_II    -> "DNA-modification systems (other)" -> Genome modification
#
# The generic system label went to degradation while its own named subunits
# went to modification -- and REase_II, the NUCLEOLYTIC subunit, ended up
# under "modification", which is wrong under either possible policy.
#
# The rule now: a specifically named system always beats a generic
# "(other)" catch-all bucket. Alphabetical order is kept only as the
# tie-break between two equally specific systems, so the function stays
# deterministic. This makes the code do what the docstring already says.
#
# Effect: MTase_II, REase_II and mREase_IV resolve to "RM Type II" /
# "RM Type IV" -> Nucleic acid degradation. oMTase (an ORPHAN
# methyltransferase, genuinely not part of an RM system) has only the
# generic hit and correctly stays in Genome modification. It also fixes
# two inconsistencies as a side effect: PbeA -> PbeABCD and mads -> MADS
# now agree with their explicit config entries and with the way
# DefenseFinder's own "MADS" label was already being categorised.
# --------------------------------------------------------------------------
def _system_rank(system):
    """Sort key: named systems first, generic '(other)' buckets last."""
    return (1 if "(other)" in str(system).lower() else 0, str(system))


def to_category(*candidates):
    """Longest-prefix, case-insensitive match against SYSTEM_TO_CATEGORY.

    Returns (category, matched_key, matched_candidate) or (None, None, None).
    """
    for cand in candidates:
        if not cand:
            continue
        low = str(cand).strip().lower()
        for k in _KEYS:
            if low.startswith(k.lower()):
                return SYSTEM_TO_CATEGORY[k], k, cand
    return None, None, None


def parse_yaml_lists(path):
    """Minimal parser for PADLOC's flat system definitions: they only ever
    contain 'key: value' and 'key:' followed by '- item', so no YAML
    dependency is needed."""
    out, key = defaultdict(list), None
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip() or s.lstrip().startswith("#"):
                continue
            if s.lstrip().startswith("- "):
                if key:
                    v = s.lstrip()[2:].strip().strip('"').strip("'")
                    if v and v != "NA":
                        out[key].append(v)
            elif ":" in s:
                k, _, v = s.partition(":")
                key = k.strip()
                v = v.strip().strip('"').strip("'")
                if v and v != "NA":
                    out[key].append(v)
                    key = None
    return out


def load_sys_meta(path):
    """yaml.name -> readable system name."""
    m = {}
    if not os.path.exists(path):
        print(f"    [warn] no sys_meta.txt at {path}")
        return m
    with open(path) as fh:
        header = [c.strip() for c in fh.readline().rstrip("\n").split("\t")]
        idx = {c: i for i, c in enumerate(header)}
        if "system" not in idx or "yaml.name" not in idx:
            print(f"    [warn] sys_meta.txt lacks system/yaml.name: {header}")
            return m
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > max(idx["system"], idx["yaml.name"]):
                m[c[idx["yaml.name"]].strip()] = c[idx["system"]].strip()
    return m


def main():
    require(F_GENE_TABLE_CLASS)
    print("=" * 78)
    print(" CATEGORY HARMONISATION")
    print("=" * 78)

    sys_meta = load_sys_meta(F_SYS_META)
    print(f"[1] sys_meta systems: {len(sys_meta)}")

    # --- components -> system, indexed lowercase --------------------------
    comp2sys = defaultdict(set)
    yamls = sorted(glob.glob(os.path.join(F_PADLOC_SYS_DIR, "*.yaml")))
    for y in yamls:
        name = os.path.basename(y)[:-5]
        system = sys_meta.get(name, name)
        d = parse_yaml_lists(y)
        for k in ("core_genes", "secondary_genes"):
            for gene in d.get(k, []):
                comp2sys[gene.strip().lower()].add(system)
    print(f"[2] yaml files: {len(yamls)}   components: {len(comp2sys)}")

    # --- profile family -> protein names ----------------------------------
    fam2prot = defaultdict(set)
    if os.path.exists(F_PADLOC_HMM_META):
        with open(F_PADLOC_HMM_META) as fh:
            header = [c.strip() for c in fh.readline().rstrip("\n").split("\t")]
            idx = {c: i for i, c in enumerate(header)}
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) <= max(idx["hmm.name"], idx["protein.name"]):
                    continue
                fam = TRAILING_NUM.sub("", c[idx["hmm.name"]])
                for p in str(c[idx["protein.name"]]).split("|"):
                    p = p.strip()
                    if p and p != "NA":
                        fam2prot[fam].add(p)
    print(f"[3] profile families: {len(fam2prot)}")

    # --- observed labels --------------------------------------------------
    with open(F_GENE_TABLE_CLASS) as fh:
        gh = fh.readline().rstrip("\n").split("\t")
    gi = {c: i for i, c in enumerate(gh)}
    seen = Counter()
    with open(F_GENE_TABLE_CLASS) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[gi["gene_class"]] == "Defence":
                seen[c[gi["category"]]] += 1
    print(f"[4] distinct labels: {len(seen)}   genes: {sum(seen.values())}")

    # --- resolve ----------------------------------------------------------
    rows, unmapped = [], []
    for raw, n in seen.most_common():
        if raw in DEFENCE_EXCLUDE:
            rows.append((raw, n, "-", "EXCLUDED (not defence)",
                         "DEFENCE_EXCLUDE"))
            continue

        system, how = None, ""
        if raw in OWN_CATEGORIES:
            system, how = OWN_CATEGORIES[raw], "SEARCH_TERMS"
        else:
            hits = set()
            for prot in fam2prot.get(raw, set()):
                hits |= comp2sys.get(prot.lower(), set())
            if not hits:
                hits |= comp2sys.get(raw.lower(), set())
            if hits:
                system, how = sorted(hits, key=_system_rank)[0], "PADLOC yaml"
            else:
                cand = [v for k, v in sys_meta.items()
                        if k.lower().replace("_", "") ==
                           raw.lower().replace("_", "")]
                if cand:
                    system, how = cand[0], "sys_meta name"

        # LABEL_OVERRIDE first (see config 7b-ter), then the normal chain:
        # system name -> protein names -> raw label (see module docstring)
        prots = sorted(fam2prot.get(raw, set()))
        cat, key, matched = override_category(raw)
        if cat is not None:
            how = "LABEL_OVERRIDE"
            system = system or matched
        else:
            cat, key, matched = to_category(system, *prots, raw)

        if cat is None:
            rows.append((raw, n, system or "-", "UNMAPPED", how or "-"))
            unmapped.append((raw, n, ";".join(prots) if prots else "-"))
        else:
            via = f"{how or 'label'} | '{matched}' -> '{key}'"
            rows.append((raw, n, system or matched, cat, via))

    with open(F_OUT, "w") as out:
        out.write("category_raw\tn_genes\tsystem\tcategory\tsource\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # --- report -----------------------------------------------------------
    total = sum(seen.values())
    mapped = sum(n for _, n, _, c, _ in rows if c != "UNMAPPED")
    print("\n" + "=" * 78)
    print(" COVERAGE")
    print("=" * 78)
    print(f"   genes resolved : {mapped} / {total} ({100.0*mapped/total:.1f}%)")
    print(f"   labels left    : {len(unmapped)}")

    agg = Counter()
    for _, n, _, c, _ in rows:
        agg[c] += n
    print("\n FINAL CATEGORIES (all genes, before the catalogue filter):")
    for c in DEFENCE_CATEGORIES + ["EXCLUDED (not defence)", "UNMAPPED"]:
        if agg.get(c):
            print(f"   {agg[c]:>6}  {c}   ({100.0*agg[c]/total:.1f}%)")

    print("\n how the 25 largest labels resolved:")
    print(f"   {'n':>6}  {'raw':<26} {'category':<34} via")
    for raw, n, system, cat, how in rows[:25]:
        print(f"   {n:>6}  {raw[:25]:<26} {cat[:33]:<34} {how}")

    if unmapped:
        print("\n" + "-" * 78)
        print(" STILL UNMAPPED — add a prefix key to SYSTEM_TO_CATEGORY")
        print("-" * 78)
        print(f"   {'n':>6}  {'label':<30} protein names from hmm_meta")
        for raw, n, prots in unmapped[:40]:
            print(f"   {n:>6}  {raw[:29]:<30} {prots[:40]}")
        if len(unmapped) > 40:
            print(f"   ... and {len(unmapped)-40} more (see the output file)")

    print("\n" + "=" * 78)
    print(f" Output: {F_OUT}")


if __name__ == "__main__":
    main()
