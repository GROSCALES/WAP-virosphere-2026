#!/usr/bin/env python3
"""
01_build_gene_table.py — Build LEVEL 2 of the pipeline: one row per ORF.

The spine is the prodigal gene call over the full vOTU catalogue
(393,330 ORFs, ONE coordinate frame). VIBRANT annotations are layered on
top as columns, never merged into a string.

Two things this fixes for good:

  1. The 18 VIBRANT columns are kept apart. Collapsing them with
     " ".join(cols[2:]) was the root cause of the false-friend collisions
     (VOG reference names bleeding into identity), of get_evalue() always
     returning the KO e-value, and of detecting the AMG flag by regex.
     The AMG flag is column 3. It is read as a field.

  2. VIBRANT renames ORFs on excised prophages to <scaffold>_fragment_F_N,
     which never matched <scaffold>_N and silently dropped 1,290 ORFs on
     21 contigs. All six vOTUs that were being patched by hand in
     06a/06b are on that list.

     N is ALREADY the ORF number in the parent contig, not an index within
     the fragment. Verified against the prophage coordinate file: for every
     one of the 27 fragments, min(N) == 'protein start', max(N) ==
     'protein stop', and count == stop - start + 1. So the translation is
     simply to drop the _fragment_F infix:

         <scaffold>_fragment_2_64   ->   <scaffold>_64

     Adding 'protein start' as an offset (an earlier attempt) shifted every
     ORF on fragments whose start was > 1. Those shifted IDs mostly still
     EXISTED in the spine, so they passed an existence check while pointing
     at the wrong gene — a silent mis-assignment. Hence the assertion in
     load_fragment_offsets(): the mapping is only sound if the fragment ORF
     numbers exactly fill the declared parent range.

Output: 01_gene_table_L2.tsv
"""

import os
import re
import sys
from collections import defaultdict

from config import (
    F_PRODIGAL_GFF, F_VIB_ANNOT, F_VIB_PROPHAGE, F_GENE_TABLE_L2,
    VIB_ANNOT_COLS, VIB_PROPHAGE_COLS, check_header, require,
)

GFF_ID_RE = re.compile(r"ID=([^;]+)")


# ==========================================================================
# STEP 1 — spine: every ORF with its coordinates
# ==========================================================================
def load_spine(path):
    """Read the prodigal GFF into an ordered per-ORF record list.

    The GFF was rewritten by 08_defensas_votus.sh so that ID=<contig>_<n>,
    matching the .faa headers and the IDs used by PADLOC and DefenseFinder.
    """
    print(f"[1] Spine from {os.path.basename(path)}")
    genes = []
    seen = set()
    bad_id = 0

    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            m = GFF_ID_RE.search(c[8])
            if not m:
                bad_id += 1
                continue
            pid = m.group(1)
            contig = c[0]
            # ID must be <contig>_<orf_number>; the contig name itself
            # contains underscores, so split on the LAST one only.
            if not pid.startswith(contig + "_"):
                bad_id += 1
                continue
            try:
                orf_num = int(pid.rsplit("_", 1)[1])
            except (IndexError, ValueError):
                bad_id += 1
                continue
            if pid in seen:
                continue
            seen.add(pid)
            genes.append({
                "protein_id": pid,
                "votu_id": contig,
                "orf_num": orf_num,
                "start": int(c[3]),
                "end": int(c[4]),
                "strand": c[6],
            })

    print(f"    ORFs            : {len(genes)}")
    print(f"    contigs         : {len(set(g['votu_id'] for g in genes))}")
    if bad_id:
        print(f"    [warn] unparsable ID lines skipped: {bad_id}")
    if not genes:
        raise SystemExit("[FATAL] empty spine — check the GFF path")
    return genes


# ==========================================================================
# STEP 2 — fragment -> parent ID translation
# ==========================================================================
def load_fragment_offsets(path):
    """Read the declared parent ORF range of each excised-prophage fragment.

    Used only to VERIFY the identity mapping, not to shift it. See the
    module docstring: the ORF number in a fragment ID is already the parent
    ORF number.
    """
    print(f"[2] Prophage fragments from {os.path.basename(path)}")
    check_header(path, VIB_PROPHAGE_COLS)
    ranges = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            try:
                ranges[c[1]] = (c[0], int(c[2]), int(c[3]))
            except ValueError:
                continue
    print(f"    fragments       : {len(ranges)}")
    print(f"    parent contigs  : {len(set(v[0] for v in ranges.values()))}")
    return ranges


FRAG_RE = re.compile(r"^(?P<scaffold>.+)_fragment_\d+_(?P<n>\d+)$")


def translate_id(pid):
    """Return (parent_id, was_fragment). Non-fragment IDs pass through.

    Drops the _fragment_F infix. No offset: N is already the parent ORF
    number.
    """
    m = FRAG_RE.match(pid)
    if not m:
        return pid, False
    return f"{m.group('scaffold')}_{m.group('n')}", True


def assert_fragment_ranges(observed, ranges):
    """The identity mapping is only sound if, for every fragment, the ORF
    numbers seen in the annotations exactly fill [protein start, protein stop].
    Anything else means the numbering is not parent-relative and every
    fragment gene would be silently mis-assigned."""
    bad = []
    for frag, nums in observed.items():
        if frag not in ranges:
            bad.append(f"{frag}: not in coordinates file")
            continue
        _, pstart, pstop = ranges[frag]
        if min(nums) != pstart or max(nums) != pstop:
            bad.append(f"{frag}: ORFs {min(nums)}-{max(nums)} "
                       f"vs declared {pstart}-{pstop}")
    if bad:
        raise SystemExit(
            "[FATAL] fragment ORF numbering is NOT parent-relative:\n  "
            + "\n  ".join(bad)
            + "\n  Do not trust any fragment gene until this is resolved.")
    print(f"    [check] all {len(observed)} fragments fill their declared "
          f"parent ORF range exactly -> identity mapping is sound")


# ==========================================================================
# STEP 3 — VIBRANT annotations, 18 columns kept apart
# ==========================================================================
def load_vibrant(path, ranges, valid_ids):
    print(f"[3] VIBRANT annotations from {os.path.basename(path)}")
    check_header(path, VIB_ANNOT_COLS)

    n_cols = len(VIB_ANNOT_COLS)
    ann = {}
    stats = defaultdict(int)
    observed = defaultdict(list)      # fragment -> ORF numbers seen

    # first pass: collect fragment ORF numbers so the mapping can be
    # verified BEFORE any gene is assigned
    with open(path) as fh:
        next(fh)
        for line in fh:
            pid = line.split("\t", 1)[0]
            m = FRAG_RE.match(pid)
            if m:
                frag = pid.rsplit("_", 1)[0]
                observed[frag].append(int(m.group("n")))
    if observed:
        assert_fragment_ranges(observed, ranges)

    with open(path) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2:
                continue
            # trailing empty fields may be stripped: pad, never index blindly
            if len(c) < n_cols:
                c += [""] * (n_cols - len(c))

            raw_id = c[0]
            pid, was_frag = translate_id(raw_id)
            stats["rows"] += 1
            if was_frag:
                stats["fragment_rows"] += 1

            if pid not in valid_ids:
                stats["unmatched"] += 1
                if was_frag:
                    stats["fragment_unmatched"] += 1
                continue
            if was_frag:
                stats["fragment_rescued"] += 1

            ann[pid] = {
                "vib_protein_raw": raw_id,
                "ko": c[2].strip(),
                "is_amg": "1" if c[3].strip() == "AMG" else "0",
                "ko_name": c[4].strip(),
                "ko_eval": c[5].strip(),
                "ko_vscore": c[7].strip(),
                "pfam": c[8].strip(),
                "pfam_name": c[9].strip(),
                "pfam_eval": c[10].strip(),
                "pfam_vscore": c[12].strip(),
                "vog": c[13].strip(),
                "vog_name": c[14].strip(),
                "vog_eval": c[15].strip(),
                "vog_vscore": c[17].strip(),
                "from_fragment": "1" if was_frag else "0",
            }

    print(f"    rows read              : {stats['rows']}")
    print(f"    matched to spine       : {len(ann)}")
    print(f"    unmatched              : {stats['unmatched']}")
    print(f"    fragment rows          : {stats['fragment_rows']}")
    print(f"    fragment RESCUED       : {stats['fragment_rescued']}")
    print(f"    fragment still unmatched: {stats['fragment_unmatched']}")
    return ann, stats


# ==========================================================================
# MAIN
# ==========================================================================
def main():
    require(F_PRODIGAL_GFF, F_VIB_ANNOT, F_VIB_PROPHAGE)

    print("=" * 66)
    print(" LEVEL 2 — one row per ORF, VIBRANT layered as columns")
    print("=" * 66)

    genes = load_spine(F_PRODIGAL_GFF)
    valid_ids = {g["protein_id"] for g in genes}

    ranges = load_fragment_offsets(F_VIB_PROPHAGE)
    ann, stats = load_vibrant(F_VIB_ANNOT, ranges, valid_ids)

    vib_cols = ["vib_protein_raw", "ko", "is_amg", "ko_name", "ko_eval",
                "ko_vscore", "pfam", "pfam_name", "pfam_eval", "pfam_vscore",
                "vog", "vog_name", "vog_eval", "vog_vscore", "from_fragment"]
    base_cols = ["protein_id", "votu_id", "orf_num", "start", "end", "strand"]

    print(f"[4] Writing {os.path.basename(F_GENE_TABLE_L2)}")
    n_annotated = 0
    with open(F_GENE_TABLE_L2, "w") as out:
        out.write("\t".join(base_cols + ["has_vibrant"] + vib_cols) + "\n")
        for g in genes:
            a = ann.get(g["protein_id"])
            if a:
                n_annotated += 1
            row = [str(g[c]) for c in base_cols]
            row.append("1" if a else "0")
            row += [(a or {}).get(c, "") for c in vib_cols]
            out.write("\t".join(row) + "\n")

    # ---------------- verification ----------------
    n_amg = sum(1 for a in ann.values() if a["is_amg"] == "1")
    pct = 100.0 * n_annotated / len(genes)

    print("\n" + "=" * 66)
    print(" VERIFICATION")
    print("=" * 66)
    print(f" ORFs in spine                   : {len(genes)}")
    print(f" ORFs with VIBRANT annotation    : {n_annotated}  ({pct:.1f}%)")
    print(f" ORFs with NO VIBRANT annotation : {len(genes) - n_annotated}")
    print(f"   (VIBRANT only annotates scaffolds it classified as viral,")
    print(f"    so this gap is expected and must be reported, not hidden)")
    print(f" ORFs flagged AMG by VIBRANT     : {n_amg}")
    print("-" * 66)
    print(f" Fragment ORFs rescued           : {stats['fragment_rescued']}")
    print(f" Fragment ORFs still unmatched   : {stats['fragment_unmatched']}")
    if stats["fragment_unmatched"] == 0 and stats["fragment_rescued"] > 0:
        print("   -> fragment translation COMPLETE.")
        print("      The manual patches in 06a/06b are now redundant.")
    elif stats["fragment_unmatched"]:
        print("   -> [CHECK] some fragment ORFs did not resolve; inspect before")
        print("      assuming the manual patches can be dropped.")
    print("=" * 66)
    print(f" Output: {F_GENE_TABLE_L2}")


if __name__ == "__main__":
    main()
