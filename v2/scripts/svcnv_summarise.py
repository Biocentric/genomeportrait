#!/usr/bin/env python3
"""svcnv_summarise.py — structural/copy-number variants: counts plus the largest calls.

Counts alone ("3593 DEL") say nothing about the genome, so this also emits a table of the
biggest PASS events with coordinates and size, which is what you can actually look up.
"""
import argparse
import gzip
import os
import sys

import pandas as pd

TOP_N = 15
# Germline SVs are typically well under a megabase; larger calls warrant a health warning.
ARTIFACT_BP = 5_000_000
# BNDs are single breakend records with no span, so they can't be size-ranked.
SIZED_TYPES = ("DEL", "DUP", "INS", "INV")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def parse_sv(path):
    """Return (counts_by_type, list of sized PASS records)."""
    counts, recs = {}, []
    if not path:
        return counts, recs
    try:
        with opener(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                chrom, pos, vid, flt, info = f[0], f[1], f[2], f[6], f[7]
                kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
                svt = kv.get("SVTYPE", "OTHER")
                counts[svt] = counts.get(svt, 0) + 1
                if svt not in SIZED_TYPES:
                    continue
                try:
                    svlen = abs(int(kv["SVLEN"])) if "SVLEN" in kv else (
                        int(kv["END"]) - int(pos) if "END" in kv else None)
                except (ValueError, TypeError):
                    svlen = None
                if not svlen:
                    continue
                recs.append({
                    "type": svt,
                    "location": f"{chrom}:{int(pos):,}",
                    "size_bp": svlen,
                    "filter": flt or ".",
                    "id": vid or ".",
                })
    except Exception as e:
        print(f"[svcnv] SV parse: {e}", file=sys.stderr)
    return counts, recs


def cnv_calls(path):
    """Parse a CNVpytor call TSV (no header): col0 = type, col1 = region, col2 = size."""
    gains = losses = 0
    recs = []
    if not path:
        return gains, losses, recs
    try:
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if not p or not p[0]:
                    continue
                t = p[0].strip().lower()
                if t.startswith("dup"):
                    gains += 1
                elif t.startswith("del"):
                    losses += 1
                else:
                    continue
                try:
                    size = int(float(p[2])) if len(p) > 2 else None
                except ValueError:
                    size = None
                recs.append({
                    "type": "duplication" if t.startswith("dup") else "deletion",
                    "location": p[1] if len(p) > 1 else ".",
                    "size_bp": size or 0,
                    "filter": "CNVpytor",
                    "id": ".",
                })
    except Exception as e:
        print(f"[svcnv] CNV parse: {e}", file=sys.stderr)
    return gains, losses, recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sv", default=None)
    ap.add_argument("--cnr", default=None)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    sv_counts, sv_recs = parse_sv(a.sv)
    dup, dele, cnv_recs = cnv_calls(a.cnr)

    rows = [{"category": "Structural variant", "type": svt, "count": n}
            for svt, n in sorted(sv_counts.items())]
    rows.append({"category": "Copy-number (CNVpytor)", "type": "duplications", "count": dup})
    rows.append({"category": "Copy-number (CNVpytor)", "type": "deletions", "count": dele})
    if not rows:
        rows.append({"category": "none", "type": "none", "count": 0})

    os.makedirs(a.outdir, exist_ok=True)
    counts = os.path.join(a.outdir, "1_counts.tsv")
    with open(counts, "w") as fh:
        fh.write(f"# structural / copy-number variant counts for {a.sample}\n")
    pd.DataFrame(rows).to_csv(counts, sep="\t", index=False, mode="a")

    # Largest events, PASS first — the part that is actually interpretable.
    allrecs = sv_recs + cnv_recs
    top = os.path.join(a.outdir, "2_largest_events.tsv")
    if allrecs:
        df = pd.DataFrame(allrecs)
        df["_pass"] = df["filter"].map(lambda v: 0 if str(v).upper() in ("PASS", "CNVPYTOR", ".") else 1)
        df = df.sort_values(["_pass", "size_bp"], ascending=[True, False]).drop(columns="_pass")
        df = df.head(TOP_N)
        # Short-read SV callers routinely emit multi-megabase events that are mapping or
        # assembly artifacts rather than real germline variation, and sorting by size puts
        # them straight at the top — so say so rather than presenting them as findings.
        df["note"] = df["size_bp"].map(
            lambda v: "implausibly large for a germline event — usually a short-read calling "
                      "artifact" if v >= ARTIFACT_BP else ".")
        df["size_bp"] = df["size_bp"].map(lambda v: f"{int(v):,}")
        df = df[["type", "location", "size_bp", "filter", "note", "id"]]
    else:
        df = pd.DataFrame([{"type": ".", "location": ".", "size_bp": ".",
                            "filter": "no sized SV/CNV calls", "note": ".", "id": "."}])
    with open(top, "w") as fh:
        fh.write(f"# largest structural / copy-number changes for {a.sample}\n")
    df.to_csv(top, sep="\t", index=False, mode="a")

    print(f"[svcnv_summarise] {sum(sv_counts.values())} SVs, {dup+dele} CNVs, "
          f"{len(allrecs)} sized -> top {len(df)}", file=sys.stderr)


if __name__ == "__main__":
    main()
