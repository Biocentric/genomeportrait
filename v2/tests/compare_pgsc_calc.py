#!/usr/bin/env python3
"""Compare the native PGS engine with a pgsc_calc run on the same genome.

    python3 compare_pgsc_calc.py <sample>.pgs_scores.tsv <pgsc_calc>/score/<sampleset>_pgs.txt.gz

Prints, for every quantity both report, the correlation and the size of the differences
across scores, then the scores that differ most. Exits 1 if Z_norm2 correlates below 0.98.
"""
import csv
import gzip
import math
import sys

import numpy as np

mine_path, pgsc_path = sys.argv[1:3]
with open(mine_path) as fh:
    mine = {r["pgs_id"]: r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}
ref = {}
with gzip.open(pgsc_path, "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["sampleset"] != "reference":
            ref[r["PGS"].split("_hmPOS")[0]] = r

ok = {k for k, m in mine.items() if m["status"] == "ok"}
both = sorted(ok & set(ref))
print(f"pgsc_calc reports {len(ref)} scores, the native engine {len(ok)}; {len(both)} in both")
only_native = sorted(ok - set(ref))
only_pgsc = sorted(set(ref) - ok)
if only_native:
    print("  native only:", ", ".join(only_native))
if only_pgsc:
    print("  pgsc_calc only:", ", ".join(f"{k} ({mine[k]['status']}, overlap {mine[k].get('overlap', '?')})"
                                         if k in mine else k for k in only_pgsc))


def col(d, c):
    return np.array([float(d[k][c]) for k in both])


print(f"\n{'quantity':28s} {'r':>8s} {'mean|d|':>9s} {'max|d|':>8s} {'rms':>8s}")
r_z2 = None
for c in ("Z_norm2", "Z_norm1", "Z_MostSimilarPop", "percentile_MostSimilarPop", "SUM"):
    x, y = col(mine, c), col(ref, c)
    d = x - y
    r = np.corrcoef(x, y)[0, 1]
    if c == "Z_norm2":
        r_z2 = r
    print(f"{c:28s} {r:8.5f} {np.abs(d).mean():9.4f} {np.abs(d).max():8.4f} {math.sqrt((d ** 2).mean()):8.4f}")

x, y = col(mine, "Z_norm2"), col(ref, "Z_norm2")
print("\nlargest Z_norm2 differences:")
for d, k, a, b in sorted(zip(np.abs(x - y), both, x, y), reverse=True)[:10]:
    print(f"  {k}  native {a:+.3f}  pgsc_calc {b:+.3f}  diff {d:.3f}  overlap {mine[k]['overlap']}")

band = lambda z: 0 if abs(z) < 0.5 else (1 if z > 0 else -1)
moved = [k for k, a, b in zip(both, x, y) if band(a) != band(b)]
print(f"\nscores crossing the |Z| = 0.5 'around average' boundary: {len(moved)} of {len(both)}"
      + (f" ({', '.join(moved)})" if moved else ""))
sys.exit(0 if r_z2 is not None and r_z2 >= 0.98 else 1)
