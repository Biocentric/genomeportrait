#!/usr/bin/env python3
"""Download harmonised GRCh38 scoring files from the PGS Catalog.

    python3 pgs_fetch.py --panel panel.tsv [--extra PGS000001,PGS000002] --outdir scorefiles

The panel is a TSV with a `pgs_id` column (comment lines start with '#'). Files that
are already present and gzip-complete are kept, so the step can be re-run after an
interrupted download. A score with no harmonised GRCh38 file is reported and skipped;
it is not an error, because one missing score must not stop a 200-score panel.

Writes <outdir>/<PGS>_hmPOS_GRCh38.txt.gz and <outdir>/fetch_report.tsv.
"""
import argparse
import concurrent.futures as cf
import csv
import gzip
import os
import sys
import time
import urllib.error
import urllib.request

URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{id}/ScoringFiles/"
       "Harmonized/{id}_hmPOS_GRCh38.txt.gz")
UA = {"User-Agent": "genomeportrait-pgs-fetch"}


def read_panel(path):
    ids = []
    with open(path, encoding="utf-8") as fh:
        rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in rows:
            if r.get("pgs_id"):
                ids.append(r["pgs_id"].strip())
    return ids


def complete(path):
    """A truncated download is a valid file until someone reads the end of it."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return os.path.getsize(path) > 0
    except (OSError, EOFError):
        return False


def fetch(pgs_id, outdir, tries=6):
    dest = os.path.join(outdir, f"{pgs_id}_hmPOS_GRCh38.txt.gz")
    if os.path.exists(dest) and complete(dest):
        return pgs_id, "cached", os.path.getsize(dest)
    url, err = URL.format(id=pgs_id), ""
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as resp, open(dest + ".part", "wb") as out:
                while True:
                    buf = resp.read(1 << 20)
                    if not buf:
                        break
                    out.write(buf)
            if complete(dest + ".part"):
                os.replace(dest + ".part", dest)
                return pgs_id, "downloaded", os.path.getsize(dest)
            err = "truncated"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return pgs_id, "no_harmonised_file", 0
            err = f"HTTP {e.code}"
        except Exception as e:  # network resets, timeouts
            err = type(e).__name__
        time.sleep(5 * (i + 1))
    if os.path.exists(dest + ".part"):
        os.remove(dest + ".part")
    return pgs_id, f"failed ({err})", 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel")
    ap.add_argument("--extra", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()

    ids = read_panel(a.panel) if a.panel else []
    ids += [x.strip() for x in a.extra.split(",") if x.strip()]
    ids = sorted(dict.fromkeys(ids))
    if not ids:
        sys.exit("no PGS IDs given")
    os.makedirs(a.outdir, exist_ok=True)

    results = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for n, res in enumerate(ex.map(lambda i: fetch(i, a.outdir), ids), 1):
            results.append(res)
            print(f"[{n}/{len(ids)}] {res[0]} {res[1]} {res[2] / 1e6:.1f} MB", file=sys.stderr, flush=True)

    with open(os.path.join(a.outdir, "fetch_report.tsv"), "w") as fh:
        fh.write("pgs_id\tstatus\tbytes\n")
        for r in results:
            fh.write("\t".join(map(str, r)) + "\n")
    failed = [r for r in results if r[1].startswith("failed")]
    ok = [r for r in results if r[1] in ("cached", "downloaded")]
    print(f"{len(ok)} of {len(ids)} scoring files present; {len(failed)} failed", file=sys.stderr)
    if failed:
        sys.exit("downloads failed: " + ", ".join(r[0] for r in failed))


if __name__ == "__main__":
    main()
