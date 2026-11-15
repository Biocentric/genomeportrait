#!/usr/bin/env python3
"""Download named files from a Zenodo record and check each against its published MD5.

    python3 zenodo_fetch.py --record 11212339 --out . hg38_mutationrate.bed hg38_ancestral.tar.gz
    python3 zenodo_fetch.py --record 13368126 --out archaic --pattern 'individuals_highcov\\.\\d+\\.bcf'

Needs only the Python standard library (it runs in images without curl or wget). A file
that already exists with the right checksum is not downloaded again; a partial or corrupt
one is replaced. Any mismatch after download is an error.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time
import urllib.request

API = "https://zenodo.org/api/records/{}"


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, dest, tries=4):
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(dest + ".part", "wb") as out:
                shutil.copyfileobj(r, out, 1 << 22)
            os.replace(dest + ".part", dest)
            return
        except Exception as e:  # network hiccups: retry with a pause
            print(f"[zenodo_fetch] {os.path.basename(dest)}: attempt {attempt} failed ({e})", file=sys.stderr)
            if attempt == tries:
                raise
            time.sleep(15 * attempt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True)
    ap.add_argument("--out", default=".")
    ap.add_argument("--pattern", help="regular expression: fetch every file whose name matches")
    ap.add_argument("files", nargs="*")
    a = ap.parse_args()

    with urllib.request.urlopen(API.format(a.record), timeout=60) as r:
        rec = json.load(r)
    listed = {f["key"]: f for f in rec["files"]}
    want = list(a.files)
    if a.pattern:
        want += sorted(k for k in listed if re.fullmatch(a.pattern, k))
    missing = [k for k in want if k not in listed]
    if missing or not want:
        sys.exit(f"[zenodo_fetch] record {a.record}: not found {missing or 'any matching file'}")

    os.makedirs(a.out, exist_ok=True)
    for k in want:
        algo, expected = listed[k]["checksum"].split(":", 1)
        assert algo == "md5", listed[k]["checksum"]
        dest = os.path.join(a.out, k)
        if os.path.exists(dest) and md5sum(dest) == expected:
            continue
        url = f"https://zenodo.org/records/{a.record}/files/{k}?download=1"
        print(f"[zenodo_fetch] {k} ({listed[k]['size'] / 1e6:.0f} MB)", file=sys.stderr)
        fetch(url, dest)
        got = md5sum(dest)
        if got != expected:
            os.remove(dest)
            sys.exit(f"[zenodo_fetch] {k}: MD5 {got}, Zenodo lists {expected}")


if __name__ == "__main__":
    main()
