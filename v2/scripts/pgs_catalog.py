#!/usr/bin/env python3
"""Build or extend the PGS Catalog metadata snapshot used for interpretation.

    python3 pgs_catalog.py --ids PGS000001,PGS000002 --outdir catalog [--seed old_catalog_dir]
    python3 pgs_catalog.py --panel panel.tsv --outdir catalog --seed assets/pgs/catalog

Writes <outdir>/scores.json and <outdir>/performance.json, keeping only the fields
the interpretation step reads. Scores already present in --seed are copied, not
re-fetched, so adding one score to a 200-score panel costs one request, not 400.

The REST API answers HTTP 429 above roughly one request per second; requests are
spaced and retried with the server's Retry-After.
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request

API = "https://www.pgscatalog.org/rest"
UA = {"User-Agent": "genomeportrait-pgs-catalog"}


def get(url):
    for i in range(8):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            wait = int(e.headers.get("Retry-After") or 0) or 10 * (i + 1)
            print(f"  HTTP {e.code}, waiting {wait}s", file=sys.stderr, flush=True)
            time.sleep(wait)
        except Exception as e:
            print(f"  {type(e).__name__}, retrying", file=sys.stderr, flush=True)
            time.sleep(5 * (i + 1))
    raise SystemExit("gave up: " + url)


def slim_score(s):
    pub = s.get("publication") or {}
    anc = s.get("ancestry_distribution") or {}
    return {
        "name": s.get("name"),
        "trait_reported": s.get("trait_reported"),
        "trait_efo": [{"id": t.get("id"), "label": t.get("label")} for t in s.get("trait_efo") or []],
        "variants_number": s.get("variants_number"),
        "weight_type": s.get("weight_type"),
        "license": s.get("license"),
        "publication": {k: pub.get(k) for k in ("firstauthor", "date_publication", "journal", "doi", "PMID")},
        "ancestry_distribution": {k: {"dist": (anc.get(k) or {}).get("dist")} for k in ("gwas", "dev", "eval")
                                  if anc.get(k)},
    }


def slim_perf(r):
    pm = r.get("performance_metrics") or {}
    samples = (r.get("sampleset") or {}).get("samples") or []

    def metrics(key):
        return [{"name_short": m.get("name_short"), "estimate": m.get("estimate")} for m in pm.get(key) or []]

    return {
        "phenotyping_reported": r.get("phenotyping_reported"),
        "comments": r.get("performance_comments"),
        "publication": ((r.get("publication") or {}).get("firstauthor") or "") + " " +
                       ((r.get("publication") or {}).get("date_publication") or "")[:4],
        "samples": [{"ancestry_broad": s.get("ancestry_broad"), "n": s.get("sample_number")} for s in samples],
        "effect_sizes": metrics("effect_sizes"),
        "class_acc": metrics("class_acc"),
        "othermetrics": metrics("othermetrics"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", default="")
    ap.add_argument("--panel")
    ap.add_argument("--seed", help="existing snapshot directory to reuse")
    ap.add_argument("--raw-scores", help="unslimmed scores JSON from an earlier fetch (dict by PGS id)")
    ap.add_argument("--raw-perf", help="unslimmed performance JSON from an earlier fetch")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    ids = [x.strip() for x in a.ids.split(",") if x.strip()]
    if a.panel:
        with open(a.panel, encoding="utf-8") as fh:
            ids += [r["pgs_id"] for r in csv.DictReader((l for l in fh if l[0] != "#"), delimiter="\t")]
    ids = sorted(dict.fromkeys(ids))

    scores, perf = {}, {}
    if a.seed and os.path.exists(os.path.join(a.seed, "scores.json")):
        scores = json.load(open(os.path.join(a.seed, "scores.json"), encoding="utf-8"))
        perf = json.load(open(os.path.join(a.seed, "performance.json"), encoding="utf-8"))
    if a.raw_scores:
        for k, v in json.load(open(a.raw_scores, encoding="utf-8")).items():
            scores.setdefault(k, slim_score(v))
    if a.raw_perf:
        for k, v in json.load(open(a.raw_perf, encoding="utf-8")).items():
            perf.setdefault(k, [slim_perf(r) for r in v])

    need = [i for i in ids if i not in scores]
    for k in range(0, len(need), 50):
        d = get(f"{API}/score/search?pgs_ids={','.join(need[k:k + 50])}&limit=50")
        for s in d["results"]:
            scores[s["id"]] = slim_score(s)
        time.sleep(1)
    for n, i in enumerate([i for i in ids if i not in perf], 1):
        res, url = [], f"{API}/performance/search?pgs_id={i}&limit=250"
        while url:
            d = get(url)
            res += d["results"]
            url = d.get("next")
            time.sleep(0.8)
        perf[i] = [slim_perf(r) for r in res]
        print(f"  [{n}] {i}: {len(res)} performance records", file=sys.stderr, flush=True)

    missing = [i for i in ids if i not in scores]
    if missing:
        print("not in the PGS Catalog: " + ", ".join(missing), file=sys.stderr)
    os.makedirs(a.outdir, exist_ok=True)
    keep = set(ids)
    json.dump({k: v for k, v in sorted(scores.items()) if k in keep},
              open(os.path.join(a.outdir, "scores.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=0)
    json.dump({k: v for k, v in sorted(perf.items()) if k in keep},
              open(os.path.join(a.outdir, "performance.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=0)
    print(f"snapshot: {len(keep & set(scores))} scores, "
          f"{sum(len(v) for k, v in perf.items() if k in keep)} performance records -> {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
