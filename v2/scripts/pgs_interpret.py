#!/usr/bin/env python3
"""Turn ancestry-adjusted polygenic scores into direction, verdict and effect size.

    python3 pgs_interpret.py --scores S.pgs_scores.tsv --panel panel.tsv --concepts concepts.tsv \
        --catalog catalog_dir --r2-exclusions r2_exclusions.tsv --sample S \
        --out-tsv S.pgs_interpreted.tsv --out-json S.pgs_interpreted.json

For every score:
  direction  what a HIGHER score predicts. Scores flagged `reversed` in the panel are
             flipped, so a positive Z always means "more of the trait as named".
  effect     quantitative traits: expected shift = sqrt(R2) * Z, with R2 from the score's
             published (preferably European) validation, converted to units where the
             concept table gives a population SD.
             binary traits: relative odds versus an average person = OR_per_SD ** Z.
  consensus  per trait concept, the median aligned Z over European-trained scores; the
             headline score is the best-validated one.

The rules are those of prs-selfscore/interpret/interpret.py; the trait knowledge that
was hard-coded there lives in the panel and concept tables.
"""
import argparse
import csv
import datetime
import json
import math
import os
import re
import statistics
import sys

sys.stdout.reconfigure(encoding="utf-8")

R2_PRIORITY = [  # (regex on the metric's short name, transform to R2)
    (r"^pgs r2\s+\(no covariates\)$", lambda v: v),
    (r"^(partial-r2|incremental r2 \(full-covars\)|incremental r2 \(full model vs\. covariates alone\)|"
     r"incremental variance explained by the prs|∆ r\^2|Δr-squared.*|adjusted r2 only by the pgs)$", lambda v: v),
    (r"^Δr2 \(%\)$", lambda v: v / 100),
    (r"^(partial correlation \(partial-r\)|pearson correlation coefficient|spearman's ρ)$", lambda v: v * v),
    (r"^r²$", lambda v: v),
]
OTHER_DOMAIN = "Other scores"


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def header_meta(path):
    """key=value pairs from the leading comment line of the engine's score table."""
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
    return dict(kv.split("=", 1) for kv in first.lstrip("# ").split() if "=" in kv) if first.startswith("#") else {}


def eur_trained(s):
    ad = s.get("ancestry_distribution") or {}
    dist = (ad.get("gwas") or ad.get("dev") or {}).get("dist") or {}
    label = ", ".join(f"{k} {v:g}%" for k, v in dist.items()) if dist else "unknown"
    return (not dist) or dist.get("EUR", 0) >= 50 or dist.get("MAE", 0) >= 50, label


def is_eur(rec):
    return bool(rec["samples"]) and all("europ" in (s.get("ancestry_broad") or "").lower() for s in rec["samples"])


def metrics(recs, kind, r2_max, pub, unreliable):
    eur = [r for r in recs if is_eur(r)]
    pool, anc = (eur, "EUR") if eur else (recs, "non-EUR")
    r2 = None
    for pat, f in R2_PRIORITY:
        vals = []
        for r in pool:
            for m in r["othermetrics"]:
                if m["estimate"] is not None and re.match(pat, (m["name_short"] or "").strip(), re.I):
                    v = f(m["estimate"])
                    if "percentage" in (r.get("comments") or "").lower():
                        v = v / 100
                    v = v / 100 if v > 1 else v
                    if 0 < v < 1:
                        vals.append(v)
        if vals:
            r2 = statistics.median(vals)
            break
    note = ""
    if r2 is not None and pub in unreliable:
        r2, note = None, f"R2 not used: {unreliable[pub]}"
    elif r2 is not None and kind == "q" and r2_max is not None and r2 > r2_max:
        r2, note = None, f"R2 {r2:.2f} exceeds the trait's heritability ceiling; not used"
    ors, aucs = [], []
    if kind == "b":
        for r in pool:
            for m in r["effect_sizes"]:
                # ORs outside 1/3-3 are not per-SD (top-vs-bottom decile contrasts reach 16-110)
                if m["name_short"] in ("OR", "HR") and m["estimate"] and 1 / 3 <= m["estimate"] <= 3:
                    ors.append(m["estimate"])
            for m in r["class_acc"]:
                if m["name_short"] in ("AUROC", "C-index") and m["estimate"]:
                    aucs.append(m["estimate"])
    return {"r2": r2, "r2_note": note, "or": statistics.median(ors) if ors else None,
            "auc": statistics.median(aucs) if aucs else None,
            "anc": anc if pool else "none", "n_eval": len(pool)}


def units(concept, shift_sd):
    unit, sd = concept["unit"], concept["population_sd"]
    if shift_sd is None or not unit or sd is None:
        return ""
    if unit == "log":
        return f"{(math.exp(shift_sd * sd) - 1) * 100:+.0f}%"
    v = shift_sd * sd
    s = f"{v:+.2f} {unit}" if abs(v) < 10 else f"{v:+.1f} {unit}"
    if unit == "mmol/L" and concept["key"] in ("ldl", "hdl", "tc"):
        s += f" ({v * 38.67:+.0f} mg/dL)"
    return s


def verdict(good, z):
    if abs(z) < 0.5 or good == 0:
        return "neutral"
    return "favourable" if (good > 0) == (z > 0) else "unfavourable"


def strength(z):
    a = abs(z)
    return "around average" if a < 0.5 else "somewhat" if a < 1 else "clearly" if a < 2 else "strongly"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--concepts", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--r2-exclusions")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--percentile", choices=["auto", "population", "adjusted"], default="auto",
                    help="population: rank within the most similar reference population; "
                         "adjusted: from the ancestry-adjusted Z; auto: population when the "
                         "similarity is >= 0.9, adjusted otherwise")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    catalog = json.load(open(os.path.join(a.catalog, "scores.json"), encoding="utf-8"))
    perf = json.load(open(os.path.join(a.catalog, "performance.json"), encoding="utf-8"))
    unreliable = {r["publication"]: r["reason"] for r in read_tsv(a.r2_exclusions)} if a.r2_exclusions else {}

    concepts, dom_order = {}, []
    for r in read_tsv(a.concepts):
        concepts[r["concept"]] = dict(
            key=r["concept"], domain=r["domain"], title=r["title"], higher=r["higher_means"],
            lower=r["lower_means"], good=int(r["good"]), kind=r["kind"], unit=r["unit"] or None,
            population_sd=float(r["population_sd"]) if r["population_sd"] else None,
            r2_max=float(r["r2_max"]) if r["r2_max"] else None, note=r["note"])
        if r["domain"] not in dom_order:
            dom_order.append(r["domain"])
    panel = {r["pgs_id"]: r for r in read_tsv(a.panel)}
    for pid, r in panel.items():
        if r["concept"] not in concepts:
            sys.exit(f"panel: {pid} refers to unknown concept '{r['concept']}'")

    meta = header_meta(a.scores)
    conf = float(meta.get("confidence", 0) or 0)
    basis = a.percentile if a.percentile != "auto" else ("population" if conf >= 0.9 else "adjusted")

    recs, excluded = [], []
    for r in read_tsv(a.scores):
        pid = r["pgs_id"]
        s = catalog.get(pid, {})
        trait = s.get("trait_reported") or (panel.get(pid) or {}).get("trait_reported") or pid
        if r["status"] != "ok":
            excluded.append(dict(acc=pid, trait=trait, status=r["status"],
                                 overlap=float(r["overlap"]) if r.get("overlap") else None))
            continue
        p = panel.get(pid)
        if p:
            key = p["concept"]
        else:
            # a score requested outside the panel: report it, but claim nothing about direction
            key = "x_" + pid
            concepts[key] = dict(key=key, domain=OTHER_DOMAIN, title=trait,
                                 higher="higher score, as coded by its authors",
                                 lower="lower score, as coded by its authors", good=0, kind="q",
                                 unit=None, population_sd=None, r2_max=None,
                                 note="Not part of the curated panel: direction and health meaning are not assigned.")
            if OTHER_DOMAIN not in dom_order:
                dom_order.append(OTHER_DOMAIN)
        c = concepts[key]
        sign = -1 if p and p["reversed"] == "1" else 1
        z_raw = float(r["Z_norm2"])
        pct_raw = float(r["percentile_MostSimilarPop"] if basis == "population" else r["percentile_norm2"])
        eur, gwas_anc = eur_trained(s)
        pubd = s.get("publication") or {}
        pub = f'{pubd.get("firstauthor") or "?"} {(pubd.get("date_publication") or "")[:4]}'.strip()
        m = metrics(perf.get(pid, []), c["kind"], c["r2_max"], pub, unreliable)
        z = z_raw * sign
        shift_sd = math.sqrt(m["r2"]) * z if (m["r2"] and c["kind"] == "q") else None
        rel_odds = m["or"] ** z_raw if m["or"] else None      # the OR is in the score's own coding
        if rel_odds and sign == -1:
            rel_odds = 1 / rel_odds
        recs.append(dict(acc=pid, key=key, trait=trait, pub=pub, n_var=int(r["n_variants"]),
                         match=round(100 * float(r["overlap"]), 1), z=z,
                         pct=pct_raw if sign == 1 else 100 - pct_raw, sign=sign,
                         verified=bool(p and p["direction_verified"] == "1"), eur=eur, gwas_anc=gwas_anc,
                         shift_sd=shift_sd, rel_odds=rel_odds, z_pop=float(r["Z_MostSimilarPop"]) * sign, **m))

    groups = {}
    for x in recs:
        groups.setdefault(x["key"], []).append(x)
    summary = []
    for key, xs in groups.items():
        c = concepts[key]
        use = [x for x in xs if x["eur"]] or xs
        zc = statistics.median(x["z"] for x in use)

        def quality(x):
            if c["kind"] == "q":
                return (x["r2"] is not None, x["anc"] == "EUR", x["r2"] or 0, x["n_var"])
            return (x["or"] is not None, x["anc"] == "EUR", abs(math.log(x["or"])) if x["or"] else 0, x["n_var"])

        zs = [x["z"] for x in use]
        summary.append(dict(key=key, dom=c["domain"], title=c["title"], higher=c["higher"], lower=c["lower"],
                            good=c["good"], kind=c["kind"], note=c["note"], n=len(xs), n_use=len(use),
                            zc=zc, zmin=min(zs), zmax=max(zs), pct_c=statistics.median(x["pct"] for x in use),
                            head=max(use, key=quality), verdict=verdict(c["good"], zc), strength=strength(zc)))
    summary.sort(key=lambda s: (dom_order.index(s["dom"]), -abs(s["zc"])))

    with open(a.out_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["domain", "trait_concept", "PGS", "trait_reported", "source", "higher_score_means",
                    "coding_flipped", "direction_verified", "aligned_percentile", "aligned_Z", "verdict",
                    "R2_eval", "expected_shift_SD", "expected_shift_units", "OR_or_HR_per_SD",
                    "relative_odds_vs_average", "AUC", "eval_ancestry", "gwas_ancestry", "european_trained",
                    "n_variants", "match_%"])
        for s in summary:
            c = concepts[s["key"]]
            for x in sorted(groups[s["key"]], key=lambda x: (not x["eur"], -(x["r2"] or 0))):
                w.writerow([s["dom"], s["title"], x["acc"], x["trait"], x["pub"], s["higher"],
                            "yes" if x["sign"] == -1 else "no", "yes" if x["verified"] else "",
                            f'{x["pct"]:.1f}', f'{x["z"]:+.2f}', verdict(s["good"], x["z"]),
                            f'{x["r2"]:.4f}' if x["r2"] else "",
                            f'{x["shift_sd"]:+.2f}' if x["shift_sd"] is not None else "",
                            units(c, x["shift_sd"]), f'{x["or"]:.2f}' if x["or"] else "",
                            f'{x["rel_odds"]:.2f}' if x["rel_odds"] else "",
                            f'{x["auc"]:.2f}' if x["auc"] else "", x["anc"], x["gwas_anc"],
                            "yes" if x["eur"] else "NO", x["n_var"], x["match"]])

    for s in summary:
        h = s["head"]
        c = concepts[s["key"]]
        eff = units(c, h["shift_sd"]) or (f'{h["shift_sd"]:+.2f} SD' if h["shift_sd"] is not None else "") \
            or (f'odds x{h["rel_odds"]:.2f}' if h["rel_odds"] else "n/a")
        print(f'{s["dom"][:18]:18s} {s["title"][:30]:30s} Zc={s["zc"]:+.2f} [{s["zmin"]:+.2f},{s["zmax"]:+.2f}] '
              f'n={s["n_use"]}/{s["n"]} {s["verdict"]:12s} head={h["acc"]} R2={h["r2"] or 0:.3f} -> {eff}')

    for x in recs:
        x["units"] = units(concepts[x["key"]], x["shift_sd"])
    out = dict(
        meta=dict(sample=a.sample, population=meta.get("most_similar_population", "?"), confidence=conf,
                  n_reference=int(meta.get("n_reference", 0) or 0), percentile_basis=basis,
                  n_scores=len(recs), n_traits=len(summary), n_excluded=len(excluded),
                  generated=datetime.date.today().isoformat(), domains=dom_order),
        summary=summary, scores=recs, excluded=excluded,
        concepts={k: v for k, v in concepts.items() if k in groups})
    json.dump(out, open(a.out_json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"[pgs_interpret] {len(recs)} scores in {len(summary)} traits; {len(excluded)} excluded", file=sys.stderr)


if __name__ == "__main__":
    main()
