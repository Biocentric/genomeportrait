#!/usr/bin/env python3
"""pharmcat_report.py — flatten the PharmCAT JSON into a gene/phenotype table.

Accepts either the reporter JSON (top-level "genes") or the phenotyper JSON
(top-level "geneReports"); both nest as {source: {gene: info}}.
"""
import argparse
import json
import sys

import pandas as pd

COLUMNS = ["gene", "what_this_gene_does", "diplotype", "phenotype",
           "what_the_phenotype_means", "activity", "source", "reference"]


# Plain-language reference for the genes PharmCAT reports and for the metaboliser calls.
# Descriptive pharmacology only — what the enzyme does and what the category means in
# general. Nothing here is specific to this person or to any prescribing decision.
GENE_FUNCTION = {
    "ABCG2":   "Transporter that pumps drugs out of cells; studied for rosuvastatin and allopurinol handling.",
    "CACNA1S": "Skeletal-muscle calcium channel; variants are studied in relation to inhaled anaesthetics.",
    "CFTR":    "Chloride channel; specific variants determine whether CF-targeting drugs can work.",
    "CYP2B6":  "Liver enzyme; main route for efavirenz and methadone.",
    "CYP2C19": "Liver enzyme; activates clopidogrel and clears several antidepressants and PPIs.",
    "CYP2C9":  "Liver enzyme; clears warfarin, phenytoin and many NSAIDs.",
    "CYP2D6":  "Liver enzyme; handles ~20-25% of common drugs (codeine, tamoxifen, many antidepressants). Needs copy-number analysis that short reads cannot provide.",
    "CYP3A4":  "The most abundant liver drug enzyme; contributes to metabolising over half of all drugs.",
    "CYP3A5":  "Liver/intestinal enzyme; most relevant for tacrolimus dosing. Most Europeans carry two non-functional *3 copies.",
    "CYP4F2":  "Vitamin-K pathway enzyme; modest influence on warfarin dose.",
    "DPYD":    "Breaks down fluoropyrimidine chemotherapies (5-FU, capecitabine); reduced activity is a well-known toxicity risk factor.",
    "G6PD":    "Red-blood-cell enzyme protecting against oxidative stress; deficiency is X-linked.",
    "IFNL3":   "Interferon lambda; historically used to predict hepatitis-C treatment response.",
    "MT-RNR1": "Mitochondrial rRNA; certain variants relate to aminoglycoside-associated hearing loss.",
    "NAT2":    "Acetylation enzyme; classic 'slow/fast acetylator' status (isoniazid, hydralazine).",
    "NUDT15":  "Clears thiopurine drugs (azathioprine, mercaptopurine).",
    "RYR1":    "Skeletal-muscle calcium-release channel; studied with anaesthetic response.",
    "SLCO1B1": "Liver uptake transporter for statins; affects how much statin reaches the liver vs the blood.",
    "TPMT":    "Second thiopurine-clearing enzyme, alongside NUDT15.",
    "UGT1A1":  "Glucuronidation enzyme; handles irinotecan and bilirubin (Gilbert's syndrome).",
    "VKORC1":  "Vitamin-K epoxide reductase — warfarin's direct target; the main genetic driver of dose.",
}

PHENOTYPE_MEANING = {
    "ultrarapid metabolizer":   "Much faster than typical — drugs cleared (or activated) more quickly than average.",
    "rapid metabolizer":        "Faster than typical.",
    "normal metabolizer":       "Typical enzyme activity — the population reference point.",
    "intermediate metabolizer": "Somewhat reduced activity — between normal and poor.",
    "poor metabolizer":         "Little to no enzyme activity.",
    "likely intermediate metabolizer": "Probably somewhat reduced activity (call is not definitive).",
    "likely poor metabolizer":  "Probably little to no activity (call is not definitive).",
    "normal function":          "Typical protein function.",
    "decreased function":       "Reduced protein function.",
    "poor function":            "Little to no protein function.",
    "increased function":       "Greater than typical function.",
    "normal":                   "Typical activity.",
    "deficient":                "Markedly reduced activity.",
    "indeterminate":            "The genotype does not map to a defined category.",
    "no result":                "Not callable from this data.",
    "uncertain susceptibility": "No established risk variant found; susceptibility not determined.",
}


def explain_phenotype(p):
    t = str(p).strip().lower()
    if t in (".", "", "nan"):
        return "."
    for key, meaning in PHENOTYPE_MEANING.items():
        if key in t:
            return meaning
    return "."


def gene_function(g):
    return GENE_FUNCTION.get(str(g).strip().upper(), ".")


def diplotype_label(dp):
    """PharmCAT diplotypes carry allele1/allele2 objects rather than a flat label."""
    for key in ("label", "name"):
        if dp.get(key):
            return str(dp[key])
    a1 = (dp.get("allele1") or {}).get("name")
    a2 = (dp.get("allele2") or {}).get("name")
    if a1 and a2:
        return f"{a1}/{a2}"
    return a1 or a2 or "."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows, note = [], None
    try:
        with open(a.json) as fh:
            data = json.load(fh)
        # reporter JSON -> "genes"; phenotyper JSON -> "geneReports"
        genes = data.get("genes") or data.get("geneReports") or {}
        if not isinstance(genes, dict) or not genes:
            note = "PharmCAT produced no gene reports"
        for source, gmap in (genes.items() if isinstance(genes, dict) else []):
            for gene, info in (gmap or {}).items():
                if not isinstance(info, dict):
                    continue
                diplotypes = info.get("recommendationDiplotypes") or info.get("sourceDiplotypes") or []
                for dp in diplotypes:
                    if not isinstance(dp, dict):
                        continue
                    sym = gene or dp.get("gene", ".")
                    pheno = "; ".join(dp.get("phenotypes") or []) or "."
                    act = dp.get("activityScore", ".")
                    rows.append({
                        "gene": sym,
                        "what_this_gene_does": gene_function(sym),
                        "diplotype": diplotype_label(dp),
                        "phenotype": pheno,
                        "what_the_phenotype_means": explain_phenotype(pheno),
                        "activity": "." if act in (None, "", "nan") else act,
                        "source": source,
                        "reference": f"https://www.pharmgkb.org/gene/{sym}" if sym != "." else ".",
                    })
    except Exception as e:
        note = f"could not parse PharmCAT JSON: {e}"

    df = pd.DataFrame(rows, columns=COLUMNS).drop_duplicates()
    if df.empty:
        # Always emit the header row: an empty frame writes a headerless file, and the final
        # report then fails with pandas' "No columns to parse from file".
        df = pd.DataFrame([{c: "." for c in COLUMNS}], columns=COLUMNS)
        df.loc[0, "phenotype"] = note or "no results"
    with open(a.out, "w") as fh:
        fh.write(f"# pharmacogenomics (pgx) for {a.sample} — informational, not a prescription\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[pharmcat_report] {len(df)} gene-phenotype rows"
          + (f" ({note})" if note else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
