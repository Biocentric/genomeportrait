#!/usr/bin/env python3
"""str_report.py — CODIS-style STR profile, Y-STR haplotype, and repeat-expansion sizes.

Each goes into its OWN table in a bundle directory. Writing them into one TSV made pandas
stack the expansion rows underneath the CODIS columns, producing a table that "bled" into
the next with nan/"." cells.

HipSTR reports genotypes as GB — the base-pair difference of each allele from the reference
locus. Forensic STR profiles are quoted in REPEAT UNITS (and use the "9.3" convention for a
partial repeat), so convert: allele_bp = len(REF) + GB, then split by the repeat period.
"""
import argparse
import gzip
import os
import sys

import pandas as pd

# Expanded CODIS core loci (the markers used in forensic identity databases)
CODIS = {"CSF1PO", "FGA", "TH01", "TPOX", "VWA", "D3S1358", "D5S818", "D7S820",
         "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", "D1S1656",
         "D2S441", "D2S1338", "D10S1248", "D12S391", "D19S433", "D22S1045"}


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def allele_repeats(ref_len, gb, period):
    """Forensic nomenclature: full repeats, plus '.N' for a partial repeat of N bases."""
    try:
        total = int(ref_len) + int(gb)
        period = int(period) or 1
        if total < 0:
            return "."
        full, rem = divmod(total, period)
        return f"{full}.{rem}" if rem else str(full)
    except (ValueError, TypeError):
        return "."


def parse_hipstr(path):
    rows = []
    try:
        with opener(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                smp = dict(zip(f[8].split(":"), f[9].split(":")))
                period = info.get("PERIOD", "4")
                ref_len = len(f[3])
                gb = smp.get("GB", "")
                alleles = []
                if gb and gb not in (".", "./."):
                    for part in gb.replace("/", "|").split("|"):
                        alleles.append(allele_repeats(ref_len, part, period))
                rows.append({
                    "marker": info.get("GENE") or f[2],
                    "chrom": f[0],
                    "pos": f[1],
                    "repeat_unit_bp": period,
                    "allele_1": alleles[0] if alleles else ".",
                    "allele_2": alleles[1] if len(alleles) > 1 else (alleles[0] if alleles else "."),
                    "reads": smp.get("DP", "."),
                })
    except Exception as e:
        print(f"[str_report] hipstr parse error: {e}", file=sys.stderr)
    return rows


def parse_eh(path):
    rows = []
    try:
        with opener(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                smp = dict(zip(f[8].split(":"), f[9].split(":")))
                rows.append({
                    "locus": info.get("REPID") or info.get("GENE") or f[2],
                    "chrom": f[0],
                    "repeat_unit": info.get("RU", "."),
                    "repeat_count": smp.get("REPCN", ".").replace("/", " / "),
                    "reads": smp.get("LC", smp.get("DP", ".")),
                })
    except Exception as e:
        print(f"[str_report] expansionhunter parse error: {e}", file=sys.stderr)
    return rows


def write(df_rows, path, header):
    with open(path, "w") as fh:
        fh.write(f"# {header}\n")
    (pd.DataFrame(df_rows) if df_rows else
     pd.DataFrame([{"note": "no loci genotyped"}])).to_csv(path, sep="\t", index=False, mode="a")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hipstr", required=True)
    ap.add_argument("--expansionhunter", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    hip = parse_hipstr(a.hipstr)
    codis = [r for r in hip if str(r["marker"]).upper() in CODIS]
    ystr = [r for r in hip if str(r["marker"]).upper().startswith(("DYS", "DYF", "DYZ"))]
    other = [r for r in hip if r not in codis and r not in ystr]
    eh = parse_eh(a.expansionhunter)

    write(codis, os.path.join(a.outdir, "1_codis_profile.tsv"),
          f"CODIS core STR profile for {a.sample} — alleles in repeat units "
          f"(x.y = y extra bases beyond x full repeats)")
    if ystr:
        write(ystr, os.path.join(a.outdir, "2_y_str_haplotype.tsv"),
              f"Y-STR haplotype for {a.sample} — paternal-line markers")
    write(eh, os.path.join(a.outdir, "3_repeat_expansions.tsv"),
          f"Repeat-expansion loci for {a.sample} (ExpansionHunter)")
    # The remaining panel markers are bulk data: kept on disk, not rendered.
    write(other, os.path.join(a.outdir, "4_other_str_markers.raw.tsv"),
          f"Other genotyped STR markers for {a.sample}")

    print(f"[str_report] {len(codis)} CODIS, {len(ystr)} Y-STR, {len(other)} other, "
          f"{len(eh)} expansion loci", file=sys.stderr)


if __name__ == "__main__":
    main()
