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

# CODIS loci whose hg38 reference tract is longer than a 2x150 bp read pair can span,
# so HipSTR cannot genotype them and omits them without a trace. ExpansionHunter covers
# these instead, from the merged catalog (assets/codis_long_loci_hg38.json).
# Canonical forensic casing (vWA, not VWA); matched case-insensitively.
CODIS_LONG = ["D21S11", "D2S1338", "FGA", "vWA", "D12S391"]


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


def _locus_of(repid):
    """Compound loci emit one record per segment (vWA_GATA_1, vWA_GACA_2) - map back."""
    up = str(repid).upper()
    for name in CODIS_LONG:
        n = name.upper()
        if up == n or up.startswith(n + "_"):
            return name
    return None


def parse_eh(path):
    """Return (codis_long_rows, expansion_rows).

    A compound locus is genotyped as several adjacent repeat segments, so sum REPCN across a
    locus's segments to get its total repeat count. EH genotypes each segment independently
    (there is no phasing between them), so the per-allele sum is an approximation and the
    per-segment detail is reported alongside it rather than in place of it.
    """
    expansions, comp = [], {}
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
                repid = info.get("REPID") or info.get("VARID") or f[2]
                locus = _locus_of(repid)
                if locus is None:
                    expansions.append({
                        "locus": repid,
                        "chrom": f[0],
                        "repeat_unit": info.get("RU", "."),
                        "repeat_count": str(smp.get("REPCN", ".")).replace("/", " / "),
                        "reads": smp.get("LC", smp.get("DP", ".")),
                    })
                    continue
                d = comp.setdefault(locus, {
                    "chrom": f[0], "pos": f[1], "period": len(info.get("RU", "")) or 4,
                    "alleles": [], "segs": [], "reads": smp.get("LC", ".")})
                counts = [c for c in str(smp.get("REPCN", "")).split("/") if c not in ("", ".")]
                d["alleles"].append([int(c) for c in counts] if all(c.isdigit() for c in counts) else [])
                d["segs"].append("%sx%s" % (info.get("RU", "."), smp.get("REPCN", ".")))
    except Exception as e:
        print("[str_report] expansionhunter parse error: %s" % e, file=sys.stderr)

    codis_long = []
    for locus, d in sorted(comp.items()):
        per = [a for a in d["alleles"] if a]
        n = max((len(a) for a in per), default=0)
        totals = [sum(a[i] if i < len(a) else a[-1] for a in per) for i in range(n)] if per else []
        codis_long.append({
            "marker": locus,
            "chrom": d["chrom"],
            "pos": d["pos"],
            "repeat_unit_bp": d["period"],
            "allele_1": totals[0] if totals else ".",
            "allele_2": totals[1] if len(totals) > 1 else (totals[0] if totals else "."),
            "reads": d["reads"],
            "source": "ExpansionHunter",
            "segments": " + ".join(d["segs"]),
        })
    return codis_long, expansions


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
    for r in hip:
        r.setdefault("source", "HipSTR")
    codis = [r for r in hip if str(r["marker"]).upper() in CODIS]
    ystr = [r for r in hip if str(r["marker"]).upper().startswith(("DYS", "DYF", "DYZ"))]
    other = [r for r in hip if r not in codis and r not in ystr]
    codis_long, eh = parse_eh(a.expansionhunter)

    # The long loci belong in the CODIS profile, not buried among disease repeat expansions.
    # HipSTR cannot reach them, so they can only come from ExpansionHunter; if HipSTR somehow
    # did emit one, keep its call and don't duplicate the marker.
    have = {str(r["marker"]).upper() for r in codis}
    added = [r for r in codis_long if str(r["marker"]).upper() not in have]
    codis = sorted(codis + added, key=lambda r: str(r["marker"]))

    missing = sorted(CODIS - {str(r["marker"]).upper() for r in codis})
    write(codis, os.path.join(a.outdir, "1_codis_profile.tsv"),
          f"CODIS core STR profile for {a.sample} - alleles in repeat units "
          f"(x.y = y extra bases beyond x full repeats). "
          f"{len(codis)}/{len(CODIS)} core loci recovered"
          + (f"; still missing: {', '.join(missing)}" if missing else "")
          + ". Loci longer than a read pair can span come from ExpansionHunter (see 'source');"
          " for those, compound repeats are summed across segments, so the total is approximate.")
    if ystr:
        write(ystr, os.path.join(a.outdir, "2_y_str_haplotype.tsv"),
              f"Y-STR haplotype for {a.sample} - paternal-line markers")
    write(eh, os.path.join(a.outdir, "3_repeat_expansions.tsv"),
          f"Repeat-expansion loci for {a.sample} (ExpansionHunter)")
    # The remaining panel markers are bulk data: kept on disk, not rendered.
    write(other, os.path.join(a.outdir, "4_other_str_markers.raw.tsv"),
          f"Other genotyped STR markers for {a.sample}")

    print("[str_report] CODIS %d/%d (%d HipSTR + %d ExpansionHunter), %d Y-STR, %d other, "
          "%d expansion loci%s"
          % (len(codis), len(CODIS), len(codis) - len(added), len(added), len(ystr), len(other),
             len(eh), ("; missing: " + ",".join(missing)) if missing else ""),
          file=sys.stderr)


if __name__ == "__main__":
    main()