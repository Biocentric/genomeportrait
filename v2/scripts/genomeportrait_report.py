#!/usr/bin/env python3
"""
genomeportrait_report.py — assemble all per-sample result tables into a single
self-contained interactive HTML "GenomePortrait" report.

Each upstream module emits a small TSV (or a directory bundle) into parts/. This
script discovers them by filename suffix, renders a themed section per analysis,
embeds any PNG plots as base64, and writes one standalone HTML file.
"""
import argparse
import base64
import datetime as dt
import glob
import html
import os
import sys

import pandas as pd

SECTIONS = [
    # (filename marker,    section id,     title)            order = order in the report
    ("genome_facts",       "facts",        "Genome at a glance"),
    ("ancestry",           "ancestry",     "Ancestry & admixture"),
    ("lineage",            "lineage",      "Maternal & paternal lineage"),
    ("pgs_fragment",       "pgs",          "Polygenic scores"),
    ("traits",             "traits",       "Popular traits"),
    ("hla",                "hla",          "HLA immune type"),
    ("pgx",                "pgx",          "Pharmacogenomics"),
    ("str_profile",        "str",          "Forensic STR profile (CODIS)"),
    ("bcftools_stats",     "variants",     "Variant overview"),
    ("svcnv_summary",      "svcnv",        "Structural & copy-number variants"),
    ("annotation",         "annotation",   "Notable annotated variants"),
    ("extras",             "extras",       "Telomeres & mitochondria"),
]

# Colour tokens are shared with scripts/pgs_overview.py: the polygenic section is an HTML
# fragment that relies on them being defined by the page it is embedded in.
CSS = """
:root { --bg:#fbfaf7; --fg:#1d1d1b; --mut:#6b6a66; --card:#ffffff; --line:#e6e3dc; --axis:#9a978f;
  --good:#1f8a70; --bad:#d1495b; --neu:#8a8780; --good-t:#1f8a7033; --bad-t:#d1495b2e; --neu-t:#8a878026;
  --dot:#3d3c39; --acc:#3b6ea5; --warn:#a56a00; --warn-t:#a56a001f; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#151513; --fg:#ecebe6; --mut:#a3a19a;
  --card:#1f1f1c; --line:#33322e; --axis:#76746d; --good:#4cc3a2; --bad:#ef7a87; --neu:#a3a19a;
  --good-t:#4cc3a233; --bad-t:#ef7a872e; --neu-t:#a3a19a22; --dot:#d8d6cf; --acc:#8ab4e8;
  --warn:#e0b25c; --warn-t:#e0b25c26; } }
* { box-sizing:border-box }
html { scroll-behavior:smooth; scroll-padding-top:56px }
body { margin:0; background:var(--bg); color:var(--fg); font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }
a { color:var(--acc) }
.top { max-width:880px; margin:0 auto; padding:28px 16px 8px }
h1 { font-size:26px; margin:0 0 4px; letter-spacing:-.01em }
.sub { color:var(--mut); margin:0 }
nav { position:sticky; top:0; z-index:5; background:var(--bg); border-bottom:1px solid var(--line) }
nav div { max-width:880px; margin:0 auto; padding:8px 16px; display:flex; gap:6px; overflow-x:auto; white-space:nowrap }
nav a { color:var(--mut); text-decoration:none; font-size:13px; padding:4px 10px; border-radius:99px; border:1px solid var(--line) }
nav a:hover { color:var(--fg); border-color:var(--axis) }
main { max-width:880px; margin:0 auto; padding:8px 16px 64px }
main > section { margin:36px 0 }
main > section > h2 { font-size:19px; margin:0 0 12px; padding-bottom:6px; border-bottom:1px solid var(--line) }
.card { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px 16px; margin:10px 0 }
.tw { overflow-x:auto }
table { border-collapse:collapse; width:100%; font-size:12.5px }
th, td { text-align:left; padding:5px 8px; border-bottom:1px solid var(--line); vertical-align:top }
th { color:var(--mut); font-weight:600; white-space:nowrap }
.badge { display:inline-block; padding:1px 8px; border-radius:99px; font-size:12px; border:1px solid currentColor; white-space:nowrap }
.b-good { color:var(--good) } .b-acc { color:var(--acc) } .b-warn { color:var(--warn) } .b-mut { color:var(--mut) }
.plot { margin:12px 0; text-align:center } .plot img { max-width:100%; border-radius:8px; border:1px solid var(--line); background:#fff }
.muted { color:var(--mut); font-size:13.5px }
.disclaim { font-size:13px; color:var(--mut); border-left:3px solid var(--warn); background:var(--warn-t); padding:8px 12px; margin:10px 0; border-radius:0 8px 8px 0 }
.facts { display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:10px }
.fact { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:12px 14px }
.fact .v { font-size:22px; font-weight:650; letter-spacing:-.01em } .fact .l { font-size:12.5px; color:var(--mut) }
.fact .n { font-size:12px; color:var(--mut); margin-top:4px }
footer { color:var(--mut); font-size:13px; max-width:880px; margin:0 auto; padding:0 16px 48px }
"""


# ---------------------------------------------------------------- annotation ranking
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
# Rough ClinVar ordering; anything unmatched sorts after these.
CLINSIG_RANK = [
    ("pathogenic/likely_pathogenic", 0), ("likely_pathogenic", 1), ("pathogenic", 0),
    ("drug_response", 2), ("risk_factor", 3), ("protective", 3),
    ("conflicting", 5), ("uncertain", 6), ("likely_benign", 8), ("benign", 9),
]
ANNOT_TOP_N = 150


def _clinsig_rank(v):
    t = str(v).strip().lower()
    if t in (".", "", "nan"):
        return 7
    for key, rank in CLINSIG_RANK:
        if key in t:
            return rank
    return 7


def call_quality(row):
    """Flag calls whose support looks weak, so ranking can't imply confidence.

    Rare frameshift/splice calls are the most error-prone class in short-read data, and a
    HIGH-impact label says nothing about whether the call is real. Depth, genotype quality
    and the variant allele fraction are what distinguish a genuine call from an artifact:
    a heterozygote should sit near 0.5 and a homozygote near 1.0."""
    notes = []
    try:
        dp = pd.to_numeric(row.get("DEPTH"), errors="coerce")
        gq = pd.to_numeric(row.get("GQ"), errors="coerce")
        vaf = pd.to_numeric(str(row.get("VAF")).split(",")[0], errors="coerce")
        gt = str(row.get("GT", ""))
        if pd.notna(dp):
            if dp < 10:
                notes.append(f"low depth ({int(dp)}x)")
            elif dp > 100:
                notes.append(f"very high depth ({int(dp)}x — repeat/CNV region?)")
        if pd.notna(gq) and gq < 20:
            notes.append(f"low genotype quality (GQ {int(gq)})")
        if pd.notna(vaf):
            het = gt in ("0/1", "1/0", "0|1", "1|0")
            hom = gt in ("1/1", "1|1")
            if het and not (0.30 <= vaf <= 0.70):
                notes.append(f"VAF {vaf:.2f} unusual for a heterozygote")
            elif hom and vaf < 0.85:
                notes.append(f"VAF {vaf:.2f} low for a homozygote")
    except Exception:
        pass
    return "; ".join(notes) if notes else "looks well supported"


def rank_annotation(df):
    """Pick the most informative variants instead of the first N rows on chr1.

    The raw table is one row per VEP consequence in genomic order, so head() only ever
    shows telomeric MODIFIER noise. Collapse to one row per variant (keeping its most
    severe consequence) and rank by predicted impact, then ClinVar significance, then
    rarity."""
    if df is None or len(df) == 0:
        return df, 0
    d = df.copy()
    d.columns = [str(c).strip() for c in d.columns]
    if "IMPACT" not in d.columns:
        return d, len(d)
    d["_imp"] = d["IMPACT"].map(lambda v: IMPACT_RANK.get(str(v).strip().upper(), 4))
    d["_cln"] = d["CLINSIG"].map(_clinsig_rank) if "CLINSIG" in d.columns else 7
    d["_af"] = pd.to_numeric(d.get("GNOMAD_AF"), errors="coerce").fillna(1.0)         if "GNOMAD_AF" in d.columns else 1.0
    d = d.sort_values(["_imp", "_cln", "_af"], kind="mergesort")
    keys = [c for c in ("CHROM", "POS", "REF", "ALT") if c in d.columns]
    if keys:
        d = d.drop_duplicates(subset=keys, keep="first")   # keep worst consequence per variant
    total = len(d)
    d = d.head(ANNOT_TOP_N).drop(columns=["_imp", "_cln", "_af"], errors="ignore")
    if "DEPTH" in d.columns or "GQ" in d.columns or "VAF" in d.columns:
        d["call_quality"] = d.apply(call_quality, axis=1)
    return d, total


def linkify(col, value):
    """Turn identifiers into outbound links so they can be looked up."""
    v = html.escape(str(value))
    if v in (".", "", "nan"):
        return v
    if col == "DBSNP":
        ids = [x for x in str(value).replace(",", "&").split("&") if x.startswith("rs")]
        if ids:
            return ", ".join(
                f'<a href="https://www.ncbi.nlm.nih.gov/snp/{html.escape(i)}" target="_blank" rel="noopener">{html.escape(i)}</a>'
                for i in ids)
        return v
    if col == "GENE":
        return f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene={v}" target="_blank" rel="noopener">{v}</a>'
    if col == "CLINSIG":
        return f'<span class="badge b-warn">{v}</span>' if _clinsig_rank(value) <= 3 else v
    if col == "IMPACT":
        cls = {"HIGH": "b-warn", "MODERATE": "b-acc"}.get(v.upper(), "b-mut")
        return f'<span class="badge {cls}">{v}</span>'
    return v


def b64_img(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


LINK_COLS = {"DBSNP", "GENE", "CLINSIG", "IMPACT"}


def df_to_html(df, max_rows=200, rich=False):
    if df is None or len(df) == 0:
        return '<p class="muted">No data produced for this section.</p>'
    df = df.head(max_rows)
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = "".join(
            f"<td>{linkify(c, v) if rich and c in LINK_COLS else html.escape(str(v))}</td>"
            for c, v in zip(df.columns, r))
        rows.append(f"<tr>{cells}</tr>")
    return f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def facts_html(df):
    """Genome facts as a grid of figures rather than a table."""
    cells = []
    for _, r in df.iterrows():
        cells.append(f'<div class="fact"><div class="v">{html.escape(str(r["value"]))}</div>'
                     f'<div class="l">{html.escape(str(r["fact"]))}</div>'
                     f'<div class="n">{html.escape(str(r["note"]))}</div></div>')
    return f'<div class="facts">{"".join(cells)}</div>'


def bcftools_stats_table(path):
    """Summary numbers from `bcftools stats` (the SN and TSTV lines)."""
    rows = []
    try:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "SN" and len(f) >= 4:
                    rows.append((f[2].rstrip(":").replace("number of ", "").capitalize(), f"{int(f[3]):,}"))
                elif f[0] == "TSTV" and len(f) >= 5:
                    rows.append(("Transition / transversion ratio", f[4]))
    except Exception as e:
        rows.append(("error", str(e)))
    return pd.DataFrame(rows, columns=["metric", "value"])


def read_part(path):
    """Read a TSV part, or collect a bundle directory (tables + plots)."""
    tables, plots = [], []
    if os.path.isdir(path):
        for tsv in sorted(glob.glob(os.path.join(path, "*.tsv"))):
            if tsv.endswith(".raw.tsv"):
                continue   # bulk data kept on disk, too large/raw to render
            try:
                # comment="#": bundle TSVs carry a leading "# ..." title line. Without this it
                # is parsed as the header row, pushing the real header down into the data.
                tables.append((os.path.basename(tsv), pd.read_csv(tsv, sep="\t", comment="#")))
            except Exception:
                pass
        for png in sorted(glob.glob(os.path.join(path, "*.png"))):
            plots.append(b64_img(png))
    else:
        try:
            tables.append((os.path.basename(path), pd.read_csv(path, sep="\t", comment="#")))
        except Exception as e:
            tables.append((os.path.basename(path), pd.DataFrame({"error": [str(e)]})))
    return tables, plots


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--parts-dir", required=True)
    ap.add_argument("--multiqc", default=None)
    ap.add_argument("--pipeline-version", default="dev")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    part_files = sorted(glob.glob(os.path.join(a.parts_dir, "*")))

    # Group discovered parts by which section their filename matches
    by_section = {sid: [] for _, sid, _ in SECTIONS}
    for pf in part_files:
        if ".raw." in os.path.basename(pf):
            continue   # bulk artifact (full annotation table etc.) — published, not rendered
        base = os.path.basename(pf)
        for suffix, sid, _ in SECTIONS:
            if suffix in base:
                by_section[sid].append(pf)
                break

    nav_links, body = [], []
    for suffix, sid, title in SECTIONS:
        files = by_section.get(sid, [])
        if not files:
            continue
        nav_links.append(f'<a href="#{sid}">{html.escape(title)}</a>')
        chunks = [f'<section id="{sid}"><h2>{html.escape(title)}</h2>']
        for pf in files:
            if sid == "pgs":
                # a ready-made HTML fragment (scripts/pgs_overview.py), embedded as it is
                with open(pf, encoding="utf-8") as fh:
                    chunks.append(fh.read())
                continue
            if sid == "variants":
                chunks.append(df_to_html(bcftools_stats_table(pf)))
                continue
            tables, plots = read_part(pf)
            for name, df in tables:
                if sid == "facts" and {"fact", "value", "note"} <= set(df.columns):
                    chunks.append(facts_html(df))
                elif sid == "annotation":
                    df, total = rank_annotation(df)
                    if total > len(df):
                        chunks.append(
                            f'<p class="muted">Showing the {len(df)} most impactful of '
                            f'{total:,} annotated variants — one row per variant (most severe '
                            f'consequence), ranked by predicted impact, then ClinVar '
                            f'significance, then rarity. rsIDs and genes link out. The full '
                            f'table is in the results folder.</p>'
                            '<div class="disclaim"><b>How to read this — it is not a health '
                            'finding.</b> "HIGH impact" is a prediction about the effect on a '
                            'protein sequence, not evidence that anything is wrong. A ClinVar '
                            '"Pathogenic" label describes a variant as recorded in a public '
                            'database, usually established in patients with a specific '
                            'condition; for recessive genes it most often indicates ordinary '
                            'carrier status, which is unremarkable — everyone carries some. '
                            'These calls are unconfirmed short-read genotypes: rare '
                            'frameshift/splice calls are the most error-prone class, and none '
                            'have been validated. Do not draw conclusions from this table.</div>')
                    chunks.append(df_to_html(df, max_rows=ANNOT_TOP_N, rich=True))
                else:
                    chunks.append(df_to_html(df))
            for img in plots:
                chunks.append(f'<div class="plot"><img src="{img}" alt=""></div>')
        if sid in ("pgx", "traits", "annotation"):
            chunks.append('<div class="disclaim">Recreational / educational interpretation only. '
                          'Not a medical or diagnostic result.</div>')
        chunks.append("</section>")
        body.append("".join(chunks))

    mqc_link = ""
    if a.multiqc:
        mqc_link = (f'<section id="qc"><h2>Technical QC</h2>'
                    f'<p>Read, alignment and coverage QC is in the bundled '
                    f'<a href="{html.escape(os.path.basename(a.multiqc))}">MultiQC report</a>.</p></section>')
        nav_links.append('<a href="#qc">Technical QC</a>')

    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GenomePortrait — {html.escape(a.sample)}</title><style>{CSS}</style></head><body>
<div class="top"><h1>GenomePortrait</h1>
<p class="sub">{html.escape(a.sample)} &middot; generated {now} &middot; genomeportrait v{html.escape(a.pipeline_version)}</p></div>
<nav><div>{''.join(nav_links)}</div></nav>
<main>
<div class="card"><p style="margin:0 0 6px">A recreational, ancestry and educational portrait of one genome. Each section is computed by an open-source tool in a reproducible container; the numbers behind every section are in the results folder.</p>
<div class="disclaim" style="margin-bottom:0"><b>Not medical advice.</b> Nothing here is a diagnosis. Polygenic scores are inherited tendencies, not measurements; pharmacogenomic and annotation sections are context, not instructions. Do not change treatment or make medical decisions from this report.</div></div>
{''.join(body)}
{mqc_link}
</main>
<footer>genomeportrait v{html.escape(a.pipeline_version)} &middot; Nextflow &middot; for recreational use</footer>
</body></html>"""

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"[genomeportrait] wrote {a.out} with {len(body)} sections", file=sys.stderr)


if __name__ == "__main__":
    main()
