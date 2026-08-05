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
    # (suffix,              section id,     title,                              emoji)
    ("variant_stats",      "variants",     "Genome & variant overview",        "🧬"),
    ("annotation",         "annotation",   "Notable annotated variants",       "📝"),
    ("svcnv_summary",      "svcnv",        "Structural & copy-number variants","🧩"),
    ("ancestry",           "ancestry",     "Ancestry & admixture",             "🌍"),
    ("lineage",            "lineage",      "Maternal & paternal lineage",      "🧭"),
    ("str_profile",        "str",          "Forensic STR profile (CODIS)",     "🔬"),
    ("prs",                "prs",          "Polygenic trait scores",           "📊"),
    ("hla",                "hla",          "HLA immune type",                  "🛡️"),
    ("pgx",                "pgx",          "Pharmacogenomics",                 "💊"),
    ("traits",             "traits",       "Popular traits",                   "✨"),
    ("extras",             "extras",       "Telomeres & mitochondria",         "⏳"),
]

CSS = """
:root{--bg:#0f1220;--card:#1a1f33;--ink:#e8ecf5;--muted:#9aa3bd;--accent:#5cc8ff;--accent2:#b388ff;--good:#5ad19a;--warn:#ffd166}
*{box-sizing:border-box}body{margin:0;background:linear-gradient(160deg,#0f1220,#161a2e);color:var(--ink);font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}
header{padding:42px 28px 26px;background:radial-gradient(1200px 300px at 20% -50%,rgba(92,200,255,.25),transparent),radial-gradient(900px 280px at 90% -40%,rgba(179,136,255,.22),transparent)}
h1{margin:0;font-size:30px;letter-spacing:.3px}h1 .sub{display:block;font-size:14px;color:var(--muted);font-weight:400;margin-top:6px}
.wrap{max-width:1080px;margin:0 auto;padding:0 22px 80px}
nav{position:sticky;top:0;z-index:5;background:rgba(15,18,32,.85);backdrop-filter:blur(8px);display:flex;flex-wrap:wrap;gap:8px;padding:12px 22px;border-bottom:1px solid #222a44}
nav a{color:var(--muted);text-decoration:none;font-size:13px;padding:5px 10px;border-radius:999px;border:1px solid transparent}
nav a:hover{color:var(--ink);border-color:#2c3556;background:#161b30}
section{background:var(--card);border:1px solid #232a45;border-radius:16px;padding:22px 24px;margin:22px 0;box-shadow:0 10px 30px rgba(0,0,0,.25)}
section h2{margin:.1em 0 .6em;font-size:20px}
table{width:100%;border-collapse:collapse;font-size:13.5px;margin:10px 0;overflow:hidden;border-radius:10px}
th,td{padding:8px 10px;text-align:left;border-bottom:1px solid #262e4d}
th{background:#202747;color:#cdd6f4;position:sticky;top:0}
tr:hover td{background:#1d2440}
.badge{display:inline-block;padding:2px 9px;border-radius:999px;font-size:12px;font-weight:600}
.b-good{background:rgba(90,209,154,.15);color:var(--good)}.b-acc{background:rgba(92,200,255,.15);color:var(--accent)}
.b-warn{background:rgba(255,209,102,.15);color:var(--warn)}.b-mut{background:#222a45;color:var(--muted)}
.plot{margin:14px 0;text-align:center}.plot img{max-width:100%;border-radius:12px;border:1px solid #2a3357}
.muted{color:var(--muted)}.disclaim{font-size:12.5px;color:var(--muted);border-left:3px solid var(--accent2);padding:6px 12px;margin-top:8px;background:#171c30;border-radius:0 8px 8px 0}
footer{color:var(--muted);font-size:12px;text-align:center;padding:30px}
.kpis{display:flex;flex-wrap:wrap;gap:14px;margin:6px 0 2px}
.kpi{flex:1 1 150px;background:#192041;border:1px solid #283254;border-radius:12px;padding:14px 16px}
.kpi .v{font-size:24px;font-weight:700}.kpi .l{font-size:12px;color:var(--muted)}
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
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def read_part(path):
    """Read a TSV part, or collect a bundle directory (tables + plots)."""
    tables, plots = [], []
    if os.path.isdir(path):
        for tsv in sorted(glob.glob(os.path.join(path, "*.tsv"))):
            if tsv.endswith(".raw.tsv"):
                continue   # bulk data kept on disk, too large/raw to render
            try:
                tables.append((os.path.basename(tsv), pd.read_csv(tsv, sep="\t")))
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
    by_section = {sid: [] for _, sid, _, _ in SECTIONS}
    for pf in part_files:
        if ".raw." in os.path.basename(pf):
            continue   # bulk artifact (full annotation table etc.) — published, not rendered
        base = os.path.basename(pf)
        for suffix, sid, _, _ in SECTIONS:
            if suffix in base:
                by_section[sid].append(pf)
                break

    nav_links, body = [], []
    for suffix, sid, title, emoji in SECTIONS:
        files = by_section.get(sid, [])
        if not files:
            continue
        nav_links.append(f'<a href="#{sid}">{emoji} {html.escape(title)}</a>')
        chunks = [f'<section id="{sid}"><h2>{emoji} {html.escape(title)}</h2>']
        for pf in files:
            tables, plots = read_part(pf)
            for name, df in tables:
                if sid == "annotation":
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
                chunks.append(f'<div class="plot"><img src="{img}"></div>')
        if sid in ("prs", "pgx", "traits", "annotation"):
            chunks.append('<div class="disclaim">Recreational / educational interpretation only. '
                          'Not a medical or diagnostic result.</div>')
        chunks.append("</section>")
        body.append("".join(chunks))

    mqc_link = ""
    if a.multiqc:
        mqc_link = (f'<section id="qc"><h2>📈 Technical QC</h2>'
                    f'<p>Full read/alignment QC is in the bundled '
                    f'<a class="badge b-acc" href="{html.escape(os.path.basename(a.multiqc))}">MultiQC report ↗</a></p></section>')
        nav_links.append('<a href="#qc">📈 Technical QC</a>')

    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GenomePortrait — {html.escape(a.sample)}</title><style>{CSS}</style></head><body>
<header><div class="wrap"><h1>GenomePortrait
<span class="sub">Sample <b>{html.escape(a.sample)}</b> · generated {now} · pipeline v{html.escape(a.pipeline_version)}</span></h1></div></header>
<nav>{''.join(nav_links)}</nav>
<div class="wrap">
<section><h2>About this report</h2>
<p class="muted">A recreational, ancestry and educational portrait of a personal genome. Every section below is computed by an open-source, peer-reviewed tool running in a reproducible Singularity container. This report deliberately avoids clinical disease-risk predictions.</p>
<div class="disclaim"><b>Not medical advice.</b> Do not use any part of this report for diagnosis, treatment, or other medical decisions. Pharmacogenomic and annotation sections are provided for context and curiosity only.</div>
</section>
{''.join(body)}
{mqc_link}
</div>
<footer>Generated by genomeportrait v{html.escape(a.pipeline_version)} · Nextflow + Singularity · for recreational use</footer>
</body></html>"""

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"[genomeportrait] wrote {a.out} with {len(body)} sections", file=sys.stderr)


if __name__ == "__main__":
    main()
