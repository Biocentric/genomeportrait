#!/usr/bin/env python3
"""Render interpreted polygenic scores as a self-contained HTML overview.

    python3 pgs_overview.py S.pgs_interpreted.json --html S.pgs_overview.html --fragment S.pgs_fragment.html

--html      a complete page that opens in any browser
--fragment  the same content as a <div class="pgs"> with scoped CSS, embedded by the
            unified GenomePortrait report (which defines the same colour tokens)

Layout and wording follow prs-selfscore/interpret/build_html.py. The hand-written
"contradictions resolved" box of that page is replaced by notes derived from the data,
so the page is correct for any genome.
"""
import argparse
import html
import json
import math

e = html.escape

TOKENS = """
:root { --bg:#fbfaf7; --fg:#1d1d1b; --mut:#6b6a66; --card:#ffffff; --line:#e6e3dc; --axis:#9a978f;
  --good:#1f8a70; --bad:#d1495b; --neu:#8a8780; --good-t:#1f8a7033; --bad-t:#d1495b2e; --neu-t:#8a878026;
  --dot:#3d3c39; --acc:#3b6ea5; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#151513; --fg:#ecebe6; --mut:#a3a19a;
  --card:#1f1f1c; --line:#33322e; --axis:#76746d; --good:#4cc3a2; --bad:#ef7a87; --neu:#a3a19a;
  --good-t:#4cc3a233; --bad-t:#ef7a872e; --neu-t:#a3a19a22; --dot:#d8d6cf; --acc:#8ab4e8; } }
:root[data-theme="dark"] { --bg:#151513; --fg:#ecebe6; --mut:#a3a19a; --card:#1f1f1c; --line:#33322e;
  --axis:#76746d; --good:#4cc3a2; --bad:#ef7a87; --neu:#a3a19a; --good-t:#4cc3a233; --bad-t:#ef7a872e;
  --neu-t:#a3a19a22; --dot:#d8d6cf; --acc:#8ab4e8; }
"""

PAGE_CSS = """
* { box-sizing:border-box }
body { margin:0; background:var(--bg); color:var(--fg); font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:880px; margin:0 auto; padding:28px 16px 64px; }
h1 { font-size:26px; margin:0 0 4px; letter-spacing:-.01em }
"""

# Every rule is scoped to .pgs so the fragment cannot restyle the report around it.
PGS_CSS = """
.pgs h2 { font-size:18px; margin:40px 0 12px; padding-bottom:6px; border-bottom:1px solid var(--line) }
.pgs h3 { font-size:16px; margin:0 }
.pgs .sub { color:var(--mut); margin:0 0 20px }
.pgs .box { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:16px 18px; margin:16px 0 }
.pgs .box h2 { border:0; margin:0 0 8px; padding:0 }
.pgs .box ul { margin:6px 0 0; padding-left:20px }
.pgs .box li { margin:4px 0 }
.pgs .trait { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px 16px; margin:10px 0 }
.pgs .trait header { display:flex; justify-content:space-between; align-items:center; gap:8px; padding:0; background:none }
.pgs .trait p { margin:6px 0 }
.pgs .bar { width:100%; max-width:520px; height:auto; display:block; margin:8px 0 0 }
.pgs .scale { display:flex; justify-content:space-between; gap:12px; max-width:520px; font-size:12px; color:var(--mut) }
.pgs .chip { font-size:12px; padding:2px 9px; border-radius:99px; white-space:nowrap; border:1px solid currentColor }
.pgs .chip.good { color:var(--good) } .pgs .chip.bad { color:var(--bad) } .pgs .chip.neu { color:var(--neu) }
.pgs .note { color:var(--mut); font-size:13.5px }
.pgs .flags { color:var(--mut); font-size:12.5px; font-style:italic }
.pgs .mono { font-family:ui-monospace,Consolas,monospace; font-size:12.5px }
.pgs details summary { cursor:pointer; color:var(--mut); font-size:13px; margin-top:4px }
.pgs .tw { overflow-x:auto }
.pgs table { border-collapse:collapse; width:100%; font-size:12.5px; margin-top:8px }
.pgs th, .pgs td { text-align:left; padding:4px 8px; border-bottom:1px solid var(--line); vertical-align:top; background:none }
.pgs th { color:var(--mut); font-weight:600 }
.pgs td.num { text-align:right; font-variant-numeric:tabular-nums; white-space:nowrap }
.pgs tr.dim td { color:var(--mut); font-style:italic }
.pgs dl { display:grid; grid-template-columns:max-content 1fr; gap:4px 14px; margin:8px 0 0 }
.pgs dt { font-weight:600 } .pgs dd { margin:0 }
.pgs .limits { color:var(--mut); font-size:13px; margin-top:40px }
@media (max-width:560px) { .pgs dl { grid-template-columns:1fr } .pgs dd { margin-bottom:6px } }
"""


def sizeword(sd):
    a = abs(sd)
    return "negligible" if a < 0.1 else "small" if a < 0.25 else "moderate" if a < 0.5 else "large"


def ordinal(p):
    p = round(p)
    suf = "th" if 10 <= p % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(p % 10, "th")
    return f"{p}{suf}"


def chip_class(v):
    return {"favourable": "good", "unfavourable": "bad"}.get(v, "neu")


def bar(s, xs):
    """Diverging scale -3..+3; the track is tinted by which side is healthier."""
    W, H, pad = 320, 30, 10

    def X(z):
        return pad + (max(-3, min(3, z)) + 3) / 6 * (W - 2 * pad)

    good = s["good"]
    left = "var(--bad-t)" if good > 0 else "var(--good-t)" if good < 0 else "var(--neu-t)"
    right = "var(--good-t)" if good > 0 else "var(--bad-t)" if good < 0 else "var(--neu-t)"
    half = (W - 2 * pad) / 2
    parts = [f'<svg viewBox="0 0 {W} {H}" class="bar" role="img" aria-label="consensus Z {s["zc"]:+.2f}">',
             f'<rect x="{pad}" y="11" width="{half}" height="8" rx="4" fill="{left}"/>',
             f'<rect x="{W / 2}" y="11" width="{half}" height="8" rx="4" fill="{right}"/>',
             f'<line x1="{W / 2}" x2="{W / 2}" y1="6" y2="24" stroke="var(--axis)" stroke-width="1"/>']
    for z in (-2, -1, 1, 2):
        parts.append(f'<line x1="{X(z)}" x2="{X(z)}" y1="9" y2="21" stroke="var(--axis)" stroke-width=".5" opacity=".6"/>')
    for x in xs:
        fill = "var(--dot)" if x["eur"] else "none"
        parts.append(f'<circle cx="{X(x["z"]):.1f}" cy="15" r="3" fill="{fill}" stroke="var(--dot)" '
                     f'stroke-width="1" opacity=".75"><title>{e(x["acc"])} {e(x["pub"])}: Z {x["z"]:+.2f}</title></circle>')
    col = {"favourable": "var(--good)", "unfavourable": "var(--bad)"}.get(s["verdict"], "var(--neu)")
    cx = X(s["zc"])
    parts.append(f'<path d="M{cx:.1f} 4 L{cx + 7:.1f} 15 L{cx:.1f} 26 L{cx - 7:.1f} 15 Z" fill="{col}" '
                 f'stroke="var(--bg)" stroke-width="1.5"/></svg>')
    return "".join(parts)


def effect_text(x):
    if x.get("shift_sd") is not None:
        return x.get("units") or f'{x["shift_sd"]:+.2f} SD'
    if x.get("rel_odds"):
        return f'odds &times;{x["rel_odds"]:.2f}'
    return ""


def trait_card(s, xs):
    h = s["head"]
    used = [x for x in xs if x["eur"]] or xs
    direction = s["higher"] if s["zc"] > 0 else s["lower"]
    if abs(s["zc"]) < 0.5:
        tendency = f"<b>Around average</b> (consensus Z {s['zc']:+.2f})"
    else:
        tendency = f"<b>{e(s['strength'].capitalize())} towards {e(direction)}</b> (consensus Z {s['zc']:+.2f})"
    src = f'<span class="mono">{e(h["acc"])}</span> ({e(h["pub"])}, your {ordinal(h["pct"])} percentile)'
    if s["kind"] == "q" and h.get("shift_sd") is not None:
        amount = f"<b>{e(h['units'])}</b>" if h.get("units") else f"<b>{h['shift_sd']:+.2f} SD</b>"
        effect = (f"Best-validated score {src} predicts {amount} versus an average person "
                  f"&mdash; a {sizeword(h['shift_sd'])} shift. It explains <b>{h['r2'] * 100:.1f}%</b> of the "
                  f"differences between people; the rest is lifestyle, age, environment and genetics it doesn't capture.")
    elif s["kind"] == "b" and h.get("rel_odds"):
        auc = f", AUC {h['auc']:.2f}" if h.get("auc") else ""
        effect = (f"Best-validated score {src}: odds about <b>{h['rel_odds']:.2f}&times;</b> an average person's "
                  f"(published OR {h['or']:.2f} per SD{auc}).")
    else:
        effect = "No published validation that lets this be translated into a size &mdash; direction only."
    flags = []
    n_flip = sum(x["sign"] == -1 for x in xs)
    if n_flip:
        how = ", verified by weight correlation" if all(x["verified"] for x in xs if x["sign"] == -1) else ""
        flags.append(f"{n_flip} score{'s' if n_flip > 1 else ''} coded in the opposite direction &mdash; flipped{how}")
    nn = sum(not x["eur"] for x in xs)
    if nn and nn < len(xs):
        flags.append(f"{nn} non-European-trained score{'s' if nn > 1 else ''} shown hollow, excluded from consensus")
    elif nn:
        flags.append("no European-trained score available; consensus uses all scores")
    if any(x.get("r2_note") for x in xs):
        flags.append("some reported R² not used (in-sample or above heritability)")
    rows = []
    for x in sorted(xs, key=lambda x: (not x["eur"], -abs(x["z"]))):
        r2_txt = "{:.1f}%".format(x["r2"] * 100) if x["r2"] else ""
        rows.append(
            f'<tr class="{"" if x["eur"] else "dim"}"><td class="mono">{e(x["acc"])}</td><td>{e(x["trait"])}</td>'
            f'<td>{e(x["pub"])}</td><td class="num">{x["pct"]:.1f}</td><td class="num">{x["z"]:+.2f}</td>'
            f'<td class="num">{r2_txt}</td>'
            f'<td class="num">{effect_text(x)}</td><td class="num">{x["n_var"]:,}</td>'
            f'<td class="num">{x["match"]:.1f}</td><td>{"&#8645;" if x["sign"] == -1 else ""}</td></tr>')
    return f"""
<article class="trait">
  <header><h3>{e(s['title'])}</h3><span class="chip {chip_class(s['verdict'])}">{e(s['verdict'])}</span></header>
  {bar(s, xs)}
  <div class="scale"><span>&larr; {e(s['lower'])}</span><span>{e(s['higher'])} &rarr;</span></div>
  <p>{tendency} &middot; median {ordinal(s['pct_c'])} percentile across {len(used)} score{'s' if len(used) > 1 else ''}.</p>
  <p>{effect}</p>
  {f'<p class="note">{e(s["note"])}</p>' if s["note"] else ''}
  {f'<p class="flags">{"; ".join(flags)}.</p>' if flags else ''}
  <details><summary>All {len(xs)} score{'s' if len(xs) > 1 else ''}</summary>
  <div class="tw"><table><thead><tr><th>PGS</th><th>Reported trait</th><th>Source</th><th>Pctile*</th><th>Z*</th>
  <th>R²</th><th>Predicts</th><th>Variants</th><th>Used %</th><th>Flipped</th></tr></thead>
  <tbody>{''.join(rows)}</tbody></table></div>
  </details>
</article>"""


def notes(data, by_key):
    """Things a reader would otherwise trip over, derived from the data."""
    out = []
    for s in data["summary"]:
        xs = by_key[s["key"]]
        used = [x for x in xs if x["eur"]] or xs
        # an unvalidated score far from the consensus is the usual source of a scary percentile
        for x in used:
            if abs(x["z"] - s["zc"]) >= 1.5 and x["r2"] is None and x["or"] is None and len(used) >= 3:
                out.append(
                    f'<li><b>{e(s["title"])}: {ordinal(x["pct"])} percentile on <span class="mono">{e(x["acc"])}</span></b> '
                    f'({e(x["pub"])}) disagrees with the other {len(used) - 1} scores (consensus Z {s["zc"]:+.2f}) '
                    f'and has no published validation. Read the consensus, not this score.</li>')
        flipped = [x for x in xs if x["sign"] == -1]
        if flipped:
            ids = ", ".join(x["acc"] for x in flipped)
            out.append(
                f'<li><b>{e(s["title"])}: {len(flipped)} of {len(xs)} scores run backwards.</b> Their source coding '
                f'makes a higher raw score mean <i>{e(s["lower"])}</i> ({e(ids)}). They are flipped here, so their raw '
                f'percentiles elsewhere will look like the opposite of what is shown.</li>')
    non_eur = [x for x in data["scores"] if not x["eur"]]
    if non_eur:
        out.append(f'<li><b>{len(non_eur)} score{"s were" if len(non_eur) > 1 else " was"} trained on '
                   f'non-European samples</b> ({e(", ".join(sorted(x["acc"] for x in non_eur)[:12]))}'
                   f'{"&hellip;" if len(non_eur) > 12 else ""}). They are shown hollow and left out of each consensus, '
                   f'because a score transfers poorly across ancestries.</li>')
    m = data["meta"]
    if m["confidence"] < 0.9:
        out.append(f'<li><b>Ancestry is not a clean match to one reference population</b> '
                   f'({e(m["population"])} at {m["confidence"] * 100:.0f}%). Percentiles here come from the '
                   f'ancestry-adjusted Z rather than a single population\'s ranking, and carry more uncertainty.</li>')
    if data["excluded"]:
        low = [x for x in data["excluded"] if x["status"] == "below_min_overlap"]
        rows = ""
        for x in sorted(data["excluded"], key=lambda x: x["acc"]):
            used = "" if x["overlap"] is None else "{:.1f}".format(x["overlap"] * 100)
            reason = x["status"].replace("_", " ")
            rows += (f'<tr><td class="mono">{e(x["acc"])}</td><td>{e(x["trait"])}</td>'
                     f'<td>{e(reason)}</td><td class="num">{used}</td></tr>')
        out.append(f'<li><b>{len(data["excluded"])} score{"s are" if len(data["excluded"]) > 1 else " is"} not reported</b>'
                   f'{f" ({len(low)} because too few of their variants could be used)" if low else ""}. '
                   f'A score computed from a fraction of its variants is not the published score.'
                   f'<details><summary>Which ones</summary><div class="tw"><table><thead><tr><th>PGS</th><th>Trait</th>'
                   f'<th>Reason</th><th>Used %</th></tr></thead><tbody>{rows}</tbody></table></div></details></li>')
    return out


def build(data):
    by_key = {}
    for x in data["scores"]:
        by_key.setdefault(x["key"], []).append(x)
    m = data["meta"]

    key_find = []
    for s in data["summary"]:
        h = s["head"]
        sized = (s["kind"] == "q" and h.get("shift_sd") is not None and abs(h["shift_sd"]) >= 0.25) or \
                (s["kind"] == "b" and h.get("rel_odds") and abs(math.log(h["rel_odds"])) >= math.log(1.25))
        if sized and abs(s["zc"]) >= 1:
            direction = s["higher"] if s["zc"] > 0 else s["lower"]
            key_find.append(f'<li><b>{e(s["title"])}</b>: towards {e(direction)} &mdash; {effect_text(h)} '
                            f'<span class="chip {chip_class(s["verdict"])}">{e(s["verdict"])}</span></li>')

    sections = []
    for d in m["domains"]:
        items = [trait_card(s, by_key[s["key"]]) for s in data["summary"] if s["dom"] == d]
        if items:
            sections.append(f'<section><h2>{e(d)}</h2>{"".join(items)}</section>')

    basis = ("rank among the reference genomes of the most similar population"
             if m["percentile_basis"] == "population" else "the ancestry-adjusted Z")
    note_items = notes(data, by_key)
    sub = (f'{e(m["sample"])} &middot; {m["n_scores"]} scores across {m["n_traits"]} traits &middot; '
           f'most similar reference population: {e(m["population"])} ({m["confidence"] * 100:.0f}%) &middot; '
           f'{m["n_reference"]:,} reference genomes &middot; {e(m["generated"])}')
    body = f"""
<p class="sub">{sub}</p>

<div class="box"><h2>How to read this</h2>
<dl>
<dt>The bar</dt><dd>Every bar runs the same way: <b>right = more of the trait as named</b> (more LDL, more sleep, more of a morning person). Scores whose coding ran backwards have been flipped, so you never have to work out direction yourself.</dd>
<dt>Colour</dt><dd><span style="color:var(--good)">Green side</span> = the healthier direction, <span style="color:var(--bad)">red side</span> = the less healthy one, grey = no health direction (height, chronotype, coffee&hellip;).</dd>
<dt>&#9670; and dots</dt><dd>The diamond is the consensus (median Z of the European-trained scores). Each dot is one score; hollow dots are non-European-trained and excluded.</dd>
<dt>Percentile vs size</dt><dd>A 1st-percentile score only matters as much as the score <i>explains</i>. The predicted shift is r&nbsp;&times;&nbsp;Z, with r&nbsp;=&nbsp;&radic;R² from the score's published validation &mdash; so a weak score at an extreme percentile can predict less than a strong score near the middle.</dd>
<dt>Units</dt><dd>Converted with approximate UK Biobank population SDs where one is known. Elsewhere effects are given in population SDs: 0.1 is negligible, 0.25 small, 0.5 moderate.</dd>
</dl></div>

<div class="box"><h2>Sizeable, well-validated findings</h2>
<p class="note">Consensus at least 1 SD from average <i>and</i> a best-validated score predicting a shift of at least 0.25 SD (or odds changed by 25% or more).</p>
<ul>{''.join(key_find) or '<li>None</li>'}</ul></div>

{f'<div class="box"><h2>Worth knowing before you read on</h2><ul>{"".join(note_items)}</ul></div>' if note_items else ''}

{''.join(sections)}

<div class="limits">
<p><b>Limits.</b> A polygenic score estimates an inherited <i>tendency</i> relative to people of similar ancestry; it is not a measurement and not a diagnosis. An actual LDL, HbA1c or resting heart rate beats any score, and for most traits here the score explains under 10% of the variation between people. Most scores were built and validated in UK Biobank, so they inherit its self-reported phenotypes and its healthy-volunteer bias. *Percentile and Z are aligned to the trait as named on each card; "Flipped" (&#8645;) marks scores whose raw coding ran the other way.</p>
<p><b>Method.</b> Scores are computed from the gVCF against the HGDP+1kGP reference panel: matched by position and alleles, strand-ambiguous variants dropped, missing genotypes set to the reference mean. Z is adjusted for ancestry by regressing the score's mean and variance on the leading principal components of the reference; percentiles are from {basis}. Effect sizes come from PGS Catalog performance records, European evaluations preferred. R² values above the trait's heritability ceiling, or from a score's own training cohort, are not used for sizing. ORs are used only for binary traits and only when plausibly per-SD (1/3&ndash;3). Every number on this page is in <span class="mono">{e(m["sample"])}.pgs_interpreted.tsv</span>.</p>
</div>"""
    return body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json")
    ap.add_argument("--html")
    ap.add_argument("--fragment")
    a = ap.parse_args()
    data = json.load(open(a.json, encoding="utf-8"))
    body = build(data)
    if a.fragment:
        with open(a.fragment, "w", encoding="utf-8") as fh:
            fh.write(f'<style>{PGS_CSS}</style>\n<div class="pgs">{body}</div>\n')
    if a.html:
        with open(a.html, "w", encoding="utf-8") as fh:
            fh.write(f'<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
                     f'<meta name="viewport" content="width=device-width,initial-scale=1">'
                     f'<title>Polygenic score overview</title>\n<style>{TOKENS}{PAGE_CSS}{PGS_CSS}</style></head>'
                     f'<body><main><h1>Polygenic score overview</h1>\n<div class="pgs">{body}</div></main></body></html>\n')
    print(f"[pgs_overview] {data['meta']['n_scores']} scores, {data['meta']['n_traits']} traits")


if __name__ == "__main__":
    main()
