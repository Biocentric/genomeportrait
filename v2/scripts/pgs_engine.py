#!/usr/bin/env python3
"""genomeportrait polygenic-score engine.

Four subcommands, run in this order (plink2 does the arithmetic between them):

    sites      scoring files            -> positions to pull from the reference panel
    match      scoring files + ref pvar -> one plink2 score file per PGS, in ALT-allele space
    target     gVCF + ref pvar          -> the sample as a VCF in the reference's allele space
    normalise  plink2 .sscore files     -> ancestry-adjusted Z and percentiles per PGS

Design notes
------------
* Everything is expressed as a weight on the reference panel's ALT allele. A scoring
  row whose effect allele is the panel REF contributes w * (2 - d_alt): the weight
  becomes -w and 2w goes into a per-score constant (`offset`). That lets several
  rows hitting the same variant be summed, and makes every score a plain dot product.
* Matching is on chr:pos plus alleles, never on position alone: a position-only
  check once reported 99.95% while no allele could match.
* Strand-ambiguous SNPs (A/T, C/G) are excluded by default, as pgsc_calc does.
* Autosomes only. chrX needs sex-aware ploidy and GWAS rarely carry it.

Only the standard library and numpy are used.
"""
import argparse
import csv
import glob
import gzip
import math
import multiprocessing as mp
import os
import subprocess
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
AUTOSOMES = {str(i) for i in range(1, 23)}


def log(*a):
    print("[pgs_engine]", *a, file=sys.stderr, flush=True)


def chrom_key(c):
    """'chr7' / '7' -> '7'; anything that is not an autosome -> None."""
    c = c[3:] if c[:3] in ("chr", "Chr", "CHR") else c
    return c if c in AUTOSOMES else None


def revcomp(a):
    try:
        return "".join(COMP[b] for b in reversed(a))
    except KeyError:
        return None


# --------------------------------------------------------------------------------------
# Scoring-file reader
# --------------------------------------------------------------------------------------
def read_scorefile(path):
    """Yield (chrom, pos, effect_allele, other_allele, weight, supported) per row.

    Harmonised files carry the author's original coordinates in chr_name /
    chr_position (often GRCh37) and the GRCh38 ones in hm_chr / hm_pos. Only the hm_
    columns are used. Rows without a usable position, allele or weight are yielded
    with chrom None so the caller can still count them in the denominator.
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        col = None
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if col is None:
                col = {k: i for i, k in enumerate(f)}
                if "hm_chr" not in col or "hm_pos" not in col:
                    raise SystemExit(f"{path}: not a harmonised scoring file (no hm_chr/hm_pos)")
                i_chr, i_pos = col["hm_chr"], col["hm_pos"]
                i_ea, i_w = col.get("effect_allele"), col.get("effect_weight")
                i_oa, i_ioa = col.get("other_allele"), col.get("hm_inferOtherAllele")
                flags = [col[k] for k in ("is_interaction", "is_recessive", "is_dominant",
                                          "is_haplotype", "is_diplotype") if k in col]
                n = len(f)
                continue
            if len(f) < n:
                f += [""] * (n - len(f))
            # Non-additive rows cannot be expressed as weight * dosage.
            if any(f[i].strip().lower() == "true" for i in flags):
                yield None, 0, "", "", 0.0, False
                continue
            c = chrom_key(f[i_chr])
            try:
                pos = int(f[i_pos])
                w = float(f[i_w])
            except (ValueError, TypeError):
                yield None, 0, "", "", 0.0, True
                continue
            ea = f[i_ea].upper() if i_ea is not None else ""
            oa = f[i_oa].upper() if i_oa is not None else ""
            if not oa and i_ioa is not None:
                # hm_inferOtherAllele may list several candidates (A/C); only use it when single
                cand = f[i_ioa].upper()
                oa = cand if cand and "/" not in cand else ""
            if c is None or not ea:
                yield None, 0, "", "", 0.0, True
                continue
            yield c, pos, ea, oa, w, True


def score_id(path):
    return os.path.basename(path).split("_hmPOS")[0].split(".")[0]


def list_scorefiles(d):
    files = sorted(glob.glob(os.path.join(d, "*_hmPOS_GRCh38.txt.gz")))
    if not files:
        raise SystemExit(f"no *_hmPOS_GRCh38.txt.gz files in {d}")
    return files


# --------------------------------------------------------------------------------------
# sites
# --------------------------------------------------------------------------------------
def _sites_one(path):
    out = {}
    for c, pos, _ea, _oa, _w, _ok in read_scorefile(path):
        if c is not None:
            out.setdefault(c, set()).add(pos)
    return out


def cmd_sites(a):
    files = list_scorefiles(a.scorefiles)
    union = {}
    with mp.Pool(a.threads) as pool:
        for n, part in enumerate(pool.imap_unordered(_sites_one, files), 1):
            for c, s in part.items():
                union.setdefault(c, set()).update(s)
            if n % 20 == 0:
                log(f"sites: {n}/{len(files)} scoring files read")
    total = 0
    with open(a.out, "w") as fh:
        for c in sorted(union, key=int):
            for p in sorted(union[c]):
                fh.write(f"{c}\t{p}\t{p}\n")     # plink2 --extract bed1: 1-based, inclusive
                total += 1
    log(f"sites: {total:,} distinct autosomal positions from {len(files)} scoring files -> {a.out}")


# --------------------------------------------------------------------------------------
# match
# --------------------------------------------------------------------------------------
REF = {}          # (chrom, pos) -> [(id, ref, alt), ...]   (module-level so forked workers share it)


def load_pvar(path, afreq=None):
    """Read a plain-text plink2 .pvar (and optionally its .afreq, same row order) into REF."""
    freqs = None
    if afreq:
        freqs = open(afreq, encoding="utf-8")
        h = freqs.readline().rstrip("\n").split("\t")
        i_id, i_af = h.index("ID"), h.index("ALT_FREQS")
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 5)
            af = 0.5
            if freqs:
                g = freqs.readline().rstrip("\n").split("\t")
                if g[i_id] != f[2]:
                    raise SystemExit(f"match: {afreq} is not in the same order as {path} ({g[i_id]} vs {f[2]})")
                try:
                    af = float(g[i_af])
                except ValueError:
                    af = 0.0
            c = chrom_key(f[0])
            if c is None:
                continue
            REF.setdefault((c, int(f[1])), []).append((f[2], f[3].upper(), f[4].rstrip("\n").upper(), af, n))
            n += 1
    if freqs:
        freqs.close()
    return n


def is_ambiguous(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and COMP.get(ref) == alt


def match_row(c, pos, ea, oa, keep_ambiguous):
    """-> (variant_id, sign, status). sign +1: effect allele is the panel ALT; -1: the REF."""
    cands = REF.get((c, pos))
    if not cands:
        return None, 0, "unmatched"
    ea_f, oa_f = revcomp(ea), (revcomp(oa) if oa else "")
    best = None
    for vid, ref, alt, af, row in cands:
        amb = is_ambiguous(ref, alt)
        if oa:
            if ea == alt and oa == ref:
                hit = ((vid, row), +1, "matched")
            elif ea == ref and oa == alt:
                hit = ((vid, row), -1, "matched")
            elif not amb and ea_f == alt and oa_f == ref:
                hit = ((vid, row), +1, "flipped")
            elif not amb and ea_f == ref and oa_f == alt:
                hit = ((vid, row), -1, "flipped")
            else:
                continue
        else:
            # No other allele published: the effect allele alone has to identify the variant.
            if ea == alt:
                hit = ((vid, row), +1, "matched")
            elif ea == ref:
                hit = ((vid, row), -1, "matched")
            elif not amb and ea_f == alt:
                hit = ((vid, row), +1, "flipped")
            elif not amb and ea_f == ref:
                hit = ((vid, row), -1, "flipped")
            else:
                continue
        if amb and not keep_ambiguous:
            hit = (None, 0, "ambiguous")
        # Several reference variants can fit one row: an allele pair such as A/AAT is both the
        # insertion A>AAT and the deletion AAT>A, and a scoring file does not say which allele
        # was the reference. Prefer a same-strand match, then the more common variant - a
        # GWAS hit is far likelier to be the common variant than a rare one at the same site.
        rank = (hit[2] != "matched", -min(af, 1.0 - af))
        if best is None or rank < best[0]:
            best = (rank, hit)
    return best[1] if best else (None, 0, "unmatched")


def _match_one(args):
    path, outdir, keep_ambiguous = args
    pid = score_id(path)
    weights, rows, offset = {}, {}, 0.0
    cnt = dict(rows=0, matched=0, flipped=0, ambiguous=0, unmatched=0, unsupported=0, merged=0)
    for c, pos, ea, oa, w, ok in read_scorefile(path):
        cnt["rows"] += 1
        if not ok:
            cnt["unsupported"] += 1
            continue
        if c is None:
            cnt["unmatched"] += 1
            continue
        hit, sign, status = match_row(c, pos, ea, oa, keep_ambiguous)
        if hit is None:
            cnt[status] += 1
            continue
        cnt["matched"] += 1
        if status == "flipped":
            cnt["flipped"] += 1
        if sign < 0:
            offset += 2.0 * w
        vid, row = hit
        if vid in weights:
            cnt["merged"] += 1
        weights[vid] = weights.get(vid, 0.0) + sign * w
        rows[vid] = row
    with open(os.path.join(outdir, f"{pid}.score"), "w") as fh:
        fh.write("ID\tA1\tW\n")
        for vid, w in weights.items():
            fh.write(f"{vid}\t{vid.rsplit(':', 1)[-1]}\t{w!r}\n")
    # The same score as row numbers into the reference .pvar plus ALT-allele weights, for
    # `score`: one sample then costs one numpy dot product per score instead of a plink2 run.
    import numpy as np
    np.savez(os.path.join(outdir, f"{pid}.npz"),
             row=np.fromiter((rows[v] for v in weights), dtype=np.int32, count=len(weights)),
             w=np.fromiter(weights.values(), dtype=np.float64, count=len(weights)))
    cnt.update(pgs_id=pid, offset=offset, n_ids=len(weights))
    return cnt


def cmd_match(a):
    files = list_scorefiles(a.scorefiles)
    n = load_pvar(a.pvar, a.afreq)
    log(f"match: {n:,} reference variants at {len(REF):,} positions")
    # The A1 column is taken from the variant ID's last field; verify the panel really
    # names variants chr:pos:ref:alt before relying on that.
    bad = sum(1 for rows in list(REF.values())[:200000] for vid, _r, alt, _af, _n in rows
              if vid.rsplit(":", 1)[-1].upper() != alt)
    if bad:
        raise SystemExit(f"match: {bad} reference variant IDs do not end in their ALT allele; "
                         "this engine expects chr:pos:ref:alt IDs")
    os.makedirs(a.outdir, exist_ok=True)
    jobs = [(f, a.outdir, a.keep_ambiguous) for f in files]
    stats = []
    with mp.Pool(a.threads) as pool:      # fork: workers share REF copy-on-write
        for k, s in enumerate(pool.imap_unordered(_match_one, jobs), 1):
            stats.append(s)
            if k % 20 == 0:
                log(f"match: {k}/{len(files)} scores")
    stats.sort(key=lambda s: s["pgs_id"])
    cols = ["pgs_id", "rows", "matched", "flipped", "ambiguous", "unmatched", "unsupported",
            "merged", "n_ids", "offset"]
    with open(a.stats, "w") as fh:
        fh.write("\t".join(cols + ["pct_matched"]) + "\n")
        for s in stats:
            pct = 100.0 * s["matched"] / s["rows"] if s["rows"] else 0.0
            fh.write("\t".join(repr(s[c]) if c == "offset" else str(s[c]) for c in cols) + f"\t{pct:.2f}\n")
    tot = sum(s["rows"] for s in stats)
    log(f"match: {len(stats)} scores, {sum(s['matched'] for s in stats):,} of {tot:,} rows matched "
        f"({sum(s['ambiguous'] for s in stats):,} ambiguous excluded, "
        f"{sum(s['flipped'] for s in stats):,} strand-flipped)")


# --------------------------------------------------------------------------------------
# target
# --------------------------------------------------------------------------------------
def trim(ref, alt):
    """Drop shared trailing bases so two spellings of an indel at one anchor compare equal."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return ref, alt


def read_tbi(path):
    """Tabix index -> {sequence name: first virtual offset of its records}.

    Only the start of each sequence is needed, so this reads the chunk lists and keeps the
    smallest begin offset per sequence (the pseudo-bin 37450 holds statistics, not data).
    """
    import struct
    d = gzip.decompress(open(path, "rb").read())
    if d[:4] != b"TBI\x01":
        raise ValueError(f"{path}: not a tabix index")
    n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack_from("<8i", d, 4)
    names = d[36:36 + l_nm].rstrip(b"\0").split(b"\0")
    o, starts = 36 + l_nm, {}
    for r in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", d, o); o += 4
        best = None
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", d, o); o += 8
            for _ in range(n_chunk):
                beg, _end = struct.unpack_from("<QQ", d, o); o += 16
                if bin_id != 37450 and (best is None or beg < best):
                    best = beg
        (n_intv,) = struct.unpack_from("<i", d, o); o += 4 + 8 * n_intv
        if best is not None:
            starts[names[r].decode()] = best
    return starts


def bgzf_lines(path, voffset):
    """Lines of a BGZF file from a virtual offset, without an external tool."""
    import struct
    import zlib
    with open(path, "rb") as fh:
        fh.seek(voffset >> 16)
        skip, tail = voffset & 0xFFFF, b""
        while True:
            head = fh.read(18)
            if len(head) < 18:
                break
            if head[:4] != b"\x1f\x8b\x08\x04":
                raise ValueError(f"{path}: not BGZF at byte {fh.tell() - 18}")
            (bsize,) = struct.unpack_from("<H", head, 16)
            block = head + fh.read(bsize + 1 - 18)
            data = zlib.decompress(block, 31)
            if skip:
                data, skip = data[skip:], 0
            if not data:
                continue
            lines = (tail + data).split(b"\n")
            tail = lines.pop()
            for line in lines:
                yield line + b"\n"
        if tail:
            yield tail


TGT = {}      # sites per chromosome and options, shared copy-on-write with forked workers


def _target_chrom(c, records, out):
    """Genotype every reference site of chromosome `c` from that chromosome's gVCF records."""
    pos_list, ids, refs, alts_ref = TGT["sites"][c]
    min_gq, trust_all = TGT["min_gq"], TGT["trust_all"]
    cnt = dict(homref=0, variant=0, other_alt=0, missing_nocov=0, missing_lowgq=0, missing_other=0)
    gq_hist = [0] * 11      # hom-ref sites from reference blocks, by GQ in steps of 5 (last bin: 50+)
    nsite, i, n_lines = len(pos_list), 0, 0

    def emit(i, gt):
        out.write(f"{c}\t{pos_list[i]}\t{ids[i]}\t{refs[i]}\t{alts_ref[i]}\t.\t.\t.\tGT\t{gt}\n")

    for raw in records:
        n_lines += 1
        if i >= nsite:
            continue
        f = raw.split(b"\t", 9)
        pos = int(f[1])
        # sites before this record were not covered by anything
        while i < nsite and pos_list[i] < pos:
            emit(i, "./."); cnt["missing_nocov"] += 1; i += 1
        if i >= nsite:
            continue
        ref = f[3].decode().upper()
        info = f[7]
        end = pos + len(ref) - 1
        k = info.find(b"END=")
        if k != -1:
            j = info.find(b";", k)
            end = int(info[k + 4:j if j != -1 else len(info)])
        if pos_list[i] > end:
            continue

        fmt = f[8].split(b":")
        smp = f[9].rstrip(b"\n").split(b":")
        gt = smp[0].decode().replace("|", "/")
        gq = None      # a caller that writes no GQ is not penalised by --min-gq
        dp = -1
        for key, val in zip(fmt, smp):
            if key == b"GQ":
                gq = int(val) if val.isdigit() else None
            elif key in (b"MIN_DP", b"DP"):
                dp = int(val) if val.isdigit() else -1
        alts = [x for x in f[4].decode().upper().split(",") if x not in ("<*>", "<NON_REF>", ".")]
        called = "." not in gt
        try:
            alleles = [int(x) for x in gt.split("/")] if called else []
        except ValueError:
            called, alleles = False, []

        if not alts:
            # reference block: every site inside it is hom-ref, if the block is trustworthy
            if not called or (dp == 0 and not trust_all):
                tag, g = "missing_nocov", "./."
            elif any(alleles):
                tag, g = "missing_other", "./."
            elif gq is not None and gq < min_gq:
                tag, g = "missing_lowgq", "./."
            else:
                tag, g = "homref", "0/0"
            n0 = i
            while i < nsite and pos_list[i] <= end:
                emit(i, g); cnt[tag] += 1; i += 1
            if tag == "homref" and gq is not None:
                gq_hist[min(gq // 5, 10)] += i - n0
        else:
            # A variant record speaks only for its own position. DeepVariant writes long
            # deletion candidates (often RefCall) that overlap the records after them; the
            # positions inside such a span have records of their own further on, so they
            # must not be consumed here.
            norm = [trim(ref, x) for x in alts]
            while i < nsite and pos_list[i] == pos:
                if not called or len(alleles) != 2 or max(alleles) > len(alts):
                    emit(i, "./."); cnt["missing_other"] += 1; i += 1
                    continue
                sref, salt = refs[i], alts_ref[i]
                d = 0
                if len(sref) == 1 and len(salt) == 1:
                    # SNP site: which base does each called allele carry at this position?
                    for x in alleles:
                        if x == 0:
                            base = ref[0]
                        else:
                            alt = alts[x - 1]
                            if len(alt) == 1 and len(ref) == 1:
                                base = alt
                            elif alt[:1] == ref[:1]:
                                base = ref[0]          # indel anchored on an unchanged base
                            else:
                                base = None            # complex allele: cannot say
                        if base is None:
                            d = None
                            break
                        d += base == salt
                else:
                    # indel site: count alleles that spell the same change at the same anchor
                    want = trim(sref, salt)
                    d = sum(1 for x in alleles if x and norm[x - 1] == want)
                if d is None:
                    emit(i, "./."); cnt["missing_other"] += 1
                elif d == 0 and not any(alleles) and gq is not None and gq < min_gq:
                    emit(i, "./."); cnt["missing_lowgq"] += 1
                else:
                    emit(i, ("0/0", "0/1", "1/1")[d])
                    cnt["variant" if d else ("other_alt" if any(alleles) else "homref")] += 1
                i += 1
    # everything not reached was never covered (also: a chromosome absent from the gVCF)
    while i < nsite:
        emit(i, "./."); cnt["missing_nocov"] += 1; i += 1
    return cnt, gq_hist, n_lines


def _target_worker(job):
    """One chromosome, read straight from its offset in the BGZF gVCF, into its own part file."""
    c, name, voffset, gvcf, part = job
    def records():
        for line in bgzf_lines(gvcf, voffset):
            if line.split(b"\t", 1)[0].decode() != name:
                return
            yield line
    with gzip.open(part, "wt", compresslevel=1) as out:
        recs = records() if voffset is not None else iter(())
        return (c,) + _target_chrom(c, recs, out)


def cmd_target(a):
    """Emit, for every reference variant, the genotype the gVCF implies.

      reference block covers the position, GT 0/0 -> 0/0
      variant record carries the panel ALT         -> its allele count
      variant record carries other alleles only    -> 0/0 (zero copies of the panel ALT)
      no record, no depth, or GQ below --min-gq     -> ./.  (mean-imputed downstream)

    A hom-ref call is only as good as the reads behind it. A block with GQ 1 and depth 0
    says nothing, and counting it as hom-ref pulls every score towards the reference
    allele; such sites are better left missing.

    With a tabix index beside the gVCF, chromosomes are processed in parallel straight
    from their offsets; without one, the file is read once from start to end.
    Output order is chromosome 1-22, then position, as in the reference .pvar.
    """
    sites = {}
    with open(a.pvar, encoding="utf-8") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 5)
            c = chrom_key(f[0])
            if c is None:
                continue
            s = sites.setdefault(c, ([], [], [], []))
            s[0].append(int(f[1])); s[1].append(f[2]); s[2].append(f[3].upper())
            s[3].append(f[4].rstrip("\n").upper())
    n_sites = sum(len(s[0]) for s in sites.values())
    log(f"target: {n_sites:,} reference variants to genotype")
    TGT.update(sites=sites, min_gq=a.min_gq, trust_all=a.trust_all_blocks)
    order = sorted(sites, key=int)
    parts = {c: f"{a.out}.part{c}" for c in order}
    results = {}

    tbi = a.gvcf + ".tbi"
    starts = None
    if os.path.exists(tbi) and a.threads > 1:
        try:
            starts = read_tbi(tbi)
        except Exception as e:  # unreadable index: fall back to a sequential read
            log(f"target: cannot use {tbi} ({e}); reading sequentially")
    if starts is not None:
        name_of = {chrom_key(n): n for n in starts if chrom_key(n)}
        jobs = [(c, name_of.get(c), starts.get(name_of.get(c)), a.gvcf, parts[c]) for c in order]
        jobs.sort(key=lambda j: -len(sites[j[0]][0]))          # longest first
        with mp.Pool(min(a.threads, len(jobs))) as pool:
            for c, cnt, hist, n in pool.imap_unordered(_target_worker, jobs):
                results[c] = (cnt, hist, n)
        log(f"target: {len(jobs)} chromosomes in parallel from the tabix index")
    else:
        import itertools
        proc = subprocess.Popen(["gzip", "-dc", a.gvcf], stdout=subprocess.PIPE, bufsize=1 << 22)
        body = (l for l in proc.stdout if l[:1] != b"#")
        seen = set()
        for name, group in itertools.groupby(body, key=lambda l: l.split(b"\t", 1)[0]):
            c = chrom_key(name.decode())
            if c is None or c not in sites:
                for _ in group:
                    pass
                continue
            if c in seen:
                raise SystemExit(f"target: gVCF is not sorted (chromosome {c} appears twice)")
            seen.add(c)
            with gzip.open(parts[c], "wt", compresslevel=1) as out:
                results[c] = _target_chrom(c, group, out)
        proc.wait()
        if proc.returncode != 0:
            raise SystemExit(f"target: gzip exited {proc.returncode} reading {a.gvcf}")
        for c in order:
            if c not in results:                          # chromosome absent from the gVCF
                with gzip.open(parts[c], "wt", compresslevel=1) as out:
                    results[c] = _target_chrom(c, iter(()), out)

    # header as one gzip member, then the chromosome parts in order (a valid multi-member gzip)
    with open(a.out, "wb") as dst:
        head = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
        head += [f"##contig=<ID={c}>" for c in order]
        head.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a.sample}")
        dst.write(gzip.compress(("\n".join(head) + "\n").encode(), compresslevel=1))
        for c in order:
            with open(parts[c], "rb") as src:
                while True:
                    buf = src.read(1 << 22)
                    if not buf:
                        break
                    dst.write(buf)
            os.remove(parts[c])

    cnt = dict.fromkeys(results[order[0]][0], 0)
    gq_hist, n_lines = [0] * 11, 0
    for c in order:
        cc, hh, nn = results[c]
        for k, v in cc.items():
            cnt[k] += v
        gq_hist = [x + y for x, y in zip(gq_hist, hh)]
        n_lines += nn
    genotyped = cnt["homref"] + cnt["variant"] + cnt["other_alt"]
    with open(a.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"reference_variants\t{n_sites}\n")
        fh.write(f"genotyped\t{genotyped}\n")
        fh.write(f"pct_genotyped\t{100.0 * genotyped / n_sites:.3f}\n")
        fh.write(f"min_gq\t{a.min_gq}\n")
        for k, v in cnt.items():
            fh.write(f"{k}\t{v}\n")
        for b, v in enumerate(gq_hist):
            fh.write(f"homref_block_sites_gq_{b * 5}{'_plus' if b == 10 else '_to_' + str(b * 5 + 4)}\t{v}\n")
    log(f"target: {n_lines:,} gVCF records read; {genotyped:,} of {n_sites:,} "
        f"({100.0 * genotyped / n_sites:.2f}%) genotyped: " + ", ".join(f"{k}={v:,}" for k, v in cnt.items()))


# --------------------------------------------------------------------------------------
# score
# --------------------------------------------------------------------------------------
def read_alt_freqs(path, n_expected=None):
    """plink2 .afreq -> numpy array of ALT frequencies, in file order."""
    import numpy as np
    with open(path) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        i_af = h.index("ALT_FREQS")
        f = np.array([float(l.split("\t")[i_af]) if l.split("\t")[i_af] not in ("NA", "nan") else 0.0
                      for l in fh], dtype=np.float64)
    if n_expected is not None and len(f) != n_expected:
        raise SystemExit(f"score: {path} has {len(f)} variants, the reference has {n_expected}")
    return f


def cmd_score(a):
    """Score one sample against every matched score in a single pass over its genotypes.

    The target VCF from `target` lists the reference variants in .pvar order. Each is read
    once into a dosage vector (ALT allele count; missing calls take 2 x the ALT frequency
    given by --impute-afreq, which is what plink2 --score --read-freq does). Every score is
    then a dot product over its <pgs>.npz rows. Output: one plink2-style .sscore per score.
    """
    import numpy as np
    ids = []
    with open(a.pvar) as fh:
        for line in fh:
            if line[0] != "#":
                f = line.split("\t", 3)
                if chrom_key(f[0]) is not None:
                    ids.append(f[2])
    n = len(ids)
    dos = np.full(n, -1, dtype=np.int8)
    sample, k, index = None, 0, None
    proc = subprocess.Popen(["gzip", "-dc", a.target], stdout=subprocess.PIPE, bufsize=1 << 22)
    for raw in proc.stdout:
        if raw[:1] == b"#":
            if raw.startswith(b"#CHROM"):
                sample = raw.rstrip(b"\n").split(b"\t")[9].decode()
            continue
        f = raw.split(b"\t", 10)
        vid = f[2].decode()
        if index is None and k < n and vid == ids[k]:
            row = k                          # the usual case: a target made by `target`
        else:
            if index is None:                # any other VCF: look variants up by ID
                index = {v: i for i, v in enumerate(ids)}
                log("score: target is not in reference order; matching variants by ID")
            row = index.get(vid)
        k += 1
        if row is None:
            continue
        gt = f[9][:3]
        if gt[:1] != b"." and gt[2:3] != b".":
            dos[row] = (gt[:1] == b"1") + (gt[2:3] == b"1")
    proc.wait()
    if proc.returncode != 0:
        raise SystemExit(f"score: gzip exited {proc.returncode} reading {a.target}")
    called = dos >= 0
    source, freq_path = "pooled", a.impute_afreq
    if a.impute != "pooled" and a.population and a.pop_afreq_dir:
        with open(a.population) as fh:
            next(fh)
            top, prob = fh.readline().rstrip("\n").split("\t")[:2]
        cand = os.path.join(a.pop_afreq_dir, f"{top}.afreq")
        if not os.path.exists(cand):
            log(f"score: no frequencies for {top} in {a.pop_afreq_dir}; using the pooled panel")
        elif a.impute == "population" or float(prob) >= a.min_confidence:
            source, freq_path = top, cand
        else:
            log(f"score: most similar population {top} at {100 * float(prob):.0f}% is below "
                f"{100 * a.min_confidence:.0f}%; using the pooled panel")
    freq = read_alt_freqs(freq_path, n)
    d = np.where(called, dos, 2.0 * freq).astype(np.float64)
    log(f"score: {int(called.sum()):,} of {n:,} reference variants called; missing ones filled from "
        f"the {source} allele frequencies ({os.path.basename(freq_path)})")
    if a.imputation_out:
        with open(a.imputation_out, "w") as fh:
            fh.write("frequencies\tfile\tcalled\treference_variants\n")
            fh.write(f"{source}\t{os.path.basename(freq_path)}\t{int(called.sum())}\t{n}\n")

    os.makedirs(a.outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.scores, "*.npz")))
    if not files:
        raise SystemExit(f"score: no .npz score files in {a.scores} (rerun `match`)")
    for path in files:
        pid = os.path.basename(path)[:-4]
        z = np.load(path)
        row, w = z["row"], z["w"]
        if len(row) == 0:
            continue
        total = float(w @ d[row])
        c = called[row]
        allele_ct = 2 * int(c.sum())
        named = float(dos[row][c].sum())
        with open(os.path.join(a.outdir, f"{pid}.sscore"), "w") as fh:
            fh.write("#IID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n")
            avg = total / allele_ct if allele_ct else float("nan")
            fh.write(f"{sample}\t{allele_ct}\t{named:g}\t{avg!r}\t{total!r}\n")
    log(f"score: {len(files)} scores written to {a.outdir}")


# --------------------------------------------------------------------------------------
# normalise
# --------------------------------------------------------------------------------------
def read_sscore(path):
    """plink2 .sscore -> (header list, {IID: row list})."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
        iid = hdr.index("IID")
        return hdr, {r[iid]: r for r in (l.rstrip("\n").split("\t") for l in fh) if len(r) > iid}


def fit_mean_variance(A, y, iters=200, tol=1e-9):
    """Heteroscedastic normal model  y ~ N(A @ beta, exp(A @ gamma)),  by maximum likelihood.

    This is the "full likelihood" adjustment pgsc_calc applies for Z_norm2 (after Khan et
    al. 2022): the mean AND the variance of a score are allowed to change along the
    principal components, so a score is not judged against a spread that belongs to
    another ancestry. Fitted by alternating two exact steps until they stop moving:
    weighted least squares for beta given the variances, and Newton steps for gamma given
    the residuals (the score equation is that of a gamma GLM with log link).
    """
    import numpy as np
    n, p = A.shape
    beta = np.linalg.lstsq(A, y, rcond=None)[0]
    gamma = np.zeros(p)
    gamma[0] = math.log(max(((y - A @ beta) ** 2).mean(), 1e-300))
    for _ in range(iters):
        w = np.exp(-(A @ gamma))
        Aw = A * w[:, None]
        beta_new = np.linalg.solve(A.T @ Aw, Aw.T @ y)
        r2 = (y - A @ beta_new) ** 2
        g = gamma.copy()
        for _ in range(50):
            ratio = r2 * np.exp(-(A @ g))
            step = np.linalg.solve((A * ratio[:, None]).T @ A, A.T @ (ratio - 1.0))
            # A full Newton step can overshoot when a few genomes sit far from the fit
            # (small scores with heavy tails); cap it so exp() stays finite.
            biggest = np.abs(step).max()
            if biggest > 1.0:
                step /= biggest
            g += step
            if biggest < 1e-12:
                break
        moved = max(np.abs(beta_new - beta).max() / (np.abs(beta).max() + 1e-300), np.abs(g - gamma).max())
        beta, gamma = beta_new, g
        if moved < tol:
            break
    if not (np.all(np.isfinite(beta)) and np.all(np.isfinite(gamma))):
        # No stable joint fit: fall back to the mean model with one common variance,
        # which makes Z_norm2 equal Z_norm1 for this score.
        beta = np.linalg.lstsq(A, y, rcond=None)[0]
        gamma = np.zeros(p)
        gamma[0] = math.log(max(((y - A @ beta) ** 2).mean(), 1e-300))
    return beta, gamma


def classify_population(ref_pcs, target_pcs, psam, pop_col="SuperPop", k=5):
    """Gaussian classifier on the leading PCs, one mean and covariance per reference population.

    Returns (best, {population: posterior probability}, labels, ref_ids, Rz, Tz) where Rz/Tz
    are the reference and target PCs standardised to the reference.
    """
    import numpy as np
    pop = {}
    with open(psam) as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
        ii, ip = hdr.index("IID"), hdr.index(pop_col)
        for l in fh:
            r = l.rstrip("\n").split("\t")
            pop[r[ii]] = r[ip]
    h, ref_pc = read_sscore(ref_pcs)
    pc_cols = [i for i, c in enumerate(h) if c.startswith("PC") and c.endswith("_AVG")]
    if len(pc_cols) < k:
        raise SystemExit(f"only {len(pc_cols)} PCs in {ref_pcs}")
    ref_ids = sorted(i for i in ref_pc if i in pop)
    R = np.array([[float(ref_pc[i][c]) for c in pc_cols] for i in ref_ids])
    h2, tgt_pc = read_sscore(target_pcs)
    pc_cols2 = [i for i, c in enumerate(h2) if c.startswith("PC") and c.endswith("_AVG")]
    tgt_id = next(iter(tgt_pc))
    T = np.array([float(tgt_pc[tgt_id][c]) for c in pc_cols2])
    mu, sd = R.mean(0), R.std(0)
    Rz, Tz = (R - mu) / sd, (T - mu) / sd
    labels = np.array([pop[i] for i in ref_ids])
    logp = {}
    for g in sorted(set(labels)):
        X = Rz[labels == g, :k]
        if len(X) < 20:
            continue
        m, cov = X.mean(0), np.cov(X.T) + np.eye(k) * 1e-6
        d = Tz[:k] - m
        logp[g] = -0.5 * (d @ np.linalg.solve(cov, d)) - 0.5 * np.linalg.slogdet(cov)[1]
    best = max(logp, key=logp.get)
    mx = max(logp.values())
    tot = sum(math.exp(v - mx) for v in logp.values())
    conf = {g: math.exp(v - mx) / tot for g, v in logp.items()}
    return best, conf, labels, ref_ids, Rz, Tz, tgt_id


def write_population(path, conf, labels, Tz):
    with open(path, "w") as fh:
        fh.write("population\tprobability\tn_reference\n")
        for g in sorted(conf, key=conf.get, reverse=True):
            fh.write(f"{g}\t{conf[g]:.6f}\t{int((labels == g).sum())}\n")
        fh.write("# target PCs (standardised to the reference): " + " ".join(f"{v:.3f}" for v in Tz) + "\n")


def cmd_population(a):
    best, conf, labels, _ids, _R, Tz, _t = classify_population(a.ref_pcs, a.target_pcs, a.psam,
                                                               a.pop_col, a.n_pcs_pop)
    write_population(a.out, conf, labels, Tz)
    log(f"population: most similar {best} ({100 * conf[best]:.1f}%)")


def cmd_normalise(a):
    import numpy as np

    best, conf, labels, ref_ids, Rz, Tz, tgt_id = classify_population(
        a.ref_pcs, a.target_pcs, a.psam, a.pop_col, a.n_pcs_pop)
    k = a.n_pcs_pop
    in_pop = labels == best
    log(f"normalise: most similar population {best} ({100 * conf[best]:.1f}%), "
        f"{int(in_pop.sum())} reference samples")

    # --- per-score adjustment ------------------------------------------------------
    stats = {}
    with open(a.match_stats) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats[r["pgs_id"]] = r
    kn = a.n_pcs_norm
    A = np.hstack([np.ones((len(ref_ids), 1)), Rz[:, :kn]])
    a_t = np.concatenate([[1.0], Tz[:kn]])

    rows = []
    for pid in sorted(stats):
        st = stats[pid]
        rf, tf = os.path.join(a.ref_scores, pid + ".sscore"), os.path.join(a.target_scores, pid + ".sscore")
        base = dict(pgs_id=pid, n_variants=int(st["rows"]), n_matched_reference=int(st["matched"]))
        if not (os.path.exists(rf) and os.path.exists(tf)):
            rows.append(dict(base, status="not_scored"))
            continue
        hr, ref_sc = read_sscore(rf)
        ht, tgt_sc = read_sscore(tf)
        cs = next(i for i, c in enumerate(hr) if c.endswith("_SUM") and "DOSAGE" not in c)
        ct = next(i for i, c in enumerate(ht) if c.endswith("_SUM") and "DOSAGE" not in c)
        ca = ht.index("ALLELE_CT")
        off = float(st["offset"])
        y = np.array([float(ref_sc[i][cs]) for i in ref_ids]) + off
        trow = tgt_sc[tgt_id]
        s = float(trow[ct]) + off
        n_geno = int(trow[ca]) // 2
        overlap = n_geno / base["n_variants"] if base["n_variants"] else 0.0
        base.update(n_genotyped=n_geno, overlap=round(overlap, 4), SUM=s)
        if overlap < a.min_overlap:
            rows.append(dict(base, status="below_min_overlap"))
            continue
        if y.std() == 0:
            rows.append(dict(base, status="no_variance"))
            continue
        beta = np.linalg.lstsq(A, y, rcond=None)[0]
        resid = y - A @ beta
        r_t = s - a_t @ beta
        z1 = r_t / resid.std(ddof=1)
        b2, g2 = fit_mean_variance(A, y)
        z2 = (s - a_t @ b2) / math.sqrt(math.exp(a_t @ g2))
        yp = y[in_pop]
        zpop = (s - yp.mean()) / yp.std()
        # Reference sums come from plink2 rounded to 6 significant digits; a sample's sum can
        # differ from a tied block of them in the last digits only. Count those as ties.
        tol = 1e-6 * y.std()
        pct_pop = 100.0 * (np.sum(yp < s - tol) + 0.5 * np.sum(np.abs(yp - s) <= tol)) / len(yp)
        pct_z2 = 50.0 * (1.0 + math.erf(z2 / math.sqrt(2.0)))
        rows.append(dict(base, status="ok", Z_MostSimilarPop=zpop, Z_norm1=z1, Z_norm2=z2,
                         percentile_MostSimilarPop=pct_pop, percentile_norm2=pct_z2))

    cols = ["pgs_id", "status", "n_variants", "n_matched_reference", "n_genotyped", "overlap", "SUM",
            "Z_MostSimilarPop", "Z_norm1", "Z_norm2", "percentile_MostSimilarPop", "percentile_norm2"]
    with open(a.out, "w") as fh:
        fh.write(f"# sample={tgt_id} most_similar_population={best} confidence={conf[best]:.4f} "
                 f"n_reference={len(ref_ids)} n_reference_in_population={int(in_pop.sum())} "
                 f"pcs_normalisation={kn} pcs_population={k}\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                "" if r.get(c) is None else (f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]))
                for c in cols) + "\n")
    write_population(a.pop_out, conf, labels, Tz)
    ok = sum(r["status"] == "ok" for r in rows)
    log(f"normalise: {ok} of {len(rows)} scores passed (min overlap {a.min_overlap})")


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("sites")
    p.add_argument("--scorefiles", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(fn=cmd_sites)

    p = sub.add_parser("match")
    p.add_argument("--scorefiles", required=True)
    p.add_argument("--pvar", required=True, help="plain-text .pvar of the reference subset")
    p.add_argument("--afreq", help="plink2 .afreq of the same fileset; breaks ties between candidate variants")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stats", required=True)
    p.add_argument("--keep-ambiguous", action="store_true")
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(fn=cmd_match)

    p = sub.add_parser("target")
    p.add_argument("--gvcf", required=True)
    p.add_argument("--pvar", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--min-gq", type=int, default=0,
                   help="hom-ref calls below this GQ become missing (0 keeps everything with depth)")
    p.add_argument("--trust-all-blocks", action="store_true",
                   help="comparison only: count every 0/0 block as hom-ref even at zero depth, "
                        "as the 2026-08 pgsc_calc target was built")
    p.add_argument("--threads", type=int, default=1,
                   help="chromosomes processed in parallel; needs a .tbi beside the gVCF")
    p.add_argument("--out", required=True)
    p.add_argument("--stats", required=True)
    p.set_defaults(fn=cmd_target)

    p = sub.add_parser("score")
    p.add_argument("--target", required=True, help="VCF from `target`")
    p.add_argument("--pvar", required=True, help="the reference .pvar the target was made from")
    p.add_argument("--scores", required=True, help="directory of <pgs>.npz from `match`")
    p.add_argument("--impute-afreq", required=True,
                   help="pooled reference .afreq; fills missing calls unless a population is chosen")
    p.add_argument("--pop-afreq-dir", help="directory of <POPULATION>.afreq for the reference populations")
    p.add_argument("--population", help="output of `population` for this sample")
    p.add_argument("--impute", choices=["auto", "population", "pooled"], default="auto",
                   help="auto: the most similar population's frequencies when its probability is at "
                        "least --min-confidence, else pooled")
    p.add_argument("--min-confidence", type=float, default=0.9)
    p.add_argument("--imputation-out", help="write which frequencies were used")
    p.add_argument("--outdir", required=True)
    p.set_defaults(fn=cmd_score)

    p = sub.add_parser("population")
    p.add_argument("--ref-pcs", required=True)
    p.add_argument("--target-pcs", required=True)
    p.add_argument("--psam", required=True)
    p.add_argument("--pop-col", default="SuperPop")
    p.add_argument("--n-pcs-pop", type=int, default=5)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_population)

    p = sub.add_parser("normalise")
    p.add_argument("--ref-scores", required=True)
    p.add_argument("--target-scores", required=True)
    p.add_argument("--ref-pcs", required=True)
    p.add_argument("--target-pcs", required=True)
    p.add_argument("--psam", required=True)
    p.add_argument("--match-stats", required=True)
    p.add_argument("--pop-col", default="SuperPop")
    p.add_argument("--n-pcs-norm", type=int, default=4)
    p.add_argument("--n-pcs-pop", type=int, default=5)
    p.add_argument("--min-overlap", type=float, default=0.75)
    p.add_argument("--out", required=True)
    p.add_argument("--pop-out", required=True)
    p.set_defaults(fn=cmd_normalise)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
