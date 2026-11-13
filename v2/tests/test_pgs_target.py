#!/usr/bin/env python3
"""Unit test for `pgs_engine.py target`: a hand-made gVCF against a hand-made reference.

    python3 tests/test_pgs_target.py

Each reference variant below exercises one rule of the gVCF -> reference-allele-space
conversion. Exits non-zero on any mismatch.
"""
import gzip
import os
import subprocess
import sys
import tempfile

ENGINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "pgs_engine.py")

PVAR = """#CHROM\tPOS\tID\tREF\tALT
1\t100\t1:100:A:G\tA\tG
1\t150\t1:150:C:T\tC\tT
1\t200\t1:200:A:G\tA\tG
1\t200\t1:200:A:T\tA\tT
1\t300\t1:300:A:AT\tA\tAT
1\t300\t1:300:AT:A\tAT\tA
1\t400\t1:400:CTT:C\tCTT\tC
1\t401\t1:401:T:G\tT\tG
1\t500\t1:500:G:A\tG\tA
1\t600\t1:600:G:A\tG\tA
1\t700\t1:700:G:A\tG\tA
1\t800\t1:800:C:T\tC\tT
1\t900\t1:900:AT:A\tAT\tA
1\t1010\t1:1010:C:T\tC\tT
1\t1015\t1:1015:G:A\tG\tA
2\t50\t2:50:A:C\tA\tC
3\t10\t3:10:A:C\tA\tC
"""


def rec(c, p, r, a, info, fmt, smp):
    return f"chr{c}\t{p}\t.\t{r}\t{a}\t0\t.\t{info}\t{fmt}\t{smp}\n"


BLOCK, VAR = "GT:GQ:MIN_DP:PL", "GT:GQ:DP:AD:VAF:PL"
GVCF = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + "".join([
    rec(1, 90, "A", "<*>", "END=120", BLOCK, "0/0:40:20:0,1,2"),             # 100: inside a confident block
    rec(1, 150, "C", "T,<*>", ".", VAR, "0/1:30:20:10,10,0:0.5,0:0,1,2"),      # 150: heterozygous
    rec(1, 200, "A", "T,<*>", ".", VAR, "1/1:30:20:0,20,0:1,0:0,1,2"),         # 200: A>T hom; A>G gets zero copies
    rec(1, 300, "A", "AT,<*>", ".", VAR, "0/1:30:20:10,10,0:0.5,0:0,1,2"),     # 300: the insertion, not the deletion
    rec(1, 400, "CTTT", "CT,<*>", ".", VAR, "0/1:30:20:10,10,0:0.5,0:0,1,2"),  # 400: CTT>C spelled CTTT>CT
    #                                                                           401: inside that deletion, no record of its own
    rec(1, 480, "G", "<*>", "END=520", BLOCK, "0/0:5:3:0,1,2"),              # 500: thin evidence (GQ 5)
    rec(1, 590, "G", "<*>", "END=610", BLOCK, "0/0:1:0:0,0,0"),              # 600: zero depth
    rec(1, 700, "G", "<*>", "END=700", BLOCK, "./.:0:2:0,0,0"),              # 700: uncalled
    #                                                                           800: no record at all
    rec(1, 900, "A", "<*>", "END=950", BLOCK, "0/0:40:20:0,1,2"),            # 900: a deletion site inside a block
    # DeepVariant writes long RefCall deletion candidates that overlap the records after
    # them. The candidate must not swallow 1010 and 1015, which have their own records.
    rec(1, 1000, "ACGTACGTACGTACGTACGT", "A,<*>", ".", VAR, "0/0:21:38:31,7,0:0.18,0:0,23,23"),
    rec(1, 1001, "C", "<*>", "END=1009", BLOCK, "0/0:40:20:0,1,2"),
    rec(1, 1010, "C", "T,<*>", ".", VAR, "1/1:30:20:0,20,0:1,0:0,1,2"),
    rec(1, 1011, "G", "<*>", "END=1030", BLOCK, "0/0:40:20:0,1,2"),
    rec(2, 40, "A", "<*>", "END=60", BLOCK, "0/0:40:20:0,1,2"),              # chr3 is absent from the gVCF
])

EXPECT = {
    0: {"1:100:A:G": "0/0", "1:150:C:T": "0/1", "1:200:A:G": "0/0", "1:200:A:T": "1/1",
        "1:300:A:AT": "0/1", "1:300:AT:A": "0/0", "1:400:CTT:C": "0/1", "1:401:T:G": "./.",
        "1:500:G:A": "0/0", "1:600:G:A": "./.", "1:700:G:A": "./.", "1:800:C:T": "./.",
        "1:900:AT:A": "0/0", "1:1010:C:T": "1/1", "1:1015:G:A": "0/0", "2:50:A:C": "0/0", "3:10:A:C": "./."},
    10: {"1:500:G:A": "./.", "1:100:A:G": "0/0", "1:150:C:T": "0/1", "1:1015:G:A": "0/0"},
}

failed = 0
with tempfile.TemporaryDirectory() as d:
    with open(os.path.join(d, "ref.pvar"), "w") as fh:
        fh.write(PVAR)
    with gzip.open(os.path.join(d, "t.g.vcf.gz"), "wt") as fh:
        fh.write(GVCF)
    for min_gq, expect in EXPECT.items():
        subprocess.run([sys.executable, ENGINE, "target", "--gvcf", os.path.join(d, "t.g.vcf.gz"),
                        "--pvar", os.path.join(d, "ref.pvar"), "--sample", "S", "--min-gq", str(min_gq),
                        "--out", os.path.join(d, "o.vcf.gz"), "--stats", os.path.join(d, "o.tsv")],
                       check=True, stderr=subprocess.DEVNULL)
        with gzip.open(os.path.join(d, "o.vcf.gz"), "rt") as fh:
            got = {l.split("\t")[2]: l.rstrip().split("\t")[9] for l in fh if l[0] != "#"}
        bad = {k: (got.get(k), v) for k, v in expect.items() if got.get(k) != v}
        n_ref = PVAR.count("\n") - 1
        if len(got) != n_ref:
            bad["<count>"] = (len(got), n_ref)
        print(f"min_gq {min_gq}: {len(got)} sites, {'OK' if not bad else 'MISMATCH (got, expected): ' + str(bad)}")
        failed += bool(bad)
sys.exit(1 if failed else 0)
