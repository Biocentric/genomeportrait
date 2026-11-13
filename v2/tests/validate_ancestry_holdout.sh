#!/usr/bin/env bash
# Does supervised ADMIXTURE on the ancestry SNP panel recover the known ancestry of
# reference genomes it was not trained on?
#
#   ADMIXTURE="singularity exec admixture.sif admixture" bash tests/validate_ancestry_holdout.sh \
#       <panel_anc prefix> <reference .psam> <workdir> [label column, default SuperPop]
#
# A handful of genomes per population are labelled '-' (unknown) instead of their
# population, so they take no part in defining the ancestral populations, and their
# estimated proportions are printed next to their known population. Unadmixed populations
# should come back almost entirely as their own super-population; admixed ones (ASW, ACB,
# MXL, PUR, CLM) as their known mixtures.
set -euo pipefail
PANEL=$1; PSAM=$2; WORK=$3; COL=${4:-SuperPop}
ADMIXTURE=${ADMIXTURE:?set ADMIXTURE to an admixture command}
HOLD=${HOLD:-"GBR:2 FIN:1 TSI:1 YRI:2 LWK:1 CHB:2 JPT:1 PJL:1 BEB:1 MXL:2 PUR:1 CLM:1 ASW:2 ACB:1 Druze:1 Bedouin:1"}
THREADS=${THREADS:-16}

mkdir -p "$WORK"
cp "$PANEL.bed" "$WORK/ho.bed"; cp "$PANEL.bim" "$WORK/ho.bim"; cp "$PANEL.fam" "$WORK/ho.fam"
cd "$WORK"
python3 - "$PSAM" "$COL" "$HOLD" ho.fam > ho.pop <<'EOF'
import sys
psam, col, hold, fam = sys.argv[1:5]
want = {p: int(n) for p, n in (x.split(":") for x in hold.split())}
rows = [l.rstrip("\n").split("\t") for l in open(psam)]
h = [c.lstrip("#") for c in rows[0]]
iid, sup, pop = h.index("IID"), h.index(col), h.index("Population")
lab = {r[iid]: (r[sup], r[pop]) for r in rows[1:]}
taken, held = {}, {}
out = []
for l in open(fam):
    s = l.split()[1]
    sp, p = lab.get(s, ("-", "-"))
    if p in want and taken.get(p, 0) < want[p]:
        taken[p] = taken.get(p, 0) + 1
        held[s] = (p, sp)
        out.append("-")
    else:
        out.append(sp)
print("\n".join(out))
with open("held.tsv", "w") as fh:
    for s, (p, sp) in held.items():
        fh.write(f"{s}\t{p}\t{sp}\n")
EOF
K=$(grep -vx -- - ho.pop | sort -u | wc -l)
echo "held out: $(wc -l < held.tsv) genomes; K=$K" >&2
$ADMIXTURE --supervised -j"$THREADS" ho.bed "$K" > admixture.log 2>&1
python3 - ho.fam ho.pop "ho.$K.Q" held.tsv <<'EOF'
import sys
fam, popf, qf, heldf = sys.argv[1:5]
ids = [l.split()[1] for l in open(fam)]
pops = [l.strip() for l in open(popf)]
Q = [list(map(float, l.split())) for l in open(qf)]
order = []                     # supervised ADMIXTURE: columns in order of first appearance
for p in pops:
    if p != "-" and p not in order:
        order.append(p)
held = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(heldf)}
print("sample\tpopulation\tsuperpop\t" + "\t".join(order))
for i, s in enumerate(ids):
    if s in held:
        print(f"{s}\t{held[s][0]}\t{held[s][1]}\t" + "\t".join(f"{q:.3f}" for q in Q[i]))
EOF
