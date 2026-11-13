#!/usr/bin/env python3
"""Check a results folder against the v2 publishing layout.

    python3 check_results_tree.py <outdir> <sample> [<sample> ...]

Every file must sit under <outdir>/<sample>/, <outdir>/multiqc/ or <outdir>/pipeline_info/;
no versions.yml; nothing that looks like a reference, index cache or intermediate.
Prints the tree summary and exits 1 on any violation.
"""
import os
import sys

outdir, samples = sys.argv[1], set(sys.argv[2:])
bad, counts = [], {}
forbidden = ('.fastq.gz', '.fq.gz', '.markdup.bam', '.sif', '.img', '.bwt', '.pac', '.0123',
             '.pgs_target.vcf.gz', '.pgs_fragment.html', 'Homo_sapiens_assembly38')
forbidden_names = ('versions.yml',)
for root, _dirs, files in os.walk(outdir):
    for f in files:
        rel = os.path.relpath(os.path.join(root, f), outdir)
        top = rel.split(os.sep)[0]
        if top not in samples | {'multiqc', 'pipeline_info'}:
            bad.append(f"outside the layout: {rel}")
        elif any(x in f for x in forbidden) or f in forbidden_names:
            bad.append(f"should not be published: {rel}")
        key = os.sep.join(rel.split(os.sep)[:2]) if top in samples else top
        counts[key] = counts.get(key, 0) + 1
for k in sorted(counts):
    print(f"{counts[k]:5d}  {k}")
for b in bad:
    print("FAIL", b)
print("OK" if not bad else f"{len(bad)} violations")
sys.exit(1 if bad else 0)
