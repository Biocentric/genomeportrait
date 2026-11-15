# Archaic introgression reference

`1000g_reference.tsv`: archaic sequence per 1000 Genomes population (median, 10th and 90th
percentile, Mb with at least one archaic copy, mean posterior >= 0.8), the yardstick the
report holds a genome against. Built from Laurits Skov's hmmix calls for the 1000 Genomes
high-coverage genomes (Zenodo record 14136628, CC-BY 4.0) by

    python3 assets/archaic/make_reference.py hg38_1000g_segments.txt > assets/archaic/1000g_reference.tsv

The published calls are per phased haplotype; merging a person's two haplotypes gives what the
pipeline's unphased run reports. On eight genomes from eight populations the two agree to
r = 0.99 (`tests/validate_archaic.sh`). African populations are not in the table: they are
hmmix's outgroup.
