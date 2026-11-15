# Blood-group assets (GRCh38)

- `sites.pvar` — the variants that define ABO, Kell, Duffy and Kidd, in the form
  `pgs_engine.py target` reads (plus-strand alleles from Ensembl). See `scripts/blood_groups.py`.
- `rhd_regions.bed` — the RHD gene (ENSG00000187010) and two 58 kb control windows of
  unique sequence 1.3 Mb either side. RHD copy number = depth over RHD / depth over the
  controls, counting reads with MAPQ >= 20.
