# Chromosome painting

`holdout.tsv`: twenty 1000 Genomes genomes (2 each of GBR, YRI, CHB and GIH, 6 each of the
admixed ASW and ACB; no trio members) that `scripts/paint_reference.sh` keeps out of both the
Beagle phasing reference and the RFMix reference, so that `tests/validate_painting.sh` can
paint them without them having seen themselves. Leaving 20 of 3,202 genomes out does not
change the reference for anyone else.

The reference itself is built at run time (`PAINTING_REFERENCE`) from the 1000 Genomes 30x
phased panel (Byrska-Bishop et al. 2022) and Beagle's GRCh38 genetic maps.
