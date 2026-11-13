#!/usr/bin/env bash
# Launch a genomeportrait v2 test run from a copy of the pipeline on local disk.
#
#   SCRATCH=/fast/local/dir REFBASE=/path/to/refs bash tests/run_test.sh <name> [nextflow options...]
#
#   smoke     -profile test: FASTQ -> CRAM -> DeepVariant -> genome facts -> report (tiny data)
#   modes     the smoke test's own outputs fed back as an aligned-only and a portrait row
#             (run `smoke` first); see tests/modes_samplesheet.example.csv
#   <other>   anything else: pass --input and -c yourself, e.g. a portrait run on a real genome
#
# The pipeline is copied to $SCRATCH first because Singularity cannot use a FUSE mount
# as its working directory and Nextflow's work directory must be on local disk.
set -euo pipefail
NAME=${1:?usage: tests/run_test.sh <smoke|modes|name> [nextflow options]}; shift
SRC=$(cd "$(dirname "$0")/.." && pwd)
SCRATCH=${SCRATCH:?set SCRATCH to a local scratch directory}
REFBASE=${REFBASE:-$SCRATCH/refs}
NEXTFLOW=${NEXTFLOW:-nextflow}
PROFILE=${PROFILE:-singularity}

mkdir -p "$SCRATCH/run-$NAME" "$SCRATCH/work-$NAME"
rsync -a --delete --exclude '.nextflow*' --exclude work --exclude /local "$SRC/" "$SCRATCH/pipeline-$NAME/"
cd "$SCRATCH/run-$NAME"
extra=()
case "$NAME" in
    smoke) extra=(-profile "test,$PROFILE") ;;
    modes)
        R=$SCRATCH/run-smoke/results
        [ -d "$R" ] || { echo "run the smoke test first" >&2; exit 1; }
        sed "s#RESULTS#$R#g" "$SRC/tests/modes_samplesheet.example.csv" > samplesheet.csv
        extra=(-profile "test,$PROFILE" --input "$SCRATCH/run-$NAME/samplesheet.csv") ;;
    *) extra=(-profile "$PROFILE") ;;
esac
exec "$NEXTFLOW" run "$SCRATCH/pipeline-$NAME" "${extra[@]}" \
    -work-dir "$SCRATCH/work-$NAME" \
    --reference_base "$REFBASE" \
    --outdir "$SCRATCH/run-$NAME/results" \
    -resume "$@"
