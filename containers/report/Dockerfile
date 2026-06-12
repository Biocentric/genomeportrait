# genomeportrait report container
# Bundles the scientific-python stack used by the report/summary bin/ scripts.
# Build & push:
#   docker build -t ghcr.io/your-org/genomeportrait-report:1.0.0 containers/report
#   docker push ghcr.io/your-org/genomeportrait-report:1.0.0
# Singularity (ORAS) image used by the singularity profile:
#   apptainer build genomeportrait-report_1.0.0.sif docker://ghcr.io/your-org/genomeportrait-report:1.0.0
#   apptainer push genomeportrait-report_1.0.0.sif oras://ghcr.io/your-org/genomeportrait-report:1.0.0-singularity
FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"
LABEL org.opencontainers.image.title="genomeportrait-report" \
      org.opencontainers.image.description="pandas/matplotlib/seaborn/jinja2 + bcftools for genomeportrait reporting" \
      org.opencontainers.image.version="1.0.0"
