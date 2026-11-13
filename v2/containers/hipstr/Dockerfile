# genomeportrait HipSTR container (Docker engine variant of containers/hipstr/Singularity.def)
# Built locally by the pipeline under -profile docker, or by hand:
#   docker build -t genomeportrait-hipstr:0.7 containers/hipstr
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential g++ make git autoconf automake \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        ca-certificates procps && \
    cd /opt && git clone --recursive https://github.com/HipSTR-Tool/HipSTR.git && \
    cd HipSTR && make && install -m 0755 HipSTR /usr/local/bin/HipSTR && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV PATH=/usr/local/bin:$PATH
LABEL org.opencontainers.image.title="genomeportrait-hipstr" \
      org.opencontainers.image.version="0.7"
ENTRYPOINT ["HipSTR"]
