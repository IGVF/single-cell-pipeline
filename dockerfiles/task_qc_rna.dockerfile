############################################################
# Dockerfile for IGVF single-cell pipeline QC RNA task
# Based on Debian slim
############################################################

FROM r-base@sha256:4f0c3a5f1681e03af47991a48d5f8c9db4d8e141a5b77b8a33526a3223a2ea6c

LABEL maintainer = "Siddarth Wekhande"
LABEL software = "IGVF single-cell pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="IGVF"
LABEL software.version.is-production="No"
LABEL software.task="qc_rna"

ENV R_VERSION=4.1.2

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    samtools \
    python3 \
    python3-dev \
    python3-full \
    python3-pip &&\
    rm -rf /var/lib/apt/lists/*

# Install python packages
RUN pip install --quiet --no-cache-dir --break-system-packages scanpy anndata pandas

# Create and setup new user
#ENV USER=shareseq
#WORKDIR /home/$USER
#RUN groupadd -r $USER &&\
 #   useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
  #  chown $USER:$USER /home/$USER

ENV R_LIBS_USER=/usr/local/lib/R

# Copy scripts
#COPY --chown=$USER:$USER src/python/qc_rna_extract_metrics.py /usr/local/bin/
#COPY --chown=$USER:$USER src/R/barcode_rank_functions.R /usr/local/bin/
#COPY --chown=$USER:$USER src/R/rna_qc_plots.R /usr/local/bin/
#COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

COPY src/python/qc_rna_extract_metrics.py /usr/local/bin/
COPY src/R/barcode_rank_functions.R /usr/local/bin/
COPY src/R/rna_qc_plots.R /usr/local/bin/
COPY src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
