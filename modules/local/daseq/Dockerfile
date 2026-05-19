FROM continuumio/miniconda3:26.3.2
LABEL authors="caroline.sullivan@yale.edu"
LABEL description="Docker image containing all requirements for the DA-seq module"

# Build with command: docker build . -t immcantation/daseq:<ver>

# Update system packages and install dependencies
RUN apt-get update \
    && apt-get install -y procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH=/opt/conda/envs/daseq/bin:$PATH

RUN Rscript -e "pak::pak('KlugerLab/DAseq')"
