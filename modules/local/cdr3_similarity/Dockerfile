FROM continuumio/miniconda3:26.3.2
LABEL authors="caroline.sullivan@yale.edu"
LABEL description="Docker image containing all requirements for the cdr3_similarity module"

# Build with command: docker build . -t immcantation/cdr3similiarity:<ver>

# Update system packages and install dependencies
RUN apt-get update \
    && apt-get install -y procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/cdr3_similarity/bin:$PATH

RUN Rscript -e "pak::pak('alakazam@1.4.3')"
RUN Rscript -e "pak::pak('shazam@1.3.2')"
RUN Rscript -e "pak::pak('scoper@1.5.0')"
