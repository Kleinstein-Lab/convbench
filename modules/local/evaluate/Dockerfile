# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:26.3.2
LABEL authors="caroline.sullivan@yale.edu"
LABEL description="Docker image containing all requirements for the evaluate module"

# Build with command (the token needs read access to ggabernet/dejabu):
# GITHUB_TOKEN=$(gh auth token) docker build --secret id=github_token,env=GITHUB_TOKEN . -t cfsullivan16/evaluate:<ver>
# Add --build-arg DEJABU_REF=<tag or commit> to install a specific version of dejabu

# Update system packages and install dependencies
RUN apt-get update \
    && apt-get install -y procps git \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Install dejabu (private repository) into the evaluate env. The token is mounted as a
# build secret so it is not stored in the image layers or metadata
ARG DEJABU_REF=main
RUN --mount=type=secret,id=github_token,required=true \
    /opt/conda/envs/evaluate/bin/pip install --no-cache-dir \
    "git+https://x-access-token:$(cat /run/secrets/github_token)@github.com/ggabernet/dejabu.git@${DEJABU_REF}"

ENV PATH=/opt/conda/envs/evaluate/bin:$PATH
# Point reticulate at the evaluate env's Python
ENV RETICULATE_PYTHON=/opt/conda/envs/evaluate/bin/python
