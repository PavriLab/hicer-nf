FROM continuumio/miniconda:4.7.12

WORKDIR /hicer-nf

COPY environment.yml /hicer-nf/environment.yml


RUN apt-get update \
    && apt-get install -y procps \
    && apt-get clean -y \
    && conda env create --name hicer-nf -f environment.yml \
    && rm -rf /opt/conda/pkgs/*

ENV PATH /opt/conda/envs/hicer-nf/bin:$PATH
