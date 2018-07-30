FROM continuumio/miniconda
MAINTAINER Kevin Menden <kevin.menden@dzne.de>
LABEL authors="kevin.menden@dzne.de" \
    description="Docker image containing all requirements for the nf-cageseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/seqwell/bin:$PATH



