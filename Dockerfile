FROM continuumio/miniconda
MAINTAINER Kevin Menden <kevin.menden@dzne.de>
LABEL authors="kevin.menden@dzne.de" \
    description="Docker image containing all requirements for the nf-cageseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/seqwell/bin:$PATH

# Install DropSeq Tools
RUN apt-get update && apt-get install -y unzip
RUN curl -fsSL mccarrolllab.com/download/1276/ -o /opt/dropseq_tools.zip
RUN cd /opt/; unzip dropseq_tools.zip;
ENV PATH $PATH:/opt/Drop-seq_tools-1.13/
