FROM continuumio/miniconda
MAINTAINER Kevin Menden <kevin.menden@dzne.de>
LABEL authors="kevin.menden@dzne.de" \
    description="Docker image containing all requirements for the nf-cageseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-cageseq/bin:$PATH

# Install paraclu
RUN apt-get update && apt-get install -y make g++ unzip
RUN curl -fsSL cbrc3.cbrc.jp/~martin/paraclu/paraclu-9.zip -o /opt/paraclu-9.zip
RUN cd /opt/; unzip paraclu-9.zip; cd paraclu-9; make
ENV PATH $PATH:/opt/paraclu-9/

## Install R
#RUN apt-get update && apt-get install -y r-base
#
## Install RECLU
#COPY /assets/reclu_nfcageseq /opt/reclu_nfcageseq
#ENV PATH $PATH:/opt/reclu_nfcageseq
