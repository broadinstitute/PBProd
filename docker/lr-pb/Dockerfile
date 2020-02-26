FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-pb/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt update
RUN apt install -y python2.7 curl
RUN curl https://sdk.cloud.google.com | bash

# copy other resources
COPY detect_run_info.py /usr/local/bin/
COPY merge_ccs_reports.py /usr/local/bin/
COPY gatk.jar /usr/local/bin/