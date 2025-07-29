#Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3
USER root

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
## Create a new environment and install packages
RUN micromamba install -y -n base -c conda-forge -c bioconda -c dnachun \
    conda-forge::r-tidyverse \
    conda-forge::datatable \
    conda-forge::r-optparse \
    conda-forge::r-enrichr  \
    conda-forge::r-gprofiler2 \
    dnachun::r-twosamplemr \
    conda-forge::r-r.utils \
    conda-forge::r-arrow 


COPY mendelian_randomization.R . 

## Default command to run when the container starts
CMD ["bash"]
