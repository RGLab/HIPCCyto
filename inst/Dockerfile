FROM bioconductor/bioconductor_docker:latest
RUN wget https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/08q6g/0/ibm-aspera-cli-3.9.6.1467.159c5b1-linux-64-release.sh -P ~/
RUN chmod +x ~/ibm-aspera-cli-3.9.6.1467.159c5b1-linux-64-release.sh
RUN ~/ibm-aspera-cli-3.9.6.1467.159c5b1-linux-64-release.sh
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("RGLab/HIPCCyto")'
# To mount Fred Hutch file systems
# https://sciwiki.fredhutch.org/compdemos/Singularity
RUN mkdir -p /fh/fast && mkdir /fh/scratch
