FROM bioconductor/bioconductor_docker:latest
RUN wget https://download.asperasoft.com/download/sw/cli/3.9.1/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh -P ~/
RUN chmod +x ~/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh
RUN ~/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("RGLab/HIPCCyto")'
