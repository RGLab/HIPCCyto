FROM bioconductor/bioconductor_docker:RELEASE_3_10
RUN wget https://download.asperasoft.com/download/sw/cli/3.9.1/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh -P ~/
RUN chmod +x ~/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh
RUN ~/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh
RUN R -q -e 'install.packages(c("remotes", "assertthat"))'
RUN R -q -e 'BiocManager::install(c("flowCore", "ncdfFlow", "flowWorkspace", "openCyto", "flowClust", "ggcyto"))'
RUN R -q -e 'remotes::install_github(c("RGLab/ImmPortR", "RGLab/HIPCCyto"))'
