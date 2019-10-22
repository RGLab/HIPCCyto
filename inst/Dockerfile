FROM bioconductor/release_cytometry2
RUN wget https://download.asperasoft.com/download/sw/cli/3.9.1/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh -P ~/
RUN chmod +x ~/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh
RUN ~/ibm-aspera-cli-3.9.1.1401.be67d47-linux-64-release.sh
RUN R -q -e 'install.packages(c("devtools", "assertthat"))'
RUN R -q -e 'devtools::install_github("RGLab/ImmPortR")'
ADD . /home/bioc/irc-network
RUN R -q -e 'devtools::install("/home/bioc/irc-network")'
