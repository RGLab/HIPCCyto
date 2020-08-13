
# HIPCCyto

<!-- badges: start -->
[![R build status](https://github.com/RGLab/HIPCCyto/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/HIPCCyto/actions)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rglab/hipccyto.svg)](https://hub.docker.com/r/rglab/hipccyto)
<!-- badges: end -->

The goal of HIPCCyto is to standardize and pre-process flow cytrometry data of [HIPC](https://www.immuneprofiling.org/) studies from [ImmPort](https://www.immport.org/).


## Local installation

You can install the development version of HIPCCyto from [GitHub](https://github.com/RGLab/HIPCCyto) with:

``` r
install.packages("remotes")
remotes::install_github("RGLab/HIPCCyto")
```

You also need to register at ImmPort and install the Aspera CLI. Please follow the instructions [here](https://rglab.github.io/ImmPortR/).


## Docker

There are many dependencies in HIPCCyto, so it takes a long time to install them all. Instead, you can use the Docker image of HIPCCyto.

``` sh
docker pull rglab/hipccyto:latest
docker run \
    -it \
    --user rstudio \
    --volume <yourLocalDirectory>:/home/rstudio \
    --env ImmPortUsername=<yourImmPortUsername> \
    --env ImmPortPassword=<yourImmPortPassword> \
    rglab/hipccyto:latest \
    R
```

Replace `<yourLocalDirectory>` with the local directory path where you'd like to store the gating sets and `<yourImmPortUsername>` and `<yourImmPortPassword>` with your Immport credential. For more information on using docker containers, please read [this documentation](https://github.com/Bioconductor/bioconductor_docker/blob/master/README.md#using-the-containers) by Bioconductor.


## Usage

### Retrieve FCS files from ImmPort

``` r
library(HIPCCyto)
file_dir <- fetch_files(study = "SDY820", output_dir = "~/")
```

### Process study

``` r
gsl <- process_study(study = "SDY820", input_dir = file_dir)
```

`process_study` function standradizes and pre-processes FCS files into gating sets.

1. Standardize marker names using ImmPort's `fcs_header_marker` table.
1. Merge metadata and batch information if available.
1. Compensate using the embeded spillover matrix in FCS files.
1. Transform the fluorescence channels using `estimateLogicle` function from `flowCore` package.
1. Pre-gate up to "Lymphocytes" using various automatic gating methods.


### Save gating sets

``` r
save_gating_sets(gsl = gsl, output_dir = file_dir)
```
