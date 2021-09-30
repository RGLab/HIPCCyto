#!/bin/bash
source /app/lmod/lmod/init/profile

ml fhR/4.1.1-foss-2020b Pandoc
Rscript --no-save --no-restore ${0%/*}/update_status.R
