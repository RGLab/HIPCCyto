#!/bin/bash
source /app/lmod/lmod/init/profile

ml fhR/4.0.2-foss-2019b Pandoc
Rscript --no-save --no-restore ${0%/*}/update_status.R
