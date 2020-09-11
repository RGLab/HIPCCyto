#!/bin/bash
source /app/lmod/lmod/init/profile

ml fhR/4.0.2-foss-2019b Pandoc
/fh/fast/gottardo_r/HIPCCyto/HIPCCyto/inst/scripts/update_status.R
