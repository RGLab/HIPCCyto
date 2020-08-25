#!/bin/bash
o=/fh/fast/gottardo_r/HIPCCyto/logs/$1_$(date +"%Y-%m-%d-%H:%M:%S")
user=$(whoami)@fredhutch.org
sbatch -o $o -J $1 --export=sdy=$1 --mail-user=$user process.sh
