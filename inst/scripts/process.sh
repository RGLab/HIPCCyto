#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=16384
#SBATCH --constraint=gizmok

ml fhR/4.0.2-foss-2019b Pandoc
Rscript --no-save --no-restore $wd/process.R $sdy 4
