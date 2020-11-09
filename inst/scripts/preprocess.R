#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

.libPaths(c(
  "/fh/fast/gottardo_r/HIPCCyto/R",
  "/app/software/fhR/4.0.2-foss-2019b",
  "/app/software/R/4.0.2-foss-2019b/lib/R/library"
))
options(bitmapType = "cairo")

print(args)
print(.libPaths())

library(HIPCCyto)

t <- Sys.time()
study <- args[1]
options(mc.cores = as.integer(args[2]))
ver <- paste0("v", packageVersion("HIPCCyto"))

input_dir <- sprintf("/fh/fast/gottardo_r/HIPCCyto/data/%s/ResultFiles/Flow_cytometry_result", study)
debug_dir <- sprintf("/fh/scratch/delete10/gottardo_r/HIPCCyto/%s", study)
output_dir <- sprintf("/fh/fast/gottardo_r/HIPCCyto/data/%s/GatingSets/%s", study, ver, recursive = TRUE)

if (!file.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)
if (!file.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gsl <- process_study(study, input_dir, debug_dir)

save_gating_sets(gsl, output_dir)

print(Sys.time() - t)
