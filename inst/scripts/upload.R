#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

.libPaths(c(
  "/fh/fast/gottardo_r/HIPCCyto/R",
  "/app/software/fhR/4.0.2-foss-2019b",
  "/app/software/R/4.0.2-foss-2019b/lib/R/library"
))

study <- args[1]
data_dir <- "/fh/fast/gottardo_r/HIPCCyto/data"
if (!dir.exists(file.path(data_dir, study))) {
  stop("No study found...")
}
sdy_dir <- file.path(data_dir, study, "GatingSets")
versions <- dir(sdy_dir, pattern = "^v\\d+.\\d+.\\d+$")
if (length(versions) == 0) {
  stop("No gating sets found...")
}
latest <- versions[order(semver::parse_version(gsub("^v", "", versions)), decreasing = TRUE)][1]
ver_dir <- file.path(sdy_dir, latest)

message(latest)
files <- dir(ver_dir, pattern = ".html", full.names = TRUE, recursive = TRUE)
system.time(res <- lapply(files, function(file) {
  object <- gsub(data_dir, "", file)
  message(object)
  try(aws.s3::put_object(
    file = file, object = object, bucket = "hipccyto", multipart = TRUE,
    region = "us-east-2", headers = list(`Content-Type` = "text/html")
  ))
}))
