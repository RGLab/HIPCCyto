#!/usr/bin/env Rscript

.libPaths(c(
  "/fh/fast/gottardo_r/HIPCCyto/R",
  "/app/software/fhR/4.1.1-foss-2020b",
  "/app/software/R/4.1.1-foss-2020b/lib/R/library"
))
options(bitmapType = "cairo")

# render status dash board
file <- "/fh/fast/gottardo_r/HIPCCyto/status.html"
rmarkdown::render(
  input = system.file("scripts/status.Rmd", package = "HIPCCyto"),
  output_file = file,
  intermediates_dir = tempdir(),
  params = list(immport = TRUE)
)
file_dev <- "/fh/fast/gottardo_r/HIPCCyto/status_dev.html"
rmarkdown::render(
  input = system.file("scripts/status.Rmd", package = "HIPCCyto"),
  output_file = file_dev,
  intermediates_dir = tempdir(),
  params = list(dev = TRUE, immport = TRUE)
)

# upload to s3 bucket
aws.s3::put_object(
  file = file,
  object = "status.html",
  bucket = "fh-pi-gottardo-r-eco-public/hipccyto",
  acl = "public-read",
  region = "us-west-2"
)
aws.s3::put_object(
  file = file_dev,
  object = "status_dev.html",
  bucket = "fh-pi-gottardo-r-eco-public/hipccyto",
  acl = "public-read",
  region = "us-west-2"
)
