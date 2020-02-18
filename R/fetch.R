#' @importFrom ImmPortR list_immport download_immport
fetch_files <- function(study, output_dir = ".") {
  # assert study
  assert_study(study)

  # list tables
  tabs <- list_immport(study)$items
  tabs <- tabs[grep("cytometry|cyto", tabs$basename, ignore.case = TRUE), ]

  # fetch tables
  output_path <- file.path(output_dir, study)
  dir.create(output_path)
  lapply(tabs$path, function(x) download_immport(x, output_path))

  # list directories
  dirs <- list_immport(file.path(study, "ResultFiles"))$items
  dirs <- dirs[grep("cytometry|cyto", dirs$basename, ignore.case = TRUE), ]

  # fetch ResultFiles directory
  output_path_resultFiles <- file.path(output_path, "ResultFiles")
  dir.create(output_path_resultFiles)
  lapply(seq_len(nrow(dirs)), function(i) {
    message(sprintf("%s: %s files (%s GB)", dirs[i, "basename"], dirs[i, "fileCount"], signif(dirs[i, "size"] / 1024 ^ 3, digits = 4)))
    download_immport(dirs[i, "path"], output_path_resultFiles)
  })

  # validate files (TODO)

  output_path
}


#' @importFrom ImmPortR query_findAllStudyAccessions
#' @importFrom assertthat assert_that
assert_study <- function(study) {
  studies <- query_findAllStudyAccessions()
  assert_that(study %in% studies, msg = paste(study, "is not a valid study."))
}
