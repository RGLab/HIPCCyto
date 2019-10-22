process_study <- function(study, input_dir, output_dir = NULL) {
  # assert
  script <- paste0("process_", study)
  assert_that(exists(script))

  # run processing script
  gs <- get(script)(input_dir)

  # save gating set
  if (is.null(output_dir)) {
    output_dir <- file.path(input_dir, "gs")
  }
  flowWorkspace::save_gs(gs, output_dir)
}
