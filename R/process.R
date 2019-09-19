process_study <- function(study, input_dir) {
  # assert
  script <- paste0("process_", study)
  assert_that(exists(script))

  # run processing script
  get(script)(input_dir)
}
