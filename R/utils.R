#' @importFrom parallel detectCores
detect_cores <- function() {
  getOption("mc.cores", detectCores())
}
