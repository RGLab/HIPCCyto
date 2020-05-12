#' @importFrom parallel detectCores
detect_cores <- function() {
  getOption("mc.cores", detectCores())
}

catf <- function(..., file = "", sep = " ", fill = TRUE, labels = NULL,
                 append = FALSE) {
  cat(..., file = file, sep = sep, fill = fill, labels = labels, append = append)
}
