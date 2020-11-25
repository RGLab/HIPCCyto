#' @importFrom parallel detectCores
detect_cores <- function() {
  getOption("mc.cores", detectCores())
}

catf <- function(..., file = "", sep = " ", fill = TRUE, labels = NULL,
                 append = FALSE) {
  time <- sprintf("[%s]", as.character(Sys.time()))
  cat(time, ..., file = file, sep = sep, fill = fill, labels = labels, append = append)
}

#' @importFrom ImmPortR query_datarelversion
get_dr <- function() {
  paste0("DR", query_datarelversion())
}

# fluorescence channel names
#' @importFrom flowWorkspace gh_pop_get_data
colnames2 <- function(gs) {
  cf <- gh_pop_get_data(gs[[1]])
  spillover <- get_spillover(cf)
  colnames(spillover)
}

#' @importFrom flowWorkspace gs_get_pop_paths
get_parent <- function(gs) {
  rev(gs_get_pop_paths(gs, path = 1))[1]
}

get_markers <- function(study, modify = TRUE) {
  headers <- ImmPortR:::query(sprintf("fcs_header_marker/%s", study))

  pns <- unique(headers[, c("pnsPreferred", "pnsReported")])
  pns <- pns[!pns$pnsReported %in% headers$pnnReported, ]

  markers <- pns$pnsPreferred
  names(markers) <- pns$pnsReported

  # use reported marker names when standardized names are not available
  markers[is.na(markers)] <- names(markers)[is.na(markers)]

  if (isTRUE(modify)) {
    toModify <- names(markers)[names(markers) %in% names(MARKERS)]
    markers[toModify] <- MARKERS[toModify]
  }
  markers
}

get_spillover <- function(x) {
  spills <- spillover(x)
  spills[!sapply(spills, is.null)][[1]]
}

get_marker_channel <- function(gs) {
  markers <- markernames(gs)
  gsub("/", "_", paste0(markers, "_", names(markers)))
}

get_nodes <- function(gs) {
  gs_get_pop_paths(gs, path = 1)[-1]
}

#' @importFrom utils packageVersion
get_version <- function() {
  paste0("v", packageVersion("HIPCCyto"))
}

#' @importFrom mime guess_type
#' @importFrom xfun base64_encode
encode_img <- function(file) {
  paste0("data:", guess_type(file), ";base64,", base64_encode(file))
}
