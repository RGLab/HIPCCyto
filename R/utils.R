#' @importFrom parallel detectCores
detect_cores <- function() {
  getOption("mc.cores", detectCores())
}

catf <- function(..., file = "", sep = " ", fill = TRUE, labels = NULL,
                 append = FALSE) {
  cat(..., file = file, sep = sep, fill = fill, labels = labels, append = append)
}

#' @importFrom ImmPortR query_datarelversion
get_dr <- function() {
  paste0("DR", query_datarelversion())
}

#' @importFrom flowCore parameters
#' @importFrom flowWorkspace gh_pop_get_data
colnames2 <- function(gs) {
  # can't use this for now
  # grep("SC-|Time", colnames(gs), invert = TRUE, value = TRUE)

  channels <- parameters(gh_pop_get_data(gs[[1]]))@data$name
  markers <- parameters(gh_pop_get_data(gs[[1]]))@data$desc

  unname(channels[!is.na(markers)])
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

  if (isTRUE(modify)) {
    toModify <- names(markers)[names(markers) %in% names(MARKERS)]
    markers[toModify] <- markers[MARKERS[toModify]]
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

get_version <- function() {
  paste0("v", packageVersion("HIPCCyto"))
}
