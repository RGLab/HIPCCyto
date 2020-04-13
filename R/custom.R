#' @importFrom flowWorkspace load_gs gs_pop_remove gs_pop_add
apply_lymphocyte_gate_custom <- function(gs, debug_dir= NULL, flowClusters = NULL, nclust = NULL, target = c(40000,5000)) {
  if (is.character(gs))
    gs <- load_gs(gs)
  if (is.character(flowClusters)) {
    flowClusters = readRDS(flowClusters)
    # to check match
    sample_match <- sampleNames(gs) %in% names(flowClusters)
    if (!all(sample_match))
      stop("Precomputed flowClusters and GatingSet don't match!")
  }
  if (is.null(flowClusters))
    if (is.null(nclust))
      flowClusters <- compute_flowClusters(gs, debug_dir)
    else
      flowClusters <- compute_flowClusters_custom(gs, debug_dir, nclust)

  # Should really keep these in a separate tabulated data structure
  # SDY301: c(40000,5000)
  # SDY144: c(60000,5000)
  targets <- rep_len(list(target), length(flowClusters))
  names(targets) <- names(flowClusters)
  gates <- create_fcEllipsoidGate(flowClusters, targets)

  message(">> Applying lymphocytes gate with flowClust by forward and side scatters (Lymphocytes)...")

  # Pull off any old/wrong Lymphocyte gate
  if (any(grepl("Lymphocytes", gs_get_pop_paths(gs))))
    gs_pop_remove(gs, "Lymphocytes")

  gs_pop_add(gs, gates, name = "Lymphocytes", parent = get_parent(gs), recompute = TRUE)

  gs
}

#' @importFrom openCyto gate_mindensity2
apply_dump_gate_custom <- function(gs, plot=FALSE, gate_range = c(1.8,2.8)) {
  if (is.character(gs))
    gs <- load_gs(gs)
  if (any(grepl("Lymphocytes", gs_get_pop_paths(gs))))
    gs_pop_remove(gs, "Lymphocytes")
  if (any(grepl("Dump", gs_get_pop_paths(gs))))
    gs_pop_remove(gs, "Dump")
  gates <- lapply(sampleNames(gs), function(sn) {
    gate_mindensity2(
      fr = gh_pop_get_data(gs[[sn]], "SCSSC"),
      channel = "Pacific Orange-A",
      filterID = "Dump",
      pivot = TRUE,
      # SDY301: c(1.8, 2.8)
      gate_range = gate_range, # Manually set
      plot = plot
    )
  })
  names(gates) <- sampleNames(gs)
  gs_pop_add(gs, gates, name = "Dump", parent = get_parent(gs), recompute = TRUE)

  gs
}

compute_flowClusters_custom <- function(gs, debug_dir = NULL, nclust = 1) {
  message(">> Computing for the optimal number of clusters (K) for each sample...")
  nc <- gs_pop_get_data(gs, get_parent(gs))
  flowClusters <- mclapply(sampleNames(nc), function(x) {
    fc <- flowClust(
      x = nc[[x]]@exprs[, c("FSC-A", "SSC-A")],
      K = nclust,
      trans = 0,
      min.count = -1,
      max.count = -1
    )
    fc@z <- matrix()
    fc@u <- matrix()
    fc
  }, mc.cores = detect_cores())
  names(flowClusters) <- sampleNames(nc)

  save_debug(flowClusters, "compute_flowClusters", debug_dir)

  flowClusters
}
