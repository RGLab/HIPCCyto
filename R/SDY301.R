#' @importFrom flowWorkspace load_gs gs_pop_remove gs_pop_add
apply_lymphocyte_gate_SDY301 <- function(gs, debug_dir= NULL, flowClusters = NULL){
  if(is.character(gs))
    gs <- load_gs(gs)
  if(is.character(flowClusters)){
    flowClusters = readRDS(flowClusters)
    # to check match
    sample_match <- sampleNames(gs) %in% names(flowClusters)
    if(!all(sample_match))
      stop("Precomputed flowClusters and GatingSet don't match!")
  }
  if(is.null(flowClusters))
    flowClusters <- compute_flowClusters(gs, debug_dir)

  # Just point to the cluster near (40000,25000)
  targets <- targets <- rep_len(list(c(40000,5000)), length(flowClusters))
  names(targets) <- names(flowClusters)
  gates <- create_fcEllipsoidGate(flowClusters, targets)

  message(">> Applying lymphocytes gate with flowClust by forward and side scatters (Lymphocytes)...")

  # Pull off any old/wrong Lymphocyte gate
  if(any(grepl("Lymphocytes", gs_get_pop_paths(gs))))
    gs_pop_remove(gs, "Lymphocytes")

  flowWorkspace::gs_pop_add(gs, gates, name = "Lymphocytes", parent = get_parent(gs))
  flowWorkspace::recompute(gs)
  gs
}

#' @importFrom openCyto gate_mindensity2
apply_dump_gate_SDY301 <- function(gs, plot=FALSE){
  if(is.character(gs))
    gs <- load_gs(gs)
  if(any(grepl("Lymphocytes", gs_get_pop_paths(gs))))
    gs_pop_remove(gs, "Lymphocytes")
  if(any(grepl("Dump", gs_get_pop_paths(gs))))
    gs_pop_remove(gs, "Dump")
  gates <- lapply(sampleNames(gs), function(sn){
    gate_mindensity2(gh_pop_get_data(gs[[sn]], "SCSSC"),
                     channel="Pacific Orange-A",
                     filterID="Dump",
                     pivot=TRUE,
                     gate_range=c(1.8,2.8), # Manually set
                     plot=plot
    )
  })
  names(gates) <- sampleNames(gs)
  flowWorkspace::gs_pop_add(gs, gates, name="Dump", parent = get_parent(gs))
  flowWorkspace::recompute(gs)
  gs
}
