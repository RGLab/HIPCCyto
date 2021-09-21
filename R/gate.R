# gate functions ---------------------------------------------------------------
#' @importFrom openCyto register_plugins
apply_quadrant_gate <- function(gs, study) {
  catf(">> Applying quadrant gate...")
  register_plugins(fun = .quadrantGate, methodName = "quadrantGate")

  gating_args <- "quadrant = 1"
  toRemove <- DATA[[study]]$toRemove
  if (!is.null(toRemove)) {
    catf(sprintf(">> Removing %s%% of FSC-A and SSC-A", signif(toRemove, 3) * 100))
    gating_args <- sprintf("%s, toRemove = %s", gating_args, toRemove)
  }

  gs_add_gating_method(
    gs = gs,
    alias = "pos",
    pop = "+",
    parent = get_parent(gs),
    dims = "FSC-A,SSC-A",
    gating_method = "quadrantGate",
    gating_args = gating_args
  )
}

apply_singlet_gate <- function(gs, channel) {
  alias <- sprintf("SC%s", channel)
  A <- sprintf("%s-A", channel)
  H <- sprintf("%s-H", channel)
  if (H %in% colnames(gs)) {
    catf(sprintf(">> Applying singlet gate by scatter channel (%s)...", alias))
    gs_add_gating_method(
      gs = gs,
      alias = alias,
      pop = "+",
      parent = get_parent(gs),
      dims = sprintf("%s,%s", A, H),
      gating_method = "singletGate",
      gating_args = "prediction_level = 0.99, wider_gate = TRUE"
    )
  }
}

#' @importFrom flowWorkspace sampleNames recompute gs_pop_get_parent gs_pop_get_stats gh_pop_set_gate
#' @importFrom flowCore exprs
#' @importFrom stats density
apply_nondebris_gate <- function(gs, study) {
  catf(">> Applying non-debris gate by forward scatter (Nondebris)...")

  gates <- lapply(sampleNames(gs), function(x) {
    fsc <- exprs(gh_pop_get_data(gs[[x]], get_parent(gs)))[, "FSC-A"]
    den <- density(fsc)
    minima <- ggpmisc:::find_peaks(-den$y)
    grad <- diff(den$y) / diff(den$x)
    deriv_minima <- ggpmisc:::find_peaks(-grad)
    dens_lim <- ifelse(any(minima), den$x[minima][1], -Inf)
    deriv_lim <- ifelse(any(deriv_minima), den$x[deriv_minima][1], -Inf)
    lim <- ifelse((is.finite(deriv_lim) && grad[deriv_minima][[1]] >= 0), deriv_lim, dens_lim)
    rectangleGate("FSC-A" = c(lim, Inf))
  })
  names(gates) <- sampleNames(gs)

  gs_pop_add(
    gs = gs,
    gate = gates,
    name = "Nondebris",
    parent = get_parent(gs)
  )
  recompute(gs, get_parent(gs))

  # Quick check for gates clipping large numbers of actual cells (usually because no debris peak)
  parent <- gs_pop_get_parent(gs, "Nondebris")
  ratios <- gs_pop_get_stats(gs, "Nondebris")$count / gs_pop_get_stats(gs, parent)$count
  names(ratios) <- sampleNames(gs)
  ratio_cutoff <- DATA[[study]]$Nondebris_ratio_cutoff
  if (is.null(ratio_cutoff)) {
    ratio_cutoff <- 0.0
  }
  bad_nondebris <- names(ratios[ratios < ratio_cutoff])
  if (length(bad_nondebris) > 0) {
    dummy_gate <- rectangleGate("FSC-A" = c(-Inf, Inf))
    # Give those samples fully-permissive gates
    for (name in bad_nondebris) {
      gh_pop_set_gate(gs[[name]], "Nondebris", dummy_gate)
    }
    recompute(gs, parent)
  }
}

apply_live_gate <- function(gs, study) {
  live <- get_live_marker(gs)
  if (!is.null(live)) {
    gating_method <- ifelse(is.null(DATA[[study]]$live_method), "mindensity", DATA[[study]]$live_method)
    gating_args <- ifelse(is.null(DATA[[study]]$live_args), NA, DATA[[study]]$live_args)
    collapseDataForGating <- !is.null(pData(gs)$batch)
    groupBy <- ifelse(collapseDataForGating, "batch", NA)

    catf(sprintf(">> Applying live/dead gate with %s (%s) by %s (Live)...", gating_method, gating_args, live))

    if (collapseDataForGating) {
      catf(sprintf(">> Collapsing data for gating by %s...", groupBy))
    }

    gs_add_gating_method(
      gs = gs,
      alias = "Live",
      pop = "-",
      parent = get_parent(gs),
      dims = live,
      gating_method = gating_method,
      gating_args = gating_args,
      groupBy = groupBy,
      collapseDataForGating = collapseDataForGating
    )
  }
}

#' @importFrom flowWorkspace gs_pop_add
apply_lymphocyte_gate <- function(gs, study, debug_dir = NULL) {
  flowClusters <- compute_flowClusters(gs, debug_dir)
  targets <- compute_targets(gs, flowClusters, study)
  gates <- create_fcEllipsoidGate(flowClusters, targets)

  catf(">> Applying lymphocytes gate with flowClust by forward and side scatters (Lymphocytes)...")
  gs_pop_add(
    gs = gs,
    gate = gates,
    name = "Lymphocytes",
    parent = get_parent(gs)
  )
  recompute(gs, get_parent(gs))
}


# quadrant gate helper functions -----------------------------------------------
#' @importFrom flowCore exprs rectangleGate
.quadrantGate <- function(fr, pp_res, channels = NA, filterId = "", toRemove = 0, quadrant, ...) {
  lim_func <- switch(quadrant,
    "1" = c(max, max),
    "2" = c(min, max),
    "3" = c(min, min),
    "4" = c(max, min)
  )
  ch1_lim <- lim_func[[1]](exprs(fr)[, channels[1]], na.rm = TRUE) * (1 - toRemove)
  ch2_lim <- lim_func[[2]](exprs(fr)[, channels[1]], na.rm = TRUE) * (1 - toRemove)

  range <- switch(quadrant,
    "1" = c(0, ch1_lim, 0, ch2_lim),
    "2" = c(-ch1_lim, 0, 0, ch2_lim),
    "3" = c(-ch1_lim, 0, -ch2_lim, 0),
    "4" = c(0, ch1_lim, -ch2_lim, 0)
  )
  mat <- matrix(range, ncol = 2, dimnames = list(c("min", "max"), c(channels[1], channels[2])))
  rectangleGate(filterId = filterId, .gate = mat)
}


# live gate helper functions ---------------------------------------------------
get_live_marker <- function(gs) {
  live <- grep("^(L|l)ive|LD|(V|v)iability|L/D$", markernames(gs), value = TRUE)

  if (length(live) == 0) {
    catf(">> There is no viability dye channel in this gating set...")
    return(NULL)
  } else if (length(live) > 1) {
    catf(">> There are more than one viability dye channel in this gating set...")
    catf(paste(live, collapse = ", "))
  }

  live
}


# lymphocyte gate helper functions ---------------------------------------------
#' @importFrom flowClust flowClust
flowclust <- function(x) {
  fcl <- flowClust(
    x = x,
    K = 1:5,
    criterion = "ICL",
    trans = 0,
    min.count = -1,
    max.count = -1
  )
  fc <- fcl@.Data[[fcl@index]]
  fc@z <- matrix()
  fc@u <- matrix()
  fc
}

#' @importFrom flowWorkspace gs_pop_get_data sampleNames
#' @importFrom flowCore exprs
#' @importFrom slurmR slurm_available Slurm_lapply opts_slurmR
compute_flowClusters <- function(gs, debug_dir = NULL) {
  catf(">> Computing for the optimal number of clusters (K) for each sample...")
  cs <- gs_pop_get_data(gs, get_parent(gs))

  if (slurm_available()) {
    catf(">> Submitting flowClust jobs to slurm...")
    ex <- lapply(sampleNames(cs), function(x) exprs(cs[[x, returnType = "cytoframe"]])[, c("FSC-A", "SSC-A")])
    names(ex) <- sampleNames(cs)
    if (is.null(debug_dir)) {
      tmp_path <- opts_slurmR$get_tmp_path()
    } else {
      tmp_path <- debug_dir
    }
    flowClusters <- Slurm_lapply(
      ex, flowclust,
      njobs = length(ex), mc.cores = 1L, tmp_path = tmp_path,
      sbatch_opt = list("constraint" = "gizmok")
    )
  } else {
    flowClusters <- mclapply(sampleNames(cs), function(x) {
      ex <- exprs(cs[[x, returnType = "cytoframe"]])[, c("FSC-A", "SSC-A")]
      flowclust(ex)
    }, mc.cores = detect_cores())
    names(flowClusters) <- sampleNames(cs)
  }

  save_debug(flowClusters, "compute_flowClusters", debug_dir)

  flowClusters
}

#' @importFrom stats dist
#' @importFrom utils tail
select_cluster <- function(fitted_means, target) {
  target_dist <- as.matrix(dist(rbind(fitted_means, target)))
  target_dist <- tail(target_dist, n = 1)[seq_len(nrow(fitted_means))]
  which.min(target_dist)
}

#' @importFrom flowClust getEstimates
#' @importFrom withr with_options
find_target <- function(flowClusters) {
  catf(">> Computing the target location of the lymphocyte clusters...")
  mus <- lapply(flowClusters, function(x) {
    est <- getEstimates(x)
    est$locations[which.max(est$proportions), ]
  })
  mus <- do.call(rbind, mus)
  colnames(mus) <- c("FSC", "SSC")
  fcl_mus <- with_options(
    list(mc.cores = 1L),
    flowClust(mus, K = 1:5, criterion = "ICL", trans = 0, min.count = -1, max.count = -1)
  )
  k_mus <- ifelse(length(fcl_mus@index) == 0, 1, fcl_mus@index)
  est_mus <- getEstimates(fcl_mus@.Data[[k_mus]])
  target <- est_mus$locations[which.max(est_mus$proportions), ]

  print(est_mus)
  catf(sprintf(">> Selecting FSC-A = %s and SSC-A = %s as target location...", target[1], target[2]))

  targets <- rep_len(list(target), length(flowClusters))
  names(targets) <- names(flowClusters)

  targets
}

compute_targets <- function(gs, flowClusters, study) {
  target <- DATA[[study]]$target
  if (!is.null(target)) {
    catf(sprintf(">> Using the predetermined target location (FSC-A = %s and SSC-A = %s)...", target[1], target[2]))
    targets <- rep_len(list(target), length(flowClusters))
    names(targets) <- names(flowClusters)
    return(targets)
  }

  batch <- pData(gs)$batch
  if (is.null(batch)) {
    targets <- find_target(flowClusters)
  } else {
    pd <- pData(gs)[, c("name", "batch")]
    pd <- split(pd, factor(pd$batch, exclude = NULL))

    list_targets <- lapply(pd, function(x) {
      catf(unique(x$batch))
      find_target(flowClusters[x$name])
    })
    targets <- unlist(unname(list_targets), recursive = FALSE)
  }

  targets
}

create_fcEllipsoidGate <- function(flowClusters, targets) {
  quantile <- 0.9
  trans <- 0
  prior <- list(NA)

  gates <- lapply(names(flowClusters), function(sample) {
    tmix_results <- flowClusters[[sample]]
    target <- targets[[sample]]

    fitted_means <- getEstimates(tmix_results)$locations
    cluster_selected <- select_cluster(fitted_means, target)

    posteriors <- list(
      mu = tmix_results@mu,
      lambda = tmix_results@lambda,
      sigma = tmix_results@sigma,
      nu = tmix_results@nu
    )

    flowClust_gate <- openCyto:::.getEllipseGate(
      filter = tmix_results,
      include = cluster_selected,
      quantile = quantile,
      trans = trans
    )

    openCyto:::fcEllipsoidGate(flowClust_gate, prior, posteriors)
  })
  names(gates) <- names(flowClusters)

  gates
}
