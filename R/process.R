#' @importFrom parallel mclapply
#' @export
process_study <- function(study, input_dir, debug_dir = NULL) {
  # summarize files
  files <- summarize_study(study, input_dir, debug_dir = debug_dir)
  files_by_panel <- split(files, files$panel)

  # create gating set for each panel
  lapply(files_by_panel, process_panel, debug_dir = debug_dir)
}


#' @importFrom ImmPortR query_filePath
#' @importFrom flowCore read.FCSheader
#' @importFrom gtools mixedsort
summarize_study <- function(study, input_dir, remove_dups = TRUE, standardize_markernames = TRUE, debug_dir = NULL) {
  files <- query_filePath(study)
  files <- files[files$fileDetail == "Flow cytometry result", ]
  files <- files[grepl(".fcs$", files$fileName), ]
  files$filePathId <- NULL
  files$sourceAccession <- NULL
  files$filePath <- file.path(input_dir, files$fileName)
  files <- unique(files)
  if (isTRUE(remove_dups)) {
    files <- files[!duplicated(files$fileInfoId), ]
  }
  rownames(files) <- files$fileName

  # check if files exist. If not, throw warning and fetch them
  fcs_exist <- file.exists(files$filePath)
  if (any(!fcs_exist)) {
    stop("These fcs files do not exist at ", input_dir, "\n", paste(files[!fcs_exist, "fileName"], collapse = "\n"))
  }

  # summarize files
  catf(sprintf(">> There are %s fcs files in %s...", nrow(files), study))

  # read headers and summarize panels
  map <- DATA[[study]]$map
  catf(">> Reading in fcs headers...")
  headers <- mclapply(files$filePath, function(file) read.FCSheader(file, channel_alias = map), mc.cores = detect_cores())

  files$tot <- sapply(headers, function(x) x[[1]]["$TOT"])
  files$par <- sapply(headers, function(x) x[[1]]["$PAR"])

  files$src <- sapply(headers, function(x) x[[1]]["$SRC"])
  files$date <- sapply(headers, function(x) x[[1]]["$DATE"])
  files$btim <- sapply(headers, function(x) x[[1]]["$BTIM"])
  files$etim <- sapply(headers, function(x) x[[1]]["$ETIM"])
  files$cyt <- sapply(headers, function(x) x[[1]]["$CYT"])

  files$creator <- sapply(headers, function(x) x[[1]]["CREATOR"])
  files$tubeName <- sapply(headers, function(x) x[[1]]["TUBE NAME"])
  files$experimentName <- sapply(headers, function(x) x[[1]]["EXPERIMENT NAME"])
  files$settings <- sapply(headers, function(x) x[[1]]["SETTINGS"])
  files$cytnum <- sapply(headers, function(x) x[[1]]["CYTNUM"])
  files$exportUserName <- sapply(headers, function(x) x[[1]]["EXPORT USER NAME"])
  files$exportTime <- sapply(headers, function(x) x[[1]]["EXPORT TIME"])
  files$asf <- sapply(headers, function(x) x[[1]]["FSC ASF"])
  files$plateName <- sapply(headers, function(x) x[[1]]["PLATE NAME"])
  files$plateId <- sapply(headers, function(x) x[[1]]["PLATE ID"])
  files$wellId <- sapply(headers, function(x) x[[1]]["WELL ID"])

  files$cstSetupStatus <- sapply(headers, function(x) x[[1]]["CST SETUP STATUS"])
  files$cstBeadsLotId <- sapply(headers, function(x) x[[1]]["CST BEADS LOT ID"])
  files$cstSetupDate <- sapply(headers, function(x) x[[1]]["CST SETUP DATE"])
  files$cstBaselineDate <- sapply(headers, function(x) x[[1]]["CST BASELINE DATE"])
  files$cytometerConfigName <- sapply(headers, function(x) x[[1]]["CYTOMETER CONFIG NAME"])
  files$cytometerConfigCreateDate <- sapply(headers, function(x) x[[1]]["CYTOMETER CONFIG CREATE DATE"])

  markers <- get_markers(study)

  panels <- sapply(headers, function(x) {
    header <- x[[1]]
    par <- as.integer(header["$PAR"])
    PNN <- unname(header[paste0("$P", seq_len(par), "N")])
    PNS <- unname(header[paste0("$P", seq_len(par), "S")])

    # standardize channel names
    if (!is.null(map)) {
      for (i in seq_len(nrow(map))) {
        PNN <- gsub(map$channels[i], map$alias[i], PNN)
      }
    }

    # standardize marker names
    if (isTRUE(standardize_markernames)) {
      marker_exist <- !is.na(PNS) & PNS %in% names(markers)
      PNS[marker_exist] <- markers[PNS[marker_exist]]
    }

    PNN <- PNN[!is.na(PNS)]
    PNS <- PNS[!is.na(PNS)]

    if (length(PNS) == 0) {
      ""
    } else {
      paste(mixedsort(paste0(PNS, " (", PNN, ")")), collapse = "; ")
    }
  })

  files$panel <- panels
  ps <- unique(panels)

  # how many gating sets will be created
  # by panels, sample type, measurement technique, experiment accession
  catf(sprintf(">> There are %s panel(s)...", length(ps)))
  panels_clean <- sapply(strsplit(ps, "; "), function(x) {
    if (length(x) == 0) {
      ""
    } else {
      paste(gsub(" \\(.+\\)$", "", x), collapse = " | ")
    }
  })
  catf(paste(mixedsort(panels_clean), collapse = "\n"))

  save_debug(files, "summarize_study", debug_dir)

  files
}


process_panel <- function(files, debug_dir = NULL) {
  study <- unique(files$studyAccession)
  panel <- unique(files$panel)
  stopifnot(length(study) == 1, length(panel) == 1)

  catf(paste(rep("=", times = 80), collapse = ""))
  catf(sprintf(">> Processing %s fcs files for this panel...", nrow(files)))
  catf(paste(strsplit(panel, split = "; ")[[1]], collapse = "\n"))

  # load files
  print(system.time(nc <- create_nc(files$filePath, study, debug_dir)))

  # merge metadata
  print(system.time(nc <- merge_metadata(nc, files, study, debug_dir)))

  # merge batch information
  print(system.time(nc <- merge_batch(nc, study, debug_dir)))

  # create a gating set
  print(system.time(gs <- create_gs(nc, study, debug_dir)))

  # pre-process
  print(system.time(gs <- standardize_markernames(gs, study, debug_dir)))
  print(system.time(gs <- compensate_gs(gs, study, debug_dir)))
  print(system.time(gs <- transform_gs(gs, study, debug_dir)))

  # gate
  gate_gs(gs, study, debug_dir)

  gs
}


#' @importFrom ncdfFlow read.ncdfFlowSet
create_nc <- function(filePath, study, debug_dir = NULL) {
  catf(">> Reading files and creating a flow set...")
  map <- DATA[[study]]$map

  nc <- suppressMessages(read.ncdfFlowSet(
    filePath,
    channel_alias = map,
    mc.cores = detect_cores()
  ))

  save_debug(nc, "create_nc", debug_dir)

  nc
}

#' @importFrom ncdfFlow phenoData phenoData<-
merge_metadata <- function(nc, files, study, debug_dir = NULL) {
  catf(">> Merging metedata...")
  phenoData(nc)$study_accession <- files[phenoData(nc)$name, ]$studyAccession
  phenoData(nc)$participant_id <- files[phenoData(nc)$name, ]$subjectAccession
  phenoData(nc)$age_reported <- files[phenoData(nc)$name, ]$ageEvent
  phenoData(nc)$gender <- files[phenoData(nc)$name, ]$gender
  phenoData(nc)$race <- files[phenoData(nc)$name, ]$race
  phenoData(nc)$study_time_collected <- files[phenoData(nc)$name, ]$studyTimeCollected
  phenoData(nc)$study_time_collected_unit <- files[phenoData(nc)$name, ]$studyTimeCollectedUnit
  phenoData(nc)$file_info_name <- files[phenoData(nc)$name, ]$fileName
  phenoData(nc)$description <- files[phenoData(nc)$name, ]$fileDetail
  phenoData(nc)$type <- files[phenoData(nc)$name, ]$biosampleType
  phenoData(nc)$subtype <- files[phenoData(nc)$name, ]$biosampleSubtype
  phenoData(nc)$cohort <- files[phenoData(nc)$name, ]$armName

  save_debug(nc, "merge_metadata", debug_dir)

  nc
}

#' @importFrom flowCore fsApply description
merge_batch <- function(nc, study, debug_dir = NULL) {
  catf(">> Merging batch column...")
  keyword <- DATA[[study]]$batch

  if (!is.null(keyword)) {
    phenoData(nc)$batch <- unlist(
      fsApply(
        x = nc,
        FUN = function(x) {
          val <- description(x)[keyword][[1]]
          if (is.null(val)) val <- NA
          val
        },
        simplify = FALSE
      )[phenoData(nc)$name],
      use.names = FALSE
    )
  }

  save_debug(nc, "merge_batch", debug_dir)

  nc
}

#' @importFrom flowWorkspace GatingSet
create_gs <- function(nc, study, debug_dir = NULL) {
  catf(">> Creating a gating set...")
  gs <- GatingSet(nc)

  save_debug(gs, "create_gs", debug_dir)

  gs
}

#' @importFrom flowWorkspace markernames markernames<-
standardize_markernames <- function(gs, study, debug_dir = NULL) {
  catf(">> Standardizing marker names...")

  # get current marker names and name them with channel names
  names_gs <- markernames(gs)
  if (is(names_gs, "list")) {
    names_gs <- names_gs[[1]]
  }
  names(names_gs) <- colnames2(gs)

  # retrieve standard marker names from ImmPort
  markers <- get_markers(study)
  marker_exist <- names_gs %in% names(markers)
  standards <- markers[names_gs[marker_exist]]
  names_gs[marker_exist] <- standards

  # assign the fixed marker names
  standards <- standards[standards != names(standards)]
  catf(paste(paste(names(standards), standards, sep = " -> "), collapse = "\n"))
  markernames(gs) <- names_gs

  save_debug(gs, "standardize_markernames", debug_dir)

  gs
}

#' @importFrom flowWorkspace gs_pop_get_data compensate
#' @importFrom flowCore spillover
compensate_gs <- function(gs, study, debug_dir = NULL) {
  catf(">> Applying compensation...")
  nc <- gs_pop_get_data(gs)
  cols <- colnames(gs)
  comp <- fsApply(nc, function(x) {
    spills <- spillover(x)
    spill <- spills[!sapply(spills, is.null)][[1]] # pick the first non-empty matrix
    keep <- colnames(spill) %in% cols
    spill[keep, keep] # remove extra channels
  }, simplify = FALSE)

  gs <- compensate(gs, comp)

  save_debug(gs, "compensate_gs", debug_dir)

  gs
}

# transform fluoresence channels with biexponential transformation
#' @importFrom flowWorkspace colnames transform
#' @importFrom flowCore estimateLogicle
transform_gs <- function(gs, study, debug_dir = NULL) {
  catf(">> Applying transformation...")
  channels <- colnames2(gs)
  if (length(channels) > 0) {
    trans <- estimateLogicle(gs[[1]], channels)
    gs <- transform(gs, trans)
  }

  save_debug(gs, "transform_gs", debug_dir)

  gs
}

#' @importFrom openCyto gatingTemplate gt_gating gs_add_gating_method
gate_gs <- function(gs, study, debug_dir = NULL) {
  filePath <- sprintf("extdata/gating_template/%s.csv", study)
  file <- system.file(filePath, package = "HIPCCyto")

  if (file != "") {
    catf(sprintf(">> Applying gating template (%s)...", file))
    gt <- gatingTemplate(file)
    gt_gating(gt, gs, mc.cores = detect_cores(), parallel_type = "multicore")
  } else {
    catf(">> Gating template does not exist for this study...")
    catf(">> Applying default gating methods...")
    print(system.time(apply_quadrant_gate(gs, study)))
    print(system.time(apply_singlet_gate(gs, "FSC")))
    print(system.time(apply_singlet_gate(gs, "SSC")))
    print(system.time(apply_live_gate(gs, study)))
    print(system.time(apply_nondebris_gate(gs)))
    print(system.time(apply_lymphocyte_gate(gs, study, debug_dir)))
  }

  save_debug(gs, "gate_gs", debug_dir)

  gs
}


# HELPER FUNCTIONS -------------------------------------------------------------

#' @importFrom ncdfFlow getFileName save_ncfs
#' @importFrom flowWorkspace save_gs
save_debug <- function(obj, func, debug_dir = NULL) {
  if (!is.null(debug_dir)) {
    path <- tempfile(paste0(func, "_"), debug_dir)
    catf(sprintf(">> Storing intermediate file to %s for debugging...", path))
    if (is(obj, "ncdfFlowSet")) {
      save_ncfs(obj, path, cdf = "copy")
    } else if (is(obj, "GatingSet")) {
      save_gs(obj, path, cdf = "copy")
    } else {
      saveRDS(obj, path)
    }
  }
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

#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom slurmR slurm_available Slurm_lapply
compute_flowClusters <- function(gs, debug_dir = NULL) {
  catf(">> Computing for the optimal number of clusters (K) for each sample...")
  nc <- gs_pop_get_data(gs, get_parent(gs))

  if (slurm_available()) {
    catf(">> Submitting flowClust jobs to slurm...")
    ex <- lapply(sampleNames(nc), function(x) nc[[x]]@exprs[, c("FSC-A", "SSC-A")])
    names(ex) <- sampleNames(nc)
    flowClusters <- Slurm_lapply(ex, flowclust, njobs = length(ex), mc.cores = 1L, sbatch_opt = list("constraint" = "gizmok"))
  } else {
    flowClusters <- mclapply(sampleNames(nc), function(x) {
      ex <- nc[[x]]@exprs[, c("FSC-A", "SSC-A")]
      flowclust(ex)
    }, mc.cores = detect_cores())
    names(flowClusters) <- sampleNames(nc)
  }

  save_debug(flowClusters, "compute_flowClusters", debug_dir)

  flowClusters
}

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
  quantile =  0.9
  trans = 0
  prior = list(NA)

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

#' @importFrom flowCore rectangleGate
.quadrantGate <- function(fr, pp_res, channels = NA, filterId = "", toRemove = 0, quadrant, ...) {
  lim_func <- switch(
    quadrant,
    "1" = c(max, max),
    "2" = c(min, max),
    "3" = c(min, min),
    "4" = c(max, min)
  )
  ch1_lim <- lim_func[[1]](fr@exprs[, channels[1]], na.rm = TRUE) * (1 - toRemove)
  ch2_lim <- lim_func[[2]](fr@exprs[, channels[1]], na.rm = TRUE) * (1 - toRemove)

  range <- switch(
    quadrant,
    "1" = c(0, ch1_lim, 0, ch2_lim),
    "2" = c(-ch1_lim, 0, 0, ch2_lim),
    "3" = c(-ch1_lim, 0, -ch2_lim, 0),
    "4" = c(0, ch1_lim, -ch2_lim, 0)
  )
  mat <- matrix(range, ncol = 2, dimnames = list(c("min", "max"), c(channels[1], channels[2])))
  rectangleGate(filterId = filterId, .gate = mat)
}

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

#' @importFrom flowWorkspace recompute
apply_nondebris_gate <- function(gs) {
  catf(">> Applying non-debris gate by forward scatter (Nondebris)...")

  gates <- lapply(sampleNames(gs), function(x) {
    fsc <- exprs(gh_pop_get_data(gs[[x]], get_parent(gs)))[, "FSC-A"]
    den <- density(fsc)
    peaks <- ggpmisc:::find_peaks(-den$y)
    if (peaks == 0) {
      grad <- diff(den$y) / diff(den$x)
      peaks <- ggpmisc:::find_peaks(-grad)
    }
    lim <- den$x[peaks][1]
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

get_live_marker <- function(gs) {
  live <- grep("^(L|l)ive|LD|(V|v)iability$", markernames(gs), value = TRUE)

  if (length(live) == 0) {
    catf(">> There is no viability dye channel in this gating set...")
    return(NULL)
  } else if (length(live) > 1) {
    catf(">> There are more than one viability dye channel in this gating set...")
    catf(paste(live, collapse = ", "))
  }

  live
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
      collapseDataForGating = collapseDataForGating,
      mc.cores = detect_cores(),
      parallel_type = "multicore"
    )
  }
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
