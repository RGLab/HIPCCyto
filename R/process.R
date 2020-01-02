#' @export
process_study <- function(study, input_dir, debug_dir = NULL) {
  # summarize files
  files <- summarize_study(study, input_dir, debug_dir)
  files_by_panel <- split(files, files$panel)

  # create gating set for each panel
  lapply(files_by_panel, process_panel, debug_dir = debug_dir)
}


#' @importFrom ImmPortR query_filePath
#' @importFrom flowCore read.FCSheader
#' @importFrom gtools mixedsort
summarize_study <- function(study, input_dir, debug_dir = NULL) {
  files <- query_filePath(study)
  files <- files[files$fileDetail == "Flow cytometry result", ]
  files <- files[grepl(".fcs$", files$fileName), ]
  files$filePath <- file.path(input_dir, files$fileName)
  # rownames(files) <- files$fileName

  # check if files exist. If not, throw warning and fetch them
  fcs_exist <- file.exists(files$filePath)
  if (any(!fcs_exist)) {
    stop("These fcs files do not exist at ", input_dir, "\n", paste(files[!fcs_exist, "fileName"], collapse = "\n"))
  }

  # summarize files
  message(sprintf(">> There are %s fcs files in %s...", nrow(files), study))

  # read headers and summarize panels
  map <- DATA[[study]]$map
  message(">> Reading in fcs headers...")
  headers <- lapply(files$filePath, function(file) read.FCSheader(file, channel_alias = map))
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
    marker_exist <- !is.na(PNS) & PNS %in% names(MARKERS)
    PNS[marker_exist] <- MARKERS[PNS[marker_exist]]

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
  message(sprintf(">> There are %s panel(s)...", length(ps)))
  panels_clean <- sapply(strsplit(ps, "; "), function(x) {
    if (length(x) == 0) {
      ""
    } else {
      paste(gsub(" \\(.+\\)$", "", x), collapse = " | ")
    }
  })
  message(paste(mixedsort(panels_clean), collapse = "\n"))

  save_debug(files, "summarize_study", debug_dir)

  files
}


process_panel <- function(files, debug_dir = NULL) {
  study <- unique(files$studyAccession)
  panel <- unique(files$panel)
  stopifnot(length(study) == 1, length(panel) == 1)

  message(paste(rep("=", times = 80), collapse = ""))
  message(sprintf(">> Processing %s fcs files for this panel...", nrow(files)))
  message(paste(strsplit(panel, split = "; ")[[1]], collapse = "\n"))

  # load files
  nc <- create_nc(files$filePath, study, debug_dir)

  # merge metadata
  nc <- merge_metadata(nc, files, study, debug_dir)

  # merge batch information
  nc <- merge_batch(nc, study, debug_dir)

  # create a gating set
  gs <- create_gs(nc, study, debug_dir)

  # pre-process
  gs <- standardize_markernames(gs, study, debug_dir)
  gs <- compensate_gs(gs, study, debug_dir)
  gs <- transform_gs(gs, study, debug_dir)

  # gate
  gate_gs(gs, study, debug_dir)
}


#' @importFrom ncdfFlow read.ncdfFlowSet
#' @importFrom parallel detectCores
create_nc <- function(filePath, study, debug_dir = NULL) {
  message(">> Reading files and creating a flow set...")
  map <- DATA[[study]]$map

  nc <- suppressMessages(read.ncdfFlowSet(
    filePath,
    channel_alias = map,
    mc.cores = detectCores()
  ))

  save_debug(nc, "create_nc", debug_dir)

  nc
}

#' @importFrom ncdfFlow phenoData phenoData<-
merge_metadata <- function(nc, files, study, debug_dir = NULL) {
  message(">> Merging metedata...")
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
  message(">> Merging batch column...")
  keyword <- DATA[[study]]$batch

  if (!is.null(keyword)) {
    phenoData(nc)$batch <- unlist(
      fsApply(
        x = nc,
        FUN = function(x) description(x)[keyword],
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
  message(">> Creating a gating set...")
  gs <- GatingSet(nc)

  save_debug(gs, "create_gs", debug_dir)

  gs
}

#' @importFrom flowWorkspace markernames markernames<-
standardize_markernames <- function(gs, study, debug_dir = NULL) {
  message(">> Standardizing marker names...")

  # get current marker names and name them with channel names
  names_gs <- markernames(gs)
  if (is(names_gs, "list")) {
    names_gs <- names_gs[[1]]
  }
  names(names_gs) <- colnames2(gs)

  # retrieve standard marker names from MARKERS
  marker_exist <- names_gs %in% names(MARKERS)
  standards <- MARKERS[names_gs[marker_exist]]
  names_gs[marker_exist] <- standards

  # assign the fixed marker names
  standards <- standards[standards != names(standards)]
  message(paste(paste(names(standards), standards, sep = " -> "), collapse = "\n"))
  markernames(gs) <- names_gs

  save_debug(gs, "standardize_markernames", debug_dir)

  gs
}

#' @importFrom flowWorkspace getData compensate
#' @importFrom flowCore spillover
compensate_gs <- function(gs, study, debug_dir = NULL) {
  message(">> Applying compensation...")
  nc <- getData(gs)
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
  message(">> Applying transformation...")
  channels <- colnames2(gs)
  if (length(channels) > 0) {
    trans <- estimateLogicle(gs[[1]], channels)
    gs <- transform(gs, trans)
  }

  save_debug(gs, "transform_gs", debug_dir)

  gs
}

#' @importFrom openCyto gatingTemplate gating add_pop remove_pop
#' @importFrom flowClust flowClust
gate_gs <- function(gs, study, debug_dir = NULL) {
  filePath <- sprintf("extdata/gating_template/%s.csv", study)
  file <- system.file(filePath, package = "HIPCCyto")

  if (file != "") {
    message(sprintf(">> Applying gating template (%s)...", file))
    gt <- gatingTemplate(file)
    gating(gt, gs, mc.cores = detectCores(), parallel_type = "multicore")
  } else {
    message(">> Gating template does not exist for this study...")
    message(">> Applying default gating methods...")
    apply_quadrant_gate(gs)
    apply_singlet_gate(gs, "FSC")
    apply_singlet_gate(gs, "SSC")
    apply_live_gate(gs)
    apply_lymphocyte_gate(gs)
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
    message(sprintf(">> Storing intermediate file to %s for debugging...", path))
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
colnames2 <- function(gs) {
  # can't use this for now
  # grep("SC-|Time", colnames(gs), invert = TRUE, value = TRUE)

  channels <- parameters(getData(gs)[[1]])@data$name
  markers <- parameters(getData(gs)[[1]])@data$desc

  unname(channels[!is.na(markers)])
}

#' @importFrom flowWorkspace gs_get_pop_paths
get_parent <- function(gs) {
  rev(gs_get_pop_paths(gs, path = 1))[1]
}

#' @importFrom flowWorkspace gs_pop_get_data
compute_K <-function(gs) {
  message(">> Computing for optimal number of clusters (K) for each sample...")
  nc <- gs_pop_get_data(gs, get_parent(gs))
  Ks <- mclapply(sampleNames(nc), function(x) {
    flowClust(nc[[x]]@exprs[, c("FSC-A", "SSC-A")], K = 1:5, criterion = "ICL")@index
  }, mc.cores = detectCores())
  print(table(unlist(Ks)))

  unique(unlist(Ks))[which.max(tabulate(match(unlist(Ks), unique(unlist(Ks)))))]
}

#' @importFrom flowWorkspace gs_pop_get_gate
compute_target <- function(gs) {
  message(">> Computing the target location of the lymphocyte cluster...")
  lymphocytes_before <- gs_pop_get_gate(gs, "Lymphocytes")
  fsc_before <- sapply(lymphocytes_before, function(x) x@mean[1])
  ssc_before <- sapply(lymphocytes_before, function(x) x@mean[2])
  density_fsc_before <- density(fsc_before)
  density_ssc_before <- density(ssc_before)

  c(
    FSC = density_fsc_before$x[which.max(density_fsc_before$y)],
    SSC = density_ssc_before$x[which.max(density_ssc_before$y)]
  )
}

#' @importFrom flowCore rectangleGate
.quadrantGate <- function(fr, pp_res, channels = NA, filterId = "", quadrant, ...) {
  range <- switch(
    quadrant,
    "1" = c(0, Inf, 0, Inf),
    "2" = c(-Inf, 0, 0, Inf),
    "3" = c(-Inf, 0, -Inf, 0),
    "4" = c(0, Inf, -Inf, 0)
  )
  mat <- matrix(range, ncol = 2, dimnames = list(c("min", "max"), c(channels[1], channels[2])))
  rectangleGate(filterId = filterId, .gate = mat)
}

#' @importFrom openCyto registerPlugins
apply_quadrant_gate <- function(gs) {
  message(">> Applying quadrant gate...")
  registerPlugins(fun = .quadrantGate, methodName = "quadrantGate")
  add_pop(
    gs = gs,
    alias = "pos",
    pop = "+",
    parent = get_parent(gs),
    dims = "FSC-A,SSC-A",
    gating_method = "quadrantGate",
    gating_args = "quadrant = 1"
  )
}

apply_singlet_gate <- function(gs, channel) {
  alias <- sprintf("SC%s", channel)
  A <- sprintf("%s-A", channel)
  H <- sprintf("%s-H", channel)
  if (H %in% colnames(gs)) {
    message(sprintf(">> Applying singlet gate by scatter channel (%s)...", alias))
    add_pop(
      gs = gs,
      alias = alias,
      pop = "+",
      parent = get_parent(gs),
      dims = sprintf("%s,%s", A, H),
      gating_method = "singletGate",
      gating_args = "prediction_level = 0.99",
      mc.cores = detectCores(),
      parallel_type = "multicore"
    )
  }
}

apply_lymphocyte_gate <- function(gs) {
  K <- compute_K(gs)
  gate_lymphocytes(gs, K)

  target <- compute_target(gs)

  message(">> Removing lymphocyte gate...")
  remove_pop(gs)

  gate_lymphocytes(gs, K, target)
}

gate_lymphocytes <- function(gs, K, target = NULL) {
  if (is.null(target)) {
    gating_args <- sprintf("randomStart = 1, K = %s", K)
  } else {
    gating_args <- sprintf("randomStart = 1, K = %s, target = c(%s, %s)", K, target["FSC"], target["SSC"])
  }

  message(">> Applying lymphocytes gate with flowClust by forward and side scatters (Lymphocytes)...")
  message(gating_args)
  add_pop_lymph(gs, gating_args)
}

get_live_marker <- function(gs) {
  live <- grep("^(L|l)ive|LD|(V|v)iability$", markernames(gs), value = TRUE)

  if (length(live) == 0) {
    message(">> There is no viability dye channel in this gating set...")
    return(NULL)
  } else if (length(live) > 1) {
    message(">> There are more than one viability dye channel in this gating set...")
    message(paste(live, collapse = ", "))
  }

  live
}

apply_live_gate <- function(gs) {
  live <- get_live_marker(gs)
  if (!is.null(live)) {
    message(sprintf(">> Applying live/dead gate with mindensity by %s (Live)...", live))
    collapseDataForGating <- !is.null(pData(gs)$batch)
    groupBy <- ifelse(collapseDataForGating, "batch", NA)
    add_pop(
      gs = gs,
      alias = "Live",
      pop = "-",
      parent = get_parent(gs),
      dims = live,
      gating_method = "tailgate",
      groupBy = groupBy,
      collapseDataForGating = collapseDataForGating,
      mc.cores = detectCores(),
      parallel_type = "multicore"
    )
  }
}

add_pop_lymph <- function(gs, gating_args) {
  add_pop(
    gs = gs,
    alias = "Lymphocytes",
    pop = "+",
    parent = get_parent(gs),
    dims = "FSC-A,SSC-A",
    gating_method = "flowClust",
    gating_args = gating_args,
    mc.cores = detectCores(),
    parallel_type = "multicore"
  )
}
