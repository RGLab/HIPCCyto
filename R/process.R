#' @export
process_study <- function(study, input_dir, debug_dir = NULL) {
  # summarize files
  files <- summarize_study(study, input_dir)
  files_by_panel <- split(files, files$panel)

  # create gating set for each panel
  lapply(files_by_panel, process_panel, debug_dir = debug_dir)
}


#' @importFrom ImmPortR query_filePath
#' @importFrom flowCore read.FCSheader
#' @importFrom gtools mixedsort
summarize_study <- function(study, input_dir) {
  files <- query_filePath(study)
  files <- files[files$fileDetail == "Flow cytometry result", ]
  files$filePath <- file.path(input_dir, files$fileName)
  rownames(files) <- files$fileName

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

    paste(mixedsort(paste0(PNS, " (", PNN, ")")), collapse = "; ")
  })

  files$panel <- panels
  ps <- unique(panels)

  # how many gating sets will be created
  # by panels, sample type, measurement technique, experiment accession
  message(sprintf(">> There are %s panel(s)...", length(ps)))
  panels_clean <- sapply(strsplit(ps, "; "), function(x) {
    paste(gsub(" \\S+", "", x), collapse = " ")
  })
  message(paste(mixedsort(panels_clean), collapse = "\n"))

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
  trans <- estimateLogicle(gs[[1]], channels)

  gs <- transform(gs, trans)

  save_debug(gs, "transform_gs", debug_dir)

  gs
}

#' @importFrom openCyto gatingTemplate gating
gate_gs <- function(gs, study, debug_dir = NULL) {
  filePath <- sprintf("extdata/gating_template/%s.csv", study)
  file <- system.file(filePath, package = "HIPCCyto")

  if (file != "") {
    message(sprintf(">> Applying gating template (%s)...", file))
    gt <- gatingTemplate(file)
    gating(gt, gs, mc.cores = detectCores(), parallel_type = "multicore")
  }

  save_debug(gs, "gate_gs", debug_dir)

  gs
}


# HELPER FUNCTIONS -------------------------------------------------------------

#' @importFrom ncdfFlow getFileName save_ncfs
#' @importFrom flowWorkspace save_gs
save_debug <- function(obj, func, debug_dir = NULL) {
  if (!is.null(debug_dir)) {
    if (is(obj, "ncdfFlowSet")) {
      path <- file.path(debug_dir, paste0(func, "_", basename(getFileName(obj))))
      save_ncfs(obj, path, cdf = "copy")
    } else if (is(obj, "GatingSet")) {
      path <- file.path(debug_dir, paste0(func, "_", basename(getFileName(getData(obj)))))
      save_gs(obj, path, cdf = "copy")
    } else {
      path <- tempfile(paste0(func, "_"), debug_dir, ".RDS")
      saveRDS(obj, path)
    }
    message(sprintf(">> Storing intermediate file to %s for debugging...", path))
  }
}

#' @importFrom flowCore parameters
colnames2 <- function(gs) {
  # grep("SC-|Time", colnames(gs), invert = TRUE, value = TRUE) # can't use this for now

  channels <- parameters(getData(gs)[[1]])@data$name
  markers <- parameters(getData(gs)[[1]])@data$desc

  unname(channels[!is.na(markers)])
}
