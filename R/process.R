#' @importFrom flowWorkspace GatingSet
#' @export
process_study <- function(study, input_dir) {
  # summariz files
  files <- summarize_study(study, input_dir)

  message("Creating gating set for each panel...")
  gsl <- lapply(unique(files$panel), function(panel) {
    message(panel)

    # load files
    nc <- create_nc(files$filePath[files$panel == panel], study)

    # merge metadata
    nc <- merge_metadata(nc, files)

    # merge acquired date (to examine batches)
    # there are two batches, but don't correspond with either cohort or
    # study_time_collected
    nc <- merge_batch(nc, study)

    # create a gating set
    message("Creating a gating set...")
    gs <- GatingSet(nc)

    # pre-process
    gs <- standardize_markernames(gs)
    gs <- compensate_gs(gs)
    gs <- transform_gs(gs)

    # gate
    gate_gs(gs, study)
  })

  gsl
}


#' @importFrom ImmPortR query_filePath
#' @importFrom flowCore read.FCSheader
summarize_study <- function(study, input_dir) {
  files <- query_filePath(study)
  files <- files[files$fileDetail == "Flow cytometry result", ]
  files$filePath <- file.path(input_dir, files$fileName)
  rownames(files) <- files$fileName

  # check if files exist. If not, throw warning and fetch them
  file.exists(files$filePath)

  # summarize files and
  # read headers and summarize panels
  map <- DATA[[study]]$map
  message("Reading in fcs headers...")
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

    paste(gtools::mixedsort(paste0(PNS, " (", PNN, ")")), collapse = "; ")
  })

  files$panel <- panels
  ps <- unique(panels)

  # how many gating sets will be created
  # by panels, sample type, measurement technique, experiment accession
  message(sprintf("There are %s panel(s)...", length(ps)))
  x <- lapply(strsplit(ps, "; "), function(x) {
    message(paste(gsub(" \\S+", "", x), collapse = " "))
  })

  files
}

#' @importFrom ncdfFlow read.ncdfFlowSet
#' @importFrom parallel detectCores
create_nc <- function(filePath, study) {
  message("Reading files and creating a flow set...")
  map <- DATA[[study]]$map

  read.ncdfFlowSet(
    filePath,
    channel_alias = map,
    mc.cores = detectCores()
  )
}

#' @importFrom ncdfFlow phenoData phenoData<-
merge_metadata <- function(nc, files) {
  message("Merging metedata...")
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

  nc
}

#' @importFrom flowCore fsApply description
merge_batch <- function(nc, study) {
  message("Merging batch column...")
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

  nc
}

#' @importFrom flowWorkspace markernames markernames<-
standardize_markernames <- function(gs) {
  message("Standardizing marker names...")

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

  gs
}

#' @importFrom flowWorkspace getData compensate
#' @importFrom flowCore spillover
compensate_gs <- function(gs) {
  message("Applying compensation...")
  nc <- getData(gs)
  cols <- colnames(gs)
  comp <- fsApply(nc, function(x) {
    spills <- spillover(x)
    spill <- spills[!sapply(spills, is.null)][[1]] # pick the first non-empty matrix
    keep <- colnames(spill) %in% cols
    spill[keep, keep] # remove extra channels
  }, simplify = FALSE)

  compensate(gs, comp)
}

# transform fluoresence channels with biexponential transformation
#' @importFrom flowWorkspace colnames transform
#' @importFrom flowCore estimateLogicle
transform_gs <- function(gs) {
  message("Applying transformation...")
  channels <- colnames2(gs)
  trans <- estimateLogicle(gs[[1]], channels)

  transform(gs, trans)
}

#' @importFrom openCyto gatingTemplate gating
gate_gs <- function(gs, study) {
  message("Gating...")
  filePath <- sprintf("extdata/gating_template/%s.csv", study)
  file <- system.file(filePath, package = "HIPCCyto")

  if (file != "") {
    gt <- gatingTemplate(file)
    gating(gt, gs, mc.cores = detectCores(), parallel_type = "multicore")
  }

  gs
}

#' @importFrom flowCore parameters
colnames2 <- function(gs) {
  # grep("SC-|Time", colnames(gs), invert = TRUE, value = TRUE) # can't use this for now

  channels <- parameters(getData(gs)[[1]])@data$name
  markers <- parameters(getData(gs)[[1]])@data$desc

  unname(channels[!is.na(markers)])
}
