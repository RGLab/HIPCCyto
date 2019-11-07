#' @importFrom flowWorkspace GatingSet
#' @export
process_study <- function(study, input_dir) {
  # summariz files
  files <- summarize_study(study, input_dir)

  # load files
  nc <- create_nc(files$filePath, DATA[[study]]$map)

  # merge metadata
  nc <- merge_metadata(nc, files)

  # merge acquired date (to examine batches)
  # there are two batches, but don't correspond with either cohort or
  # study_time_collected
  nc <- merge_batch(nc, DATA[[study]]$batch)

  # create a gating set
  gs <- GatingSet(nc)

  # pre-process
  gs <- standardize_markernames(gs, study)
  gs <- compensate_gs(gs)
  gs <- transform_gs(gs)

  # gate
  gate_gs(gs, study)
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
  headers <- lapply(files$filePath, function(file) read.FCSheader(file, channel_alias = map))
  panels <- sapply(headers, function(x) {
    header <- x[[1]]
    par <- as.integer(header["$PAR"])
    PNS <- unname(header[paste0("$P", seq_len(par), "S")])
    PNN <- unname(header[paste0("$P", seq_len(par), "N")])

    for (i in seq_len(nrow(map))) {
      PNN <- gsub(map$channels[i], map$alias[i], PNN)
    }

    # TODO: standardize marker names too

    PNN <- PNN[!is.na(PNS)]
    PNS <- PNS[!is.na(PNS)]

    paste(sort(paste(PNN, PNS, sep = ":")), collapse = "; ")
  })

  files$panel <- panels
  ps <- unique(panels)

  # how many gating sets will be created
  # by panels, sample type, measurement technique, experiment accession
  cat(sprintf("There are %s panels...\n", length(ps)))
  x <- lapply(strsplit(ps, "; "), function(x) {
    cat("-------------------------------------\n")
    cat(x, sep = "\n")
    cat("-------------------------------------\n")
  })

  files
}

#' @importFrom ncdfFlow read.ncdfFlowSet
#' @importFrom parallel detectCores
create_nc <- function(filePath, map) {
  message("Reading files and creating a flow set...")
  read.ncdfFlowSet(
    filePath,
    channel_alias = map,
    mc.cores = detectCores()
  )
}

#' @importFrom ncdfFlow phenoData
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
merge_batch <- function(nc, keyword) {
  message("Merging batch column...")
  phenoData(nc)$batch <- unlist(
    fsApply(
      x = nc,
      FUN = function(x) description(x)[keyword],
      simplify = FALSE
    )[phenoData(nc)$name],
    use.names = FALSE
  )

  nc

}

#' @importFrom flowWorkspace markernames
standardize_markernames <- function(gs, study) {
  message("Standardizing marker names...")
  markernames(gs) <- DATA[[study]]$markers

  gs
}

#' @importFrom flowWorkspace getData compensate
#' @importFrom flowCore spillover
compensate_gs <- function(gs) {
  message("Applying compensation...")
  nc <- getData(gs)
  comp <- fsApply(nc, function(x) spillover(x)$SPILL, simplify = FALSE)

  compensate(gs, comp)
}

# transform fluoresence channels with biexponential transformation
#' @importFrom flowWorkspace colnames transform
#' @importFrom flowCore estimateLogicle
transform_gs <- function(gs) {
  message("Applying transformation...")
  channels <- grep("SC-|Time", colnames(gs), invert = TRUE, value = TRUE)
  trans <- estimateLogicle(gs[[1]], channels)

  transform(gs, trans)
}

#' @importFrom openCyto gatingTemplate gating
gate_gs <- function(gs, study) {
  message("Gating...")
  filePath <- sprintf("extdata/gating_template/%s.csv", study)
  file <- system.file(filePath, package = "HIPCCyto")
  gt <- gatingTemplate(file)
  gating(gt, gs, mc.cores = detectCores(), parallel_type = "multicore")
}
