process_SDY820 <- function(input_dir) {
  # create a flowset
  files <- ImmPortR::query_filePath("SDY820")
  files <- files[files$fileDetail == "Flow cytometry result", ]
  files$filePath <- file.path(input_dir, gsub("/SDY820/", "", files$filePath))
  map <- data.frame(
    alias = c("APC-eFluor 780-A", "eFluor 450-A"),
    channels = c("APC-eFluor780-A", "eFluor450-A")
  )
  nc <- ncdfFlow::read.ncdfFlowSet(
    files = files$filePath,
    channel_alias = map,
    mc.cores = parallel::detectCores()
  )

  # merge metadata
  data.table::setDT(files, fileName)
  ncdfFlow::phenoData(nc)$participant_id <- files[ncdfFlow::phenoData(nc)$name, subjectAccession]
  ncdfFlow::phenoData(nc)$age_reported <- files[ncdfFlow::phenoData(nc)$name, ageEvent]
  ncdfFlow::phenoData(nc)$gender <- files[ncdfFlow::phenoData(nc)$name, gender]
  ncdfFlow::phenoData(nc)$race <- files[ncdfFlow::phenoData(nc)$name, race]
  ncdfFlow::phenoData(nc)$study_time_collected <- files[ncdfFlow::phenoData(nc)$name, studyTimeCollected]
  ncdfFlow::phenoData(nc)$study_time_collected_unit <- files[ncdfFlow::phenoData(nc)$name, studyTimeCollectedUnit]
  ncdfFlow::phenoData(nc)$file_info_name <- files[ncdfFlow::phenoData(nc)$name, fileName]
  ncdfFlow::phenoData(nc)$description <- files[ncdfFlow::phenoData(nc)$name, fileDetail]
  ncdfFlow::phenoData(nc)$type <- files[ncdfFlow::phenoData(nc)$name, biosampleType]
  ncdfFlow::phenoData(nc)$subtype <- files[ncdfFlow::phenoData(nc)$name, biosampleSubtype]
  ncdfFlow::phenoData(nc)$cohort <- files[ncdfFlow::phenoData(nc)$name, armName]

  # merge acquired date (to examine batches)
  # there are two batches, but don't correspond with either cohort or study time collected
  ncdfFlow::phenoData(nc)$acquired_date <- unlist(flowCore::fsApply(nc, function(x) flowCore::description(x)$`CST SETUP DATE`, simplify = FALSE)[ncdfFlow::phenoData(nc)$name])

  # create a gating set
  gs <- flowWorkspace::GatingSet(nc)

  # standardize marker names
  markers <- c(
    "APC-A" = "CD14",
    "Alexa Fluor 700-A" = "CD3E",
    "APC-eFluor 780-A" = "CD4",
    "V500-A" = "LD",
    "FITC-A" = "IL2RA",
    "PE-A" = "NCAM1",
    "PerCP-Cy5-5-A" = "CCR7",
    "PE-Cy7-A" = "CD19",
    "eFluor 450-A" = "PTPRC",
    "BV650-A" = "CD8A"
  )
  flowWorkspace::markernames(gs) <- markers

  # compensate
  comp <- flowCore::fsApply(flowWorkspace::getData(gs), function(x) flowCore::spillover(x)$SPILL, simplify = FALSE)
  gs <- flowWorkspace::compensate(gs, comp)

  # transform fluoresence channels with biexponential transformation
  channels <- grep("SC-|Time", flowWorkspace::colnames(gs), invert = TRUE, value = TRUE)
  trans <- flowCore::estimateLogicle(gs[[1]], channels)
  gs <- flowWorkspace::transform(gs, trans)

  # add gates manually
  openCyto::add_pop(gs = gs, alias = "Cells", pop = "+", parent = "root", dims = "FSC-A,SSC-A", gating_method = "singletGate", gating_args = "prediction_level = 0.99, maxit = 100", mc.cores = parallel::detectCores(), parallel_type = "multicore")
  openCyto::add_pop(gs = gs, alias = "SCFSC", pop = "+", parent = "Cells", dims = "FSC-A,FSC-H", gating_method = "singletGate", gating_args = "prediction_level = 0.99, maxit = 100", mc.cores = parallel::detectCores(), parallel_type = "multicore")
  openCyto::add_pop(gs = gs, alias = "SCSSC", pop = "+", parent = "SCFSC", dims = "SSC-A,SSC-H", gating_method = "singletGate", gating_args = "prediction_level = 0.99, maxit = 100", mc.cores = parallel::detectCores(), parallel_type = "multicore")
  openCyto::add_pop(gs = gs, alias = "Live", pop = "-", parent = "SCSSC", dims = "V500-A", gating_method = "mindensity", collapseDataForGating = TRUE, groupBy = "acquired_date", mc.cores = parallel::detectCores(), parallel_type = "multicore")

  # save gating set
  output_dir <- file.path(input_dir, "gs")
  flowWorkspace::save_gs(gs, output_dir)
}
