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
  print(system.time(cs <- create_cytoset(files$filePath, study, debug_dir)))

  # merge metadata
  print(system.time(cs <- merge_metadata(cs, files, study, debug_dir)))

  # merge batch information
  print(system.time(cs <- merge_batch(cs, study, debug_dir)))

  # create a gating set
  print(system.time(gs <- create_gs(cs, study, debug_dir)))

  # pre-process
  print(system.time(gs <- standardize_markernames(gs, study, debug_dir)))
  print(system.time(gs <- compensate_gs(gs, study, debug_dir)))
  print(system.time(gs <- transform_gs(gs, study, debug_dir)))

  # gate
  gate_gs(gs, study, debug_dir)

  gs
}


# processing functions ---------------------------------------------------------
#' @importFrom flowWorkspace load_cytoset_from_fcs cytoset
#' @importFrom cytoqc cqc_load_fcs cqc_check cqc_match cqc_match_update cqc_match_remove cqc_fix
create_cytoset <- function(filePath, study, debug_dir = NULL) {
  catf(">> Reading files and creating a cytoset...")

  cs <- suppressMessages(cqc_load_fcs(
    filePath,
    num_threads = detect_cores(),
    is_h5 = TRUE
  ))
  channel_check <- cqc_check(cs, "channel")

  # Check if inconsistent
  if (length(unique(channel_check$group_id)) > 1) {
    # Get custom control of channel reference, re-mapping, and fuzzy-matching control from study info in DATA
    study_info <- DATA[[study]]
    max.distance <- study_info$max.distance
    channel_ref <- study_info$channel_ref
    map <- study_info$map

    channel_match <- custom_match_cytoset(channel_check, max.distance, channel_ref, map)

    # We need to handle possibility of extra channels in reference. By design, cytoqc will not automatically delete these
    # but instead will just throw a warning. For our purposes, if there are extra channels in the reference (e.g. Time, extra scatter channels),
    # we can explicitly delete them to ensure consistency for the resulting cytoset
    missing_channels <- do.call(c, lapply(channel_match$match_result, function(group) {group$missing}))
    if(length(missing_channels) > 0){
      # Drop those extra channels from the reference
      channel_ref <- channel_ref[!channel_ref %in% missing_channels]
      # And re-run the match (now the suggested fix will delete them)
      channel_match <- custom_match_cytoset(channel_check, max.distance, channel_ref, map)
    }

    cqc_fix(channel_match)
  }
  # cs with consistent channels
  cs <- cytoset(cs)
  save_debug(cs, "create_cs", debug_dir)

  cs
}

# A simple wrapper to handle the HIPCCyto matching logic of:
# 1) Specification of reference channels manually or automatically by most abundant group in check
# 2) Automatic match with optional fuzziness by max.distance
# 3) Manual updates to override automatic match using map
custom_match_cytoset <- function(check_result, max.distance, channel_ref, map){
  # 1) If not specified, use the panel with the greatest consensus (most abundant group in check)
  if (is.null(channel_ref))
    channel_ref <- colnames(cs[[as.data.frame(channel_check)[which.max(channel_check$nObject), "object"]]])

  # If not specified, no fuzzy match
  if (is.null(max.distance))
    max.distance <- 0.0
  # 2) First try automatic match
  channel_match <- cqc_match(channel_check, ref = channel_ref, max.distance = max.distance)

  # 3) Allow manual updating to override automatic match
  if (!is.null(map)) {
    # Remove any existing match
    tryCatch(channel_match <- cqc_match_remove(channel_match, map$channels), error = function(e) {})
    update_ref <- map$alias
    names(update_ref) <- map$channels
    channel_match <- cqc_match_update(channel_match, map = update_ref)
  }
  channel_match
}

#' @importFrom flowWorkspace phenoData phenoData<- cf_keyword_insert
merge_metadata <- function(cs, files, study, debug_dir = NULL) {
  catf(">> Merging metedata...")
  phenoData(cs)$study_accession <- files[phenoData(cs)$name, ]$studyAccession
  phenoData(cs)$participant_id <- files[phenoData(cs)$name, ]$subjectAccession
  phenoData(cs)$age_reported <- files[phenoData(cs)$name, ]$ageEvent
  phenoData(cs)$gender <- files[phenoData(cs)$name, ]$gender
  phenoData(cs)$race <- files[phenoData(cs)$name, ]$race
  phenoData(cs)$study_time_collected <- files[phenoData(cs)$name, ]$studyTimeCollected
  phenoData(cs)$study_time_collected_unit <- files[phenoData(cs)$name, ]$studyTimeCollectedUnit
  phenoData(cs)$file_info_name <- files[phenoData(cs)$name, ]$fileName
  phenoData(cs)$description <- files[phenoData(cs)$name, ]$fileDetail
  phenoData(cs)$type <- files[phenoData(cs)$name, ]$biosampleType
  phenoData(cs)$subtype <- files[phenoData(cs)$name, ]$biosampleSubtype
  phenoData(cs)$cohort <- files[phenoData(cs)$name, ]$armName

  ver <- get_version()
  dr <- get_dr()
  for (i in seq_along(cs)) {
    cf_keyword_insert(cs[[i, returnType = "cytoframe"]], "HIPCCyto_version", ver)
    cf_keyword_insert(cs[[i, returnType = "cytoframe"]], "ImmPort_data_release", dr)
  }

  save_debug(cs, "merge_metadata", debug_dir)

  cs
}

#' @importFrom flowCore description
merge_batch <- function(cs, study, debug_dir = NULL) {
  catf(">> Merging batch column...")
  keyword <- DATA[[study]]$batch

  if (!is.null(keyword)) {
    phenoData(cs)$batch <- unlist(
      lapply(
        cs,
        FUN = function(x) {
          val <- description(x)[keyword][[1]]
          if (is.null(val)) val <- NA
          val
        }
      )[phenoData(cs)$name],
      use.names = FALSE
    )
  }

  save_debug(cs, "merge_batch", debug_dir)

  cs
}

#' @importFrom flowWorkspace GatingSet
create_gs <- function(cs, study, debug_dir = NULL) {
  catf(">> Creating a gating set...")
  gs <- GatingSet(cs)

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

#' @importFrom flowWorkspace gs_pop_get_data compensate lapply
#' @importFrom flowCore spillover
compensate_gs <- function(gs, study, debug_dir = NULL) {
  catf(">> Applying compensation...")
  cs <- gs_pop_get_data(gs)
  cols <- colnames(gs)
  comp <- lapply(cs, function(x) {
    spills <- spillover(x)
    spill <- spills[!sapply(spills, is.null)][[1]] # pick the first non-empty matrix
    keep <- colnames(spill) %in% cols
    spill[keep, keep] # remove extra channels
  })

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


# helper functions -------------------------------------------------------------
#' @importFrom flowWorkspace cs_get_h5_file_path save_cytoset
#' @importFrom flowWorkspace save_gs
save_debug <- function(obj, func, debug_dir = NULL) {
  if (!is.null(debug_dir)) {
    path <- tempfile(paste0(func, "_"), debug_dir)
    catf(sprintf(">> Storing intermediate file to %s for debugging...", path))
    if (is(obj, "cytoset")) {
      save_cytoset(obj, path)
    } else if (is(obj, "GatingSet")) {
      save_gs(obj, path, backend_opt = "copy")
    } else {
      saveRDS(obj, path)
    }
  }
}
