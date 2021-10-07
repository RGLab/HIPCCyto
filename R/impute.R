#' @importFrom assertthat assert_that is.dir is.writeable not_empty
#' @importFrom flowWorkspace load_gs pData sampleNames gs_pop_get_gate pData<- recompute save_gs
impute_gates <- function(gs_dir, samples_to_impute, by_batch = TRUE, method = "nearest") {
  assert_that(is.dir(gs_dir))
  assert_that(is.writeable(gs_dir))
  assert_that(dir.exists(file.path(gs_dir, "gs")))

  gs <- load_gs(file.path(gs_dir, "gs"))
  pd <- pData(gs)

  assert_that(is.character(samples_to_impute))
  assert_that(not_empty(samples_to_impute))
  assert_that(all(samples_to_impute %in% sampleNames(gs)))
  assert_that(is.logical(by_batch))
  if (isTRUE(by_batch)) {
    assert_that(!is.null(pd$batch), msg = "There is no batch group. Try again with `by_batch = FALSE`.")
  }
  assert_that(method %in% c("nearest", "consensus"), msg = "Use 'nearest' or 'consensus'.")

  gate <- "Lymphocytes"

  # save original gates, plot, and QC report
  catf("Backing up original gates and QC files")
  original_gates <- gs_pop_get_gate(gs, gate)
  saveRDS(original_gates, file.path(gs_dir, sprintf("%s_original.RDS", gate)))
  file.copy(file.path(gs_dir, "markers.png"), file.path(gs_dir, "markers_original.png"))
  file.copy(file.path(gs_dir, "QC.html"), file.path(gs_dir, "QC_original.html"))

  pd$to_impute <- pd$name %in% samples_to_impute
  if (isTRUE(by_batch)) {
    batch_groups <- split(pd, pd$batch)
  } else {
    batch_groups <- list(pd)
  }

  # impute gates
  catf(sprintf("Imputing %s gate by '%s' method", gate, method))
  imputed <- lapply(batch_groups, function(batch_pd) {
    if (!is.null(batch_pd$batch)) {
      catf(sprintf("batch: %s", unique(batch_pd$batch)))
    }

    if (all(batch_pd$to_impute == TRUE)) {
      catf("No template gates remaining to use for imputation")
      return()
    } else if (all(batch_pd$to_impute == FALSE)) {
      catf("No samples to impute")
      return()
    }

    batch_gs <- gs[batch_pd$name]
    batch_samples_to_impute <- batch_pd$name[batch_pd$name %in% samples_to_impute]
    if (method == "consensus") {
      impute_by_consensus(batch_gs, batch_samples_to_impute, gate)
    } else { # "nearest" as default
      impute_by_nearest(batch_gs, batch_samples_to_impute, gate)
    }
  })
  pd$imputed <- pd$name %in% unlist(imputed)
  pData(gs) <- pd
  catf(sprintf("%s out of %s samples were imputed in %s gate", sum(pd$imputed), sum(pd$to_impute), gate))
  catf("Recomputing the cell events")
  recompute(gs)

  # re-save gating set and QC files
  catf("Saving new gates and creating new QC files")
  save_gs(gs, path = file.path(gs_dir, "gs"))
  render_comparison_report(gs_dir)
  render_qc_report(gs_dir)

  gs
}

#' @importFrom flowWorkspace gs_pop_get_gate gs_pop_set_gate
#' @importFrom methods is
#' @importFrom flowCore ellipsoidGate
impute_by_consensus <- function(gs, samples_to_impute, gate) {
  # Get the template gates
  good_gates <- gs_pop_get_gate(gs[samples_to_impute], gate)

  if (is(good_gates[[1]], "ellipsoidGate")) {
    # Getting consensus center:
    good_centers <- do.call(rbind, lapply(good_gates, function(gate) {
      gate@mean
    }))
    consensus_location <- colMeans(good_centers)

    # Getting consensus Mahalanobis distance:
    consensus_distance <- mean(sapply(good_gates, function(gate) gate@distance))

    # Getting consensus shape (covariance matrix):
    good_eigs <- lapply(good_gates, function(gate) {
      eigen(gate@cov)
    })
    good_eigVals <- do.call(rbind, lapply(good_eigs, function(this_eigs) {
      this_eigs$values
    }))
    good_radii <- sqrt(good_eigVals)
    consensus_radii <- colMeans(good_radii)
    consensus_eigVals <- consensus_radii^2

    # Get consensus principal eigenvector as vector sum
    consensus_eigVec1 <- colMeans(do.call(rbind, lapply(good_eigs, function(this_eigs) {
      this_eigs$vectors[, 1]
    })))

    # Just get the second consensus eigenvector from the normal to the first
    consensus_eigVec2 <- c(consensus_eigVec1[2], -consensus_eigVec1[1])
    magnitude <- sqrt(sum(consensus_eigVec1^2))

    # Renormalize vector sums
    consensus_eigVec1 <- consensus_eigVec1 / magnitude
    consensus_eigVec2 <- consensus_eigVec2 / magnitude

    consensus_eigVecs <- cbind(consensus_eigVec1, consensus_eigVec2)
    consensus_cov <- consensus_eigVecs %*% diag(consensus_eigVals) %*% solve(consensus_eigVecs)
    rownames(consensus_cov) <- colnames(consensus_cov) <- names(consensus_location)

    # Wrapping it together in consensus ellipsoidGate object:
    consensus_eg <- ellipsoidGate(gate = consensus_cov, mean = consensus_location, distance = consensus_distance, filterId = "consensus_ellipsoidGate")
    imputed_gates <- lapply(seq_along(samples_to_impute), function(dummy) consensus_eg)
    names(imputed_gates) <- samples_to_impute
    imputed_gates
  } else {
    warning("Unsupported gate type for imputation.")
    return()
  }

  gs_pop_set_gate(gs[samples_to_impute], gate, imputed_gates)

  samples_to_impute
}

#' @importFrom flowIncubator nearestSamples regateNearestSamples
impute_by_nearest <- function(gs, samples_to_impute, gate) {
  passed <- sampleNames(gs)[!sampleNames(gs) %in% samples_to_impute]
  res <- nearestSamples(gs, node = gate, failed = samples_to_impute, passed = passed)
  regateNearestSamples(gs, res, gate)

  samples_to_impute
}

#' @importFrom flowWorkspace load_gs gs_pop_set_gate pData<- save_gs
recover_original_gates <- function(gs_dir) {
  gate <- "Lymphocytes"
  gs <- load_gs(file.path(gs_dir, "gs"))
  gates <- readRDS(file.path(gs_dir, sprintf("%s_original.RDS", gate)))
  gs_pop_set_gate(gs, gate, gates)
  pData(gs)$to_impute <- NULL
  pData(gs)$imputed <- NULL

  catf("Recovering origianl gates and QC files")
  save_gs(gs, path = file.path(gs_dir, "gs"))
  file.copy(file.path(gs_dir, "markers_original.png"), file.path(gs_dir, "markers.png"))
  file.copy(file.path(gs_dir, "QC_original.html"), file.path(gs_dir, "QC.html"))

  gs
}

#' @importFrom rmarkdown render
render_comparison_report <- function(gs_dir) {
  catf("Compiling gate comparison report")
  file_path <- file.path(gs_dir, "comparison.html")

  catf(sprintf("output_file: %s", file_path))
  input <- file.path(gs_dir, "comparison.Rmd")
  file.copy(system.file("qc/comparison.Rmd", package = "HIPCCyto"), input, overwrite = TRUE)
  render(
    input = input,
    output_file = file_path,
    params = list(gs_dir = gs_dir)
  )

  file_path
}
