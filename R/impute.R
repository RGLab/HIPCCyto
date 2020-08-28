#' @importFrom flowWorkspace lapply
#' @importFrom flowIncubator nearestSamples regateNearestSamples
impute_gates <- function(gs, gate = "Lymphocytes", outliers = NULL, batch = NULL, method = "nearest", plot = FALSE) {
  to_impute <- outliers
  if (is.null(to_impute) || length(to_impute) <= 0) {
    message("No samples passed to to_impute arg of impute_gate.")
  } else if (is.numeric(to_impute)) {
    to_impute_idx <- to_impute
    to_impute_sn <- sampleNames(gs[to_impute])
  } else {
    to_impute_sn <- intersect(to_impute, sampleNames(gs))
    if (length(to_impute_sn) != length(to_impute)) {
      warning("Some sample names passed to to_impute arg of inpute_gate do not match GatingSet")
    }
    to_impute_idx <- match(to_impute_sn, sampleNames(gs))
  }

  # For flexibility of imputation level
  if (isTRUE(batch)) {
    batch <- "batch"
  } else if (isFALSE(batch)) {
    batch <- NULL
  }

  if (is.null(batch)) {
    batch_groups <- list(pData(gs))
  } else {
    batch_groups <- split(pData(gs), pData(gs)$batch)
  }

  for (idx in seq_along(batch_groups)) {
    batch_pd <- batch_groups[[idx]]
    to_impute_sn_batch <- to_impute_sn[to_impute_sn %in% rownames(batch_pd)]
    to_impute_idx_batch <- match(to_impute_sn_batch, sampleNames(gs))
    # Get the template gates
    good_gates <- lapply(gs[-to_impute_idx_batch], function(gh) {
      gh_pop_get_gate(gh, gate)
    })
    if (length(good_gates) <= 0) {
      warning("No template gates remaining in batch to use for imputation.")
      break
    }
    if (is(good_gates[[1]], "ellipsoidGate")) {
      if (method == "consensus") {
        # Building consensus ellipse

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
        imputed_gates <- lapply(seq_along(to_impute_sn_batch), function(dummy) consensus_eg)
        names(imputed_gates) <- to_impute_sn_batch

        if (plot) {
          # Currently ggcyto does not directly support gate-only plots without data
          # So this hacky workaround is necessary
          dummy <- matrix(as.numeric(c(0, 0)), nrow = 1, ncol = 2)
          gate_dims <- colnames(consensus_eg@cov)
          colnames(dummy) <- gate_dims
          fr <- flowFrame(dummy)

          # Scale the window to fit all gates
          max_bounds <- lapply(good_gates, function(gate) {
            pg <- as(gate, "polygonGate")
            maxima <- apply(pg@boundaries, 2, max)
          })
          max_bounds <- apply(do.call(rbind, max_bounds), 2, max)

          # Blank canvas
          combined <- ggcyto(fr, aes_(gate_dims[[1]], gate_dims[[2]])) +
            ggcyto_par_set(limits = list(x = c(0, max_bounds[[1]]), y = c(0, max_bounds[[2]])))

          # Add all the good template gates
          good_gate_geoms <- lapply((good_gates), geom_gate)
          combined <- Reduce(`+`, c(list(combined), good_gate_geoms))

          # Add the consensus gate and plot
          print(combined + geom_gate(consensus_eg, colour = "blue"))
        }

        gs_pop_set_gate(gs[to_impute_sn_batch], gate, imputed_gates)
      } else if (method == "nearest") {
        suppressWarnings(res <- nearestSamples(gs, node = gate, failed = to_impute_sn_batch))
        suppressWarnings(regateNearestSamples(gs, res, gate))
      } else {
        warning('Only "consensus" or "nearest" methods are currently implemented for imputing ellipsoidGates')
        break
      }
    } else {
      warning("Unsupported gate type for imputation.")
    }
  }

  recompute(gs)
}
