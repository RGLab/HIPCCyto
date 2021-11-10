#' Save gating sets and create QC reports and plots
#'
#' @param gsl A list. List of gating sets to save.
#' @param output_dir A character. Output directory.
#' @param qc A logical. Create QC reports/plots.
#'
#' @return A character. Output path.
#' @examples
#' \dontrun{
#' save_gating_sets(gsl, "SDY820/gs_dir")
#' }
#'
#' @importFrom flowWorkspace save_gs
#' @importFrom jsonlite toJSON
#' @importFrom ggcyto autoplot ggcyto
#' @importFrom ggplot2 ggsave facet_wrap facet_null aes_ aes
#' @importFrom flowCore getChannelMarker
#' @export
save_gating_sets <- function(gsl, output_dir, qc = TRUE) {
  gsl_path <- lapply(seq_len(length(gsl)), function(i) {
    catf(sprintf("Saving gating set #%s", i))
    gs <- gsl[[i]]
    gs_accession <- paste0("gs", i)

    path <- file.path(output_dir, gs_accession)
    gs_path <- file.path(path, "gs")
    if (dir.exists(gs_path)) {
      unlink(gs_path, recursive = TRUE)
    }
    dir.create(gs_path, showWarnings = FALSE, recursive = TRUE)

    catf(sprintf("gs_path = %s", gs_path))
    save_gs(gs, gs_path, overwrite = TRUE, backend_opt = "copy")

    if (isTRUE(qc)) {
      try(create_qc_files(gs, gs_accession, path))
      try(render_qc_report(path))
    }

    path
  })

  if (isTRUE(qc)) {
    try(render_study_report(output_dir))
  }

  gsl_path
}

create_qc_files <- function(gs, gs_accession, output_dir, full = FALSE) {
  catf("Creating QC summary and plots")

  save_gating_set_summary(gs, gs_accession, output_dir)
  save_spillover_heatmaps(gs, output_dir)
  save_gate_plots(gs, output_dir)
  save_marker_plots(gs, output_dir)

  if (isTRUE(full)) {
    for (node in gs_get_pop_paths(gs, path = 1)[-1]) {
      save_gate_plots_by_gate(gs, node, output_dir)
    }

    for (channel in colnames2(gs)) {
      marker <- markernames(gs)[channel]
      save_density_plots_by_marker(gs, marker, output_dir)
    }

    for (channel in colnames2(gs)) {
      marker <- markernames(gs)[channel]
      save_density_plots_by_gate(gs, marker, output_dir)
    }
  }

  output_dir
}


# Summarize --------------------------------------------------------------------

summarize_gating_set <- function(gs, gs_accession) {
  pd <- pData(gs)
  list(
    study_accession = unique(pd$study_accession),
    gs_accession = gs_accession,
    number_of_samples = length(gs),
    number_of_markers = length(markernames(gs)),
    markers = mixedsort(markernames(gs)),
    nodes = gs_get_pop_paths(gs, path = 1),
    sample_types = unique(pd$type),
    timepoints = unique(paste(pd$study_time_collected, pd$study_time_collected_unit)),
    cohorts = unique(pd$cohort),
    batches = as.list(table(pd$batch)),
    participants = unique(pd$participant_id),
    version = get_version(),
    commit_hash = get_commit_hash(),
    data_release = get_dr()
  )
}

save_gating_set_summary <- function(gs, gs_accession, output_dir) {
  catf("Saving gating set summary")

  summary <- summarize_gating_set(gs, gs_accession)
  json <- toJSON(summary, pretty = TRUE, auto_unbox = TRUE)
  file <- sprintf("%s/summary", output_dir)
  cat(json, file = file)

  file
}


# qc gate functions  -----------------------------------------------------------
#' @importFrom flowWorkspace gh_pop_get_gate gs_pop_get_gate
#' @importFrom methods is
qc_gates <- function(gs, gate) {
  gt <- gh_pop_get_gate(gs[[1]], gate)

  if (is(gt, "rectangleGate")) {
    n <- length(gt@min)
    if (n == 1) {
      p <- qc_1d_gates(gs, gate)
    } else if (n == 2) {
      p <- qc_2d_gates(gs, gate)
    } else {
      stop("Can't have more than 2 channels")
    }
  } else if (is(gt, "ellipsoidGate") | is(gt, "polygonGate")) {
    p <- qc_polygon_gates(gs, gate)
  } else {
    stop("Gate should be one of `rectagleGate`, `ellipsoidGate` or `polygonGate` object")
  }

  p
}

#' @importFrom flowCore exprs
#' @importFrom ggplot2 ggplot geom_vline xlim labs facet_grid
qc_1d_gates <- function(gs, gate) {
  # retrieve gates
  gates <- gs_pop_get_gate(gs, gate)
  channel <- names(gates[[1]]@min)[1]
  batch <- pData(gs)$batch
  if (is.null(batch)) {
    batch <- NA
  }

  # build data
  mins <- sapply(gates, function(x) x@min)
  dt <- data.frame(
    sample = names(mins),
    i = seq_along(mins),
    channel = mins,
    batch = batch,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  xmax <- max(exprs(gh_pop_get_data(gs[[1]]))[, channel])

  # plot
  title <- sprintf("%s gates by sample", gate)
  p <- ggplot(dt) +
    geom_vline(
      aes(xintercept = channel, group = sample, text = i),
      color = "red",
      alpha = 0.4
    ) +
    xlim(0, xmax) +
    labs(title = title) +
    xlab(channel)

  if (!is.null(batch)) {
    p <- p +
      facet_wrap("batch") +
      labs(title = sprintf("%s (faceted by batch)", title))
  }

  p
}

#' @importFrom methods is as
#' @importFrom flowCore exprs
#' @importFrom ggplot2 geom_polygon ylim theme
qc_polygon_gates <- function(gs, gate) {
  # retrieve gates
  gates <- gs_pop_get_gate(gs, gate)
  batch <- pData(gs)$batch
  if (is.null(batch)) {
    batch <- NA
  }

  outlier <- pData(gs)$outlier

  # build data
  tmp <- lapply(seq_along(gates), function(i) {
    sample <- names(gates)[i]
    gt <- gates[[sample]]
    if (is(gt, "ellipsoidGate")) {
      gt <- as(gt, "polygonGate")
    }
    bd <- gt@boundaries
    df <- data.frame(
      sample = sample,
      i = i,
      bd,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (is.null(outlier)) {
      df$outlier <- FALSE
    } else {
      df$outlier <- outlier[i]
    }
    if (!is.null(batch)) {
      df$batch <- batch[i]
    }

    df
  })
  dt <- do.call(rbind, tmp)
  channels <- colnames(dt)[3:4]
  xmax <- max(exprs(gh_pop_get_data(gs[[1]]))[, channels[1]])
  ymax <- max(exprs(gh_pop_get_data(gs[[1]]))[, channels[2]])

  title <- sprintf("%s gates by sample", gate)
  p <- ggplot(dt, aes_(x = as.name(channels[1]), y = as.name(channels[2]))) +
    geom_polygon(aes(group = sample, text = i, col = outlier), fill = NA, alpha = 0.25) +
    xlim(0, xmax) +
    ylim(0, ymax) +
    labs(title = title)

  if (is.null(outlier)) {
    p <- p + theme(legend.position = "none")
  }

  if (!is.null(batch)) {
    p <- p +
      facet_wrap("batch") +
      labs(title = sprintf("%s (faceted by batch)", title))
  }

  p
}

# plot utility functions -------------------------------------------------------
get_channel <- function(gs, marker) {
  if (is.null(names(marker))) {
    cf <- gh_pop_get_data(gs[[1]], node)
    channel <- getChannelMarker(cf, marker)$name
  } else {
    channel <- names(marker)
  }
}

make_label <- function(gs, channel, node) {
  marker <- unname(markernames(gs)[channel])
  label <- ifelse(is.na(marker), channel, paste(channel, marker))
  paste(label, "@", node)
}

get_range <- function(cs, channel) {
  c(
    min(unlist(lapply(cs, function(x) min(exprs(x[, channel])[, 1])))),
    max(unlist(lapply(cs, function(x) max(exprs(x[, channel])[, 1]))))
  )
}

# plot generating functions ----------------------------------------------------
plot_density_by_marker <- function(gs, marker, node = get_parent(gs), by = "batch") {
  channel <- get_channel(gs, marker)
  plot_density_by_channel(gs, channel, node)
}

plot_ecdf_by_marker <- function(gs, marker, node = get_parent(gs), by = "batch") {
  channel <- get_channel(gs, marker)
  plot_ecdf_by_channel(gs, channel, node)
}

#' @importFrom stats density
#' @importFrom flowCore exprs
#' @importFrom ggplot2 xlab ylab geom_path ggtitle
plot_density_by_channel <- function(gs, channel, node = get_parent(gs), by = "batch") {
  cs <- gs_pop_get_data(gs, node)
  pd <- pData(gs)
  label <- make_label(gs, channel, node)

  densities <- lapply(sampleNames(cs), function(x) {
    tmp <- density(exprs(cs[[x]][, channel]))
    df <- data.frame(
      sample = x,
      x = tmp$x,
      y = tmp$y,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (by %in% colnames(pd)) {
      df$by <- pd[x, by]
    }
    df
  })
  dt <- do.call(rbind, densities)

  p <- ggplot(dt) +
    xlab(label) +
    ylab("density") +
    theme_light()

  alpha <- 0.2
  if (by %in% colnames(pd)) {
    p <- p +
      geom_path(aes(x = x, y = y, group = sample, color = by), alpha = alpha) +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = by, override.aes = list(alpha = 1)))
  } else {
    p <- p +
      geom_path(aes(x = x, y = y, group = sample), alpha = alpha)
  }

  p
}

#' @importFrom stats ecdf
#' @importFrom ggplot2 stat_function aes_
plot_ecdf_by_channel <- function(gs, channel, node = get_parent(gs), by = "batch") {
  cs <- gs_pop_get_data(gs, node)
  pd <- pData(gs)
  label <- make_label(gs, channel, node)

  ecdf_list <- lapply(sampleNames(cs), function(x) {
    ecdf(exprs(cs[[x]][, channel])[, 1])
  })
  channel_range <- data.frame(x = get_range(cs, channel))

  if (by %in% colnames(pd)) {
    mapping <- lapply(pd[, by], function(x) aes_(color = x))
  } else {
    mapping <- lapply(sampleNames(gs), function(x) aes())
  }

  stat_list <- Map(
    f = stat_function,
    fun = ecdf_list, geom = "path", alpha = 0.6, mapping = mapping
  )

  p <- ggplot(channel_range, aes(x = x)) +
    stat_list +
    xlab(label) +
    ylab("F(x)") +
    theme_light()

  if (by %in% colnames(pd)) {
    p <- p +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = by, override.aes = list(alpha = 1)))
  }

  p
}

#' @importFrom stats density
#' @importFrom flowCore exprs
#' @importFrom ggplot2 facet_wrap theme_light theme guides guide_legend
plot_markers <- function(gs) {
  nc <- gs_pop_get_data(gs, get_parent(gs))
  pd <- pData(gs)
  markers <- markernames(gs)

  densities <- lapply(sampleNames(gs), function(x) {
    l <- lapply(names(markers), function(channel) {
      dat <- exprs(nc[[x]])[, channel]
      dens <- density(dat)
      df <- data.frame(
        sample = x,
        marker = unname(markers[channel]),
        x = dens$x,
        y = dens$y,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      if (!is.null(pd$batch)) {
        df$batch <- pd[x, "batch"]
      }
      df
    })
    do.call(rbind, l)
  })
  dt <- do.call(rbind, densities)

  p <- ggplot(dt) +
    xlab("intensity") +
    ylab("density") +
    facet_wrap("marker")

  alpha <- 0.2
  if (!is.null(pd$batch)) {
    p <- p +
      geom_path(aes(x = x, y = y, group = sample, color = batch), alpha = alpha)
  } else {
    p <- p +
      geom_path(aes(x = x, y = y, group = sample), alpha = alpha)
  }

  p +
    theme_light() +
    theme(legend.position = "top") +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
}

#' @importFrom flowWorkspace gh_get_pop_paths
#' @importFrom flowCore exprs
#' @importFrom ggplot2 geom_density stat
plot_density_by_gate <- function(gh, channel) {
  gates <- gh_get_pop_paths(gh, path = 1)

  s <- lapply(gates, function(gate) {
    expr <- exprs(gh_pop_get_data(gh, gate))[, channel]
    tmp <- data.frame(gate = gate, expr = expr)
    colnames(tmp)[2] <- channel
    tmp
  })

  dt <- do.call(rbind, s)
  dt$gate <- factor(dt$gate, gates)

  ggplot(dt, aes_(x = as.name(channel))) +
    geom_density(aes(y = stat(count), color = gate))
}


# qc file saving functions -----------------------------------------------------
#' @importFrom gtools mixedorder
#' @importFrom pheatmap pheatmap
save_spillover_heatmaps <- function(gs, output_dir) {
  catf("Saving spillover matrix heatmaps")

  output_dir <- file.path(output_dir, "spillover")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_names <- sampleNames(gs)
  sapply(sample_names, function(x) {
    filename <- file.path(output_dir, paste0(x, ".png"))
    mat <- get_spillover(gh_pop_get_data(gs[[x]]))
    rownames(mat) <- colnames(mat)
    ord <- mixedorder(colnames(mat))
    mat <- mat[ord, ord]
    p <- pheatmap(
      mat = mat,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      display_numbers = TRUE,
      main = "Spillover Matrix",
      filename = filename
    )

    filename
  })
}

save_gate_plots <- function(gs, output_dir) {
  catf("Saving QC gate plots")

  output_dir <- file.path(output_dir, "gates")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_names <- sampleNames(gs)
  sapply(sample_names, function(x) {
    p <- autoplot(gs[[x]], bins = 64)
    filename <- file.path(output_dir, paste0(x, ".png"))
    suppressWarnings(ggsave(filename, ggcyto_arrange(p)))
    filename
  })
}

#' @importFrom flowWorkspace sampleNames
#' @importFrom ggcyto autoplot ggcyto_arrange
#' @importFrom ggplot2 ggsave
save_gate_plots_by_gate <- function(gs, gate, output_dir) {
  catf(sprintf("Saving QC plot of %s gate", gate))

  output_dir <- file.path(output_dir, "gates", gate)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_names <- sampleNames(gs)
  sapply(sample_names, function(x) {
    p <- autoplot(gs[[x]], gate, bins = 64)
    filename <- file.path(output_dir, paste0(x, ".png"))
    suppressWarnings(ggsave(filename, ggcyto_arrange(p), dpi = 80))
    filename
  })
}

save_marker_plots <- function(gs, output_dir) {
  catf("Saving QC marker plots")

  p <- plot_markers(gs)
  filename <- file.path(output_dir, "markers.png")
  ggsave(filename, p)

  filename
}

save_density_plots_by_marker <- function(gs, marker, output_dir) {
  catf(sprintf("Saving density plot of %s", marker))

  output_dir <- file.path(output_dir, "markers")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  marker_channel <- paste0(marker, "_", names(marker))
  filename <- sprintf("%s/%s.png", output_dir, gsub("/", "_", marker_channel))
  catf(sprintf("Saving QC plot of %s at the terminal node", marker_channel))
  p <- plot_density_by_marker(gs, marker)
  suppressWarnings(ggsave(filename, p, dpi = 80))

  filename
}

save_density_plots_by_gate <- function(gs, marker, output_dir) {
  catf(sprintf("Saving desnity plot of %s by gate", marker))

  marker_channel <- paste0(marker, "_", names(marker))
  output_dir <- file.path(output_dir, "density", marker_channel)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_names <- sampleNames(gs)
  sapply(sample_names, function(x) {
    p <- plot_density_by_gate(gs[[x]], names(marker))
    filename <- file.path(output_dir, paste0(x, ".png"))
    suppressWarnings(ggsave(filename, p, dpi = 80))

    filename
  })
}


# report rendering functions ---------------------------------------------------
#' @importFrom rmarkdown render
render_qc_report <- function(gs_dir) {
  catf("Compiling QC report")
  file_path <- file.path(gs_dir, "QC.html")

  catf(sprintf("output_file: %s", file_path))
  input <- file.path(gs_dir, "QC.Rmd")
  file.copy(system.file("qc/QC.Rmd", package = "HIPCCyto"), input, overwrite = TRUE)
  render(
    input = input,
    output_file = file_path,
    params = list(gs_dir = gs_dir)
  )

  file_path
}

render_study_report <- function(study_dir) {
  catf("Compiling study report")

  input <- file.path(study_dir, "study.Rmd")
  file.copy(system.file("qc/study.Rmd", package = "HIPCCyto"), input, overwrite = TRUE)
  render(
    input = input,
    params = list(study_dir = study_dir)
  )

  file_path <- file.path(study_dir, "study.html")
  catf(sprintf("output_file: %s", file_path))

  file_path
}


# outlier detection functions --------------------------------------------------
find_outliers <- function(gs, byBatch = TRUE, cl = 0.99, step = 3) {
  batch <- unique(pData(gs)$batch)

  if (isFALSE(byBatch) | length(batch) <= 1) {
    test_outliers(gs, cl, step)
  } else {
    unlist(lapply(batch, function(x) {
      catf(x)
      test_outliers(subset(gs, batch == x), cl, step)
    }))
  }
}

#' @importFrom diptest dip.test
#' @importFrom flowCore exprs
#' @importFrom stats cov t.test
test_outliers <- function(gs, cl = 0.99, step = 3) {
  if (length(cl) == 1) {
    cl <- rep(cl, 3)
  }

  # step 1: dip test
  cs <- gs_pop_get_data(gs, "Lymphocytes")
  pv <- lapply(cs, function(x) dip.test(exprs(x)[, "FSC-A"])$p.value)
  outliers <- names(which(pv < (1 - cl[1])))
  catf("step #1: dip test")
  catf(outliers)
  if (step == 1) {
    return(outliers)
  }

  if (sum(!sampleNames(gs) %in% outliers) > 5) {
    # step 2: location test
    gates <- gs_pop_get_gate(gs[!sampleNames(gs) %in% outliers], "Lymphocytes")
    mat <- t(sapply(gates, function(x) {
      c(
        ux = unname(x@mean[1]),
        uy = unname(x@mean[2]),
        sx = x@cov[1, 1],
        sy = x@cov[2, 2],
        sxy = x@cov[1, 2],
        det = abs(det(x@cov))
      )
    }))

    mu <- colMeans(mat[, 1:2])
    sigma <- cov(mat[, 1:2])
    converted <- convert_points_to_ellipsoid(mat[, 1:2], mu, sigma)
    inside <- check_inside_ellipse(converted, p = cl[2])
    catf("step #2: location test")
    catf(names(which(!inside)))
    outliers <- c(outliers, names(which(!inside)))
    if (step == 2) {
      return(outliers)
    }

    # step 3: size test
    size <- mat[inside, "det"]
    res <- t.test(size, conf.level = cl[3])
    inside <- size > res$conf.int[1] & size < res$conf.int[2]
    catf("step #3: size test")
    catf(names(which(!inside)))
    outliers <- c(outliers, names(which(!inside)))
  } else {
    catf("insufficient samples for steps #2 and #3")
  }

  outliers
}

# ellipsoid functions adopted from SIBER ---------------------------------------
# https://github.com/AndrewLJackson/SIBER/blob/master/R/pointsToEllipsoid.R
# https://github.com/AndrewLJackson/SIBER/blob/master/R/ellipsoidTransform.R
convert_points_to_ellipsoid <- function(points, mu, sigma) {
  if (ncol(sigma) != nrow(sigma)) {
    stop("`sigma` must be a square matrix")
  }
  if (ncol(points) != ncol(sigma)) {
    stop("number of columns in `points` must be of same dimension as `sigma`")
  }
  if (length(mu) != ncol(sigma)) {
    stop("length of `mu` must be of same dimension as `sigma`")
  }

  # inverse of sigma
  eig <- eigen(sigma)
  sigma_inverse <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

  # transform the points
  t(apply(points, 1, function(point) solve(sigma_inverse, point - mu)))
}

# https://github.com/AndrewLJackson/SIBER/blob/master/R/ellipseInOut.R
#' @importFrom stats qchisq
check_inside_ellipse <- function(points, p = 0.95, radius = NULL) {
  # if radius is NULL, calculate it based on p
  if (is.null(radius)) {
    radius <- qchisq(p, df = ncol(points))
  }

  # determine if each point is inside this radius or not
  rowSums(points ^ 2) < radius
}
