#' @importFrom flowWorkspace save_gs
#' @importFrom jsonlite toJSON
#' @importFrom ggcyto autoplot ggcyto
#' @importFrom ggplot2 ggsave facet_wrap facet_null aes_ aes
#' @importFrom ggridges geom_density_ridges
#' @importFrom flowCore getChannelMarker
#' @export
save_gating_sets <- function(gsl, output_dir) {
  lapply(seq_len(length(gsl)), function(i) {
    gs <- gsl[[i]]

    path <- file.path(output_dir, paste0("gs", i))
    message(sprintf(">> Saving gating set #%s...", i))
    message(path)
    save_gs(gs, path, overwrite = TRUE, cdf = "copy")

    message(">> Saving gating set summary...")
    pd <- pData(gs)
    summary <- list(
      number_of_samples = length(gs),
      number_of_markers = length(markernames(gs)),
      markers = mixedsort(markernames(gs)),
      nodes = gs_get_pop_paths(gs, path = 1),
      sample_types = unique(pd$type),
      timepoints = unique(paste(pd$study_time_collected, pd$study_time_collected_unit)),
      cohorts = unique(pd$cohort),
      batches = as.list(table(pd$batch)),
      participants = unique(pd$participant_id)
    )
    json <- toJSON(summary, pretty = TRUE, auto_unbox = TRUE)
    cat(json, file = sprintf("%s/summary", path))

    if (!dir.exists(file.path(path, "gates"))) {
      dir.create(file.path(path, "gates"))
    }
    if (!dir.exists(file.path(path, "markers"))) {
      dir.create(file.path(path, "markers"))
    }

    for (node in gs_get_pop_paths(gs, path = 1)[-1]) {
      message(sprintf(">> Saving QC plot of %s gate...", node))
      p <- autoplot(gs, node)
      ggsave(
        filename = sprintf("%s/gates/%s.pdf", path, node),
        plot = p,
        width = 50,
        height = 50,
        dpi = 320,
        limitsize = FALSE
      )
    }

    nc <- gs_pop_get_data(gs, get_parent(gs))
    for (channel in colnames2(gs)) {
      p <- ggcyto(nc, aes_(x = as.name(channel))) +
        geom_density_ridges(aes(y = name))

      if (!is.null(pData(gs)$batch)) {
        p <- p + facet_wrap(~ batch)
      } else {
        p <- p + facet_null()
      }

      marker <- paste(rev(getChannelMarker(nc[[1]], channel)[1, ]), collapse = "_")
      message(sprintf(">> Saving QC plot of %s at the terminal node...", marker))
      ggsave(
        filename = sprintf("%s/markers/%s.pdf", path, gsub("/", "_", marker)),
        plot = p,
        width = 50,
        height = 50,
        dpi = 320,
        limitsize = FALSE
      )
    }

    path
  })
}
