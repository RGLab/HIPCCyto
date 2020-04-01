#' @importFrom flowWorkspace save_gs
#' @importFrom jsonlite toJSON
#' @importFrom ggcyto autoplot ggcyto
#' @importFrom ggplot2 ggsave facet_wrap facet_null aes_ aes
#' @importFrom ggridges geom_density_ridges
#' @export
save_gating_sets <- function(gsl, output_dir) {
  lapply(seq_len(length(gsl)), function(i) {
    gs <- gsl[[i]]

    path <- file.path(output_dir, paste0("gs", i))
    message(sprintf(">> Saving gating set #%s...", i))
    message(path)
    save_gs(gs, path, overwrite = TRUE, cdf = "copy")

    message(">> Saving gating set summary...")
    summary <- list(
      number_of_samples = length(gs),
      number_of_markers = length(markernames(gs)),
      markers = mixedsort(markernames(gs)),
      nodes = gs_get_pop_paths(gs, path = 1),
      batches = as.list(table(pData(gs)$batch)),
      participants = unique(pData(gs)$participant_id)
    )
    json <- toJSON(summary, pretty = TRUE, auto_unbox = TRUE)
    cat(json, file = sprintf("%s/summary", path))

    message(">> Saving QC plots...")
    if (!dir.exists(file.path(path, "gates"))) {
      dir.create(file.path(path, "gates"))
    }
    if (!dir.exists(file.path(path, "markers"))) {
      dir.create(file.path(path, "markers"))
    }

    for (node in gs_get_pop_paths(gs, path = 1)[-1]) {
      p <- autoplot(gs, node)
      ggsave(
        filename = sprintf("%s/gates/%s.pdf", path, node),
        plot =  p,
        width = 50,
        height = 50,
        dpi = 320,
        limitsize = FALSE
      )
    }

    nc <- gs_pop_get_data(gs, get_parent(gs))
    for (marker in markernames(gs)) {
      p <- ggcyto(nc, aes_(x = as.name(marker))) +
        geom_density_ridges(aes(y = name))

      if (!is.null(pData(gs)$batch)) {
        p <- p + facet_wrap(~ batch)
      } else {
        p <- p + facet_null()
      }

      ggsave(
        filename = sprintf("%s/markers/%s.pdf", path, gsub("/", "_", marker)),
        plot =p,
        width = 50,
        height = 50,
        dpi = 320,
        limitsize = FALSE
      )
    }

    path
  })
}
