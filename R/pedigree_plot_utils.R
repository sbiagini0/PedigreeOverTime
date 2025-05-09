# ──────────────────────────────────────────────────────────────────────────────
# Helper: plot_pedigree()
# -----------------------------------------------------------------------------
# Purpose
#   Create a pedigree diagram for a given family (`fam_id`) at its current
#   genotyping state, shading typed members and annotating each symbol with
#   the first genotyping date found in the global dataframe `df`.
#
# Features
#   • Colours / fills the symbols of typed individuals (`member_selected`)
#     using the palette provided in `members_colors`.  
#   • Extracts genotyping dates from `df` and writes them underneath each
#     individual.  
#   • Relies on `pedtools::plot()` for consistent pedigree layout.  
#   • Saves the diagram as `<fam_id>_pedigree.jpeg` inside
#     `path_info$pedigree_plot_path`.  
#   • Performs a sanity check: stops with an informative error if `df` or
#     `path_info` are missing from the global environment.
#
# Returns
#   A single character string: the full path to the saved JPEG file.
# -----------------------------------------------------------------------------
plot_pedigree <- function(ped_base, ped_completo, fam_id, member_selected, members_colors) {
  if (!exists("df")) stop("The dataframe 'df' is not defined in the global environment")
  if (!exists("path_info")) stop("path_info must be defined before calling plot_pedigree()")
  
  ped_labels <- labels(ped_base)
  annotation_vector <- sapply(ped_labels, function(id) {
    date <- df %>% 
      dplyr::filter(grepl(fam_id, Pedigree) & `Member_id` == id) %>% 
      dplyr::pull(`Typed_date`)
    if (length(date) == 0) "" else as.character(date[1])
  })
  
  out_file <- paste0(path_info$pedigree_plot_path, fam_id, "_pedigree.jpeg")
  jpeg(out_file, width = 2100, height = 2160, res = 300)
  plot(
    ped_base,
    cex   = 0.8,
    hatched = typedMembers,
    margin  = c(5),
    branch  = 1,
    packed  = FALSE,
    align   = FALSE,
    col     = setNames(as.list(member_selected), members_colors),
    carrier    = member_selected,
    textAnnot = list(bottom = list(annotation_vector, cex = 0.5, col = 1, offset = 1.5)),
    keep.par = TRUE,
    arrows   = FALSE
  )
  title(fam_id, cex.main = 1.2)
  dev.off()
  out_file
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: plot_simulation_results()
# -----------------------------------------------------------------------------
# Purpose
#   Visualise EP/IP simulation outputs for a pedigree by building two
#   powerPlot panels— EP vs IP and Expected Mismatch vs log-LR—and
#   combining them into a single quality JPEG.
#
# Features
#   • Accepts pre-computed `EP` and `IP` lists plus a `labels` vector so each
#     simulation stage is clearly identified.  
#   • Colours points with `members_colors`, adds in-plot text labels showing
#     exact coordinates, and suppresses the default legend to reduce clutter.  
#   • Aligns both ggplot-based grobs to a common width, then stacks them under
#     a title grob that states the number of simulations (`nsims`) and
#     seed (`pedigree_seed`).  
#   • Saves the final panel as `<fam_id>_simulation.jpeg` in
#     `path_info$simulation_plot_path`.  
#   • Validates that `path_info` exists in the global environment; stops with
#     an informative error if not found.
#
# Returns
#   A character string giving the full path to the saved JPEG file.
# -----------------------------------------------------------------------------
plot_simulation_results <- function(EP, IP, labels, members_colors, fam_id,
                                    pedigree_seed, nsims) {
  if (!exists("path_info")) stop("path_info must be defined before calling plot_simulation_results()")
  
  power_colors <- c("black", members_colors)
  
  plot4 <- powerPlot(ep = EP, ip = IP, type = 1, labs = labels,
                     col = power_colors, majorpoints = TRUE, minorpoints = TRUE)
  plot5 <- powerPlot(ep = EP, ip = IP, type = 3, labs = labels,
                     col = power_colors, majorpoints = TRUE, minorpoints = TRUE)
  
  data_points_p4 <- data.frame(
    x = purrr::map_dbl(EP, ~ .x$EPtotal[1]),
    y = purrr::map_dbl(IP, ~ .x$IP[1]),
    label = labels
  )
  data_points_p5 <- data.frame(
    x = purrr::map_dbl(EP, ~ .x$expectedMismatch[1]),
    y = purrr::map_dbl(IP, ~ .x$meanLogLR[1]),
    label = labels
  )
  
  plot4 <- plot4 + ggplot2::geom_text(
    data = data_points_p4,
    aes(x = x, y = y, label = paste0("(", round(x, 2), ", ", round(y, 2), ")")),
    vjust = -0.5, size = 3, color = "black"
  ) + ggplot2::theme(legend.position = "none")
  
  plot5 <- plot5 + ggplot2::geom_text(
    data = data_points_p5,
    aes(x = x, y = y, label = paste0("(", round(x, 2), ", ", round(y, 2), ")")),
    vjust = -0.5, size = 3, color = "black"
  )
  
  # Align panels
  g1 <- ggplotGrob(plot4)
  g2 <- ggplotGrob(plot5)
  maxWidth <- grid::unit.pmax(g1$widths, g2$widths)
  g1$widths <- maxWidth; g2$widths <- maxWidth
  
  title_grob <- textGrob(
    paste("Number of Simulations:", nsims, "| Seed:", pedigree_seed),
    gp = gpar(fontsize = 10, fontface = "bold"),
    just = "left", x = 0.1
  )
  
  out_file <- paste0(path_info$simulation_plot_path, fam_id, "_simulation.jpeg")
  jpeg(out_file, width = 3840, height = 2160, res = 300)
  grid.arrange(title_grob, arrangeGrob(plot4, plot5, ncol = 2), heights = c(0.05, 1))
  dev.off()
  
  out_file
}
