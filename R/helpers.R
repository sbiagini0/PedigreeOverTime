# ──────────────────────────────────────────────────────────────────────────────
# Helper: setup_directories()
# -----------------------------------------------------------------------------
# Purpose
#   Create (if necessary) a standard output-folder hierarchy and return a list
#   (`path_info`) with fully-qualified paths used throughout the pipeline.
#
# Features
#   • Accepts a single argument `base_path` (default `"Output"`).  
#   • Automatically creates the base directory and five sub-directories:  
#       • `Pedigree plot/` – individual pedigree JPEGs  
#       • `Simulation plot/` – panel JPEGs with EP/IP results  
#       • `EP/` – RData files containing EP objects  
#       • `IP/` – RData files containing IP objects  
#       • `RData peds/` – saved `*_all_ped.rdata` files  
#   • Uses `recursive = TRUE`, so the entire tree is built even if parent
#     folders don’t exist.  
#   • Returns every path with a trailing slash for easy concatenation.
#
# Returns
#   A named list of character strings:  
#     $base, $pedigree_plot_path, $simulation_plot_path,  
#     $EP_path, $IP_path, $ped_path, and $output_path.
# -----------------------------------------------------------------------------
setup_directories <- function(base_path = "output") {
  
  # ensure base folder exists
  if (!dir.exists(base_path)) dir.create(base_path, recursive = TRUE)
  
  # sub-folders to create
  subdirs <- c(
    pedigree_plot_path   = "Pedigree plot",
    simulation_plot_path = "Simulation plot",
    EP_path              = "EP",
    IP_path              = "IP",
    ped_path             = "RData peds"
  )
  
  # build full paths and create them if necessary
  full_paths <- file.path(base_path, subdirs)
  invisible(lapply(full_paths, function(p) if (!dir.exists(p)) dir.create(p)))
  
  # return list with trailing slashes for convenience
  c(
    list(base = paste0(base_path, "/")),
    setNames(paste0(full_paths, "/"), names(subdirs)),
    output_path = paste0(base_path, "/")
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: read_famfile()
# -----------------------------------------------------------------------------
# Purpose
#   Read a multi-family `.fam` file in DVI (Disaster Victim Identification) mode
#   and return the pedigrees as a named list.
#
# Features
#   • Calls `pedFamilias::readFam()` with `useDVI = TRUE`, ensuring that all
#     family IDs inside the file are preserved.  
#   • Adds the prefix `"EXTRA"` to any individuals introduced during DVI
#     reconstruction (`prefixAdded = "EXTRA"`).  
#   • Automatically simplifies redundant pedigree structures (`simplify1 = TRUE`)
#     and removes duplicate entries (`deduplicate = TRUE`).  
#   • Wraps the call in `tryCatch()`—if reading fails, issues a warning and
#     returns `NULL` instead of stopping the entire pipeline.
#
# Returns
#   A named list of pedigree objects keyed by their Group Family ID, or
#   `NULL` if the file cannot be read.
# -----------------------------------------------------------------------------
read_famfile <- function(file) {
  tryCatch({
    pedFamilias::readFam(
      file,
      useDVI     = TRUE,
      prefixAdded = "EXTRA",
      simplify1  = TRUE,
      deduplicate = TRUE,
      verbose    = FALSE
    )
  }, error = function(e) {
    warning("Error reading the file: ", file, "\n", e)
    NULL
  })
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: process_pedigree()
# -----------------------------------------------------------------------------
# Purpose
#   Clean the list returned by `read_famfile()`:
#   • Extract the core `_comp1` pedigree object from each MPI entry.
#   • Remove or report NULL pedigrees.
#   • Set `famid()` for every pedigree to match its list name.
#
# Features
#   • Prints a message listing any NULL pedigrees it encounters.  
#   • Optional `drop_null` argument (default TRUE) lets you keep or discard
#     NULL entries.  
#   • Uses `purrr::map()` for concise extraction of `_comp1`.
#
# Returns
#   A cleaned named list of pedigree objects—each with `famid()` already set—
#   ready to be passed to `simulatePedigreeOverTime()`.
# -----------------------------------------------------------------------------
process_pedigree <- function(mpi, drop_null = TRUE) {
  
  # 1) Extract the main `_comp1` component from every entry
  mpi <- purrr::map(mpi, ~ .x[["_comp1"]])
  
  # 2) Identify NULL pedigrees and notify the user
  null_peds <- names(mpi)[vapply(mpi, is.null, logical(1))]
  if (length(null_peds) > 0) {
    message("NULL pedigrees detected (will be skipped): ",
            paste(null_peds, collapse = ", "))
    if (drop_null) mpi <- mpi[!names(mpi) %in% null_peds]
  }
  
  # 3) Assign each pedigree its family ID (famid) = list name
  for (i in seq_along(mpi))
    famid(mpi[[i]]) <- names(mpi)[i]
  
  mpi
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: process_df()
# -----------------------------------------------------------------------------
# Purpose : Locate, read, and clean the timeline Excel file that lists when
#           each pedigree member was genotyped.
#
# Features
#   • Auto-detects the first `.xlsx` file in the specified `data_dir`
#     (defaults to `"data"`).  
#   • Verifies that the required columns — Pedigree, Member_id, Typed_date —
#     are present; stops execution if any are missing.  
#   • Trims column names, converts `Typed_date` to `Date`, removes rows with
#     missing dates, and sorts chronologically within each pedigree.
#
# Returns : A clean `data.frame` ready for downstream analysis.
# -----------------------------------------------------------------------------
process_df <- function(data_dir = "data") {
  
  # 1) Locate an Excel file
  xlsx_files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
  
  if (length(xlsx_files) == 0)
    stop("No .xlsx file found in the ", data_dir, " folder.")
  if (length(xlsx_files) > 1)
    warning("Multiple .xlsx files detected; using the first one: ",
            basename(xlsx_files[1]))
  
  excel_path <- xlsx_files[1]
  
  # 2) Read and validate
  df <- readxl::read_excel(excel_path) |>
    dplyr::rename_with(trimws)
  
  required_cols <- c("Pedigree", "Member_id", "Typed_date")
  if (!all(required_cols %in% names(df)))
    stop("Excel must contain columns: ", paste(required_cols, collapse = ", "))
  
  # 3) Clean and order
  df |>
    dplyr::select(all_of(required_cols)) |>
    dplyr::mutate(Typed_date = as.Date(Typed_date)) |>
    dplyr::filter(!is.na(Typed_date)) |>
    dplyr::arrange(Pedigree, Typed_date)
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: pedigrees_to_run()
# -----------------------------------------------------------------------------
# Purpose
#   Analyse the current timeline (`df`) against previously saved results
#   to decide which pedigrees should be simulated, re-simulated, or skipped.
#
# Features
#   • Counts members currently typed (`n_now`) for every pedigree in the Excel.  
#   • Reads each `<fam_id>_all_ped.rdata` (if present) to recover the number of
#     members simulated last time (`n_prev`).  
#   • Flags pedigrees whose final JPEG is missing, even if the member count
#     is unchanged.  
#   • Classifies pedigrees into:
#       − **new_peds**      : never simulated before.  
#       − **updated_peds**  : previously simulated but now outdated (extra
#         members or missing image).  
#       − **to_run**        : union of the two, ready for processing.
#
# Returns
#   A list with three character vectors:
#     $new_peds, $updated_peds, $to_run
#   These are used by `simulatePedigreeOverTime()` to control workflow.
# -----------------------------------------------------------------------------
pedigrees_to_run <- function(df, path_info) {
  
  ## a) Members present right now
  members_now <- df |>
    dplyr::group_by(Pedigree) |>
    dplyr::summarise(n_now = dplyr::n(), .groups = "drop")
  
  ## b) Members simulated last time
  get_prev_count <- function(fam) {
    rdata_path <- file.path(path_info$ped_path, paste0(fam, "_all_ped.rdata"))
    if (!file.exists(rdata_path)) return(0L)
    e <- new.env()
    load(rdata_path, envir = e)
    as.integer(length(e$all_ped) - 1L)
  }
  
  members_prev <- data.frame(
    Pedigree = members_now$Pedigree,
    n_prev   = vapply(members_now$Pedigree, get_prev_count, integer(1))
  )
  
  ## c) Does the final JPEG exist?
  jpeg_missing <- function(fam) {
    pat <- paste0("^", fam, ".*\\.jpe?g$")
    !any(grepl(
      pat,
      list.files(path_info$output_path, pattern = "\\.jpe?g$", full.names = FALSE)
    ))
  }
  missing_img <- vapply(members_now$Pedigree, jpeg_missing, logical(1))
  
  ## d) Combine criteria
  out <- dplyr::inner_join(members_now, members_prev, by = "Pedigree") |>
    dplyr::mutate(
      new_ped      = n_prev == 0L,                       # never simulated
      needs_update = (n_now > n_prev) | missing_img      # outdated or missing
    )
  
  new_peds     <- out$Pedigree[out$new_ped]
  updated_peds <- out$Pedigree[!out$new_ped & out$needs_update]
  to_run       <- c(new_peds, updated_peds)
  
  list(
    new_peds     = new_peds,
    updated_peds = updated_peds,
    to_run       = to_run
  )
}

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

