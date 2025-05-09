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
