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
