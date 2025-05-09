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
