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
