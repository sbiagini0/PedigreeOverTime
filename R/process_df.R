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
