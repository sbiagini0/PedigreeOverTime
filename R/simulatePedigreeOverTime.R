# ──────────────────────────────────────────────────────────────────────────────
# Main function: simulatePedigreeOverTime()
# -----------------------------------------------------------------------------
# Purpose
#   Drive the end-to-end workflow for each pedigree:
#     • Add typed members chronologically.  
#     • Generate intermediate pedigree versions.  
#     • Run EP/IP simulations at every step.  
#     • Produce and save high-resolution plots and RData objects.
#
# Features
#   • Parallel back-end handled internally via {future}/{furrr}; if `ncores`
#     is not supplied, the function auto-detects *total_cores − 1*.  
#   • Re-simulation logic: calls `pedigrees_to_run()` to skip pedigrees that
#     are already up-to-date.  
#   • Saves for each family (fam_id):  
#       − `*_all_ped.rdata` (all incremental pedigrees)  
#       − `*_EP.rdata` - Expected Probability simulations  
#       − `*_IP.rdata` - Identification Probability simulations  
#       − `<fam_id>_pedigree.jpeg` (diagram)  
#       − `<fam_id>_simulation.jpeg` (EP/IP panels)  
#       − `<fam_id>.jpeg` (diagram + panels side-by-side).  
#   • Restores the user’s original future plan on exit, even if an error occurs.
#
# Returns
#   Invisibly returns NULL. All tangible results are written to disk inside the
#   folders defined by `path_info`.
# -----------------------------------------------------------------------------
simulatePedigreeOverTime <- function(df, 
                                     pedigree, 
                                     threshold = 10000, 
                                     nsims = nsims, 
                                     seed = 2000, 
                                     ncores = NULL,
                                     path = path_info
                                     ) {
  # --- Check pedigrees to simulate ---  
  run_info <- pedigrees_to_run(df, path_info)
  
  if (length(run_info$to_run) == 0) {
    message("All pedigrees are up to date — no new simulations to run.")
    return(invisible(NULL))
  }
  
  if (length(run_info$new_peds) > 0)
    message("New pedigrees to simulate: ",
            paste(run_info$new_peds, collapse = ", "))
  
  if (length(run_info$updated_peds) > 0)
    message("Pedigrees to re-simulate (new members): ",
            paste(run_info$updated_peds, collapse = ", "))
  
  # Keep only the rows we actually need
  df <- dplyr::filter(df, Pedigree %in% run_info$to_run)
  
  # --- Set up parallel plan ---
  if (is.null(ncores))
    ncores <- max(1, parallel::detectCores() - 1)
  
  old_plan <- future::plan()
  future::plan(multisession, workers = ncores)
  on.exit(future::plan(old_plan), add = TRUE)
  
  # --- Set seed ---
  set.seed(seed)
  
  # --- df selection ---
  if (!exists("df")) stop("The dataframe 'df' is not defined in the global environment.")
  
  filtered_df <- df
  
  unique_peds <- unique(filtered_df$Pedigree)
  fam_seeds <- setNames(seed + seq_along(unique_peds) - 1, unique_peds)
  
  for (fam_id in unique(filtered_df$Pedigree)) {
    cat("Proccesing:", fam_id, "\n")
    
    members_df <- filtered_df %>%
      dplyr::filter(Pedigree == fam_id) %>%
      arrange(`Typed_date`)
    
    members_to_add <- members_df$`Member_id`
    full_ped <- pedigree[[fam_id]]
    if (is.null(full_ped)) {
      warning(paste("No pedigree found for ", fam_id))
      next
    }
    
    # --- Create base ped ---
    base_ped <- setAlleles(full_ped, ids = members_to_add, alleles = 0)
    
    # --- Appende peds ---
    all_ped <- list(base_ped)
    for (i in seq_along(members_to_add)) {
      partial_ids <- members_to_add[1:i]
      ped_i <- transferMarkers(full_ped, base_ped, ids = partial_ids, erase = FALSE)
      all_ped[[i + 1]] <- ped_i
    }
    
    labels <- c("Basal", purrr::map_chr(seq_along(members_to_add), 
                                           ~ paste(members_to_add[1:.x], collapse = ", ")))
    
    # --- Save ped rdata ---
    save(all_ped, file = paste0(path_info$ped_path, fam_id, "_all_ped.rdata"))
    
    # --- Plot colors ---
    member_palette <- c("lightgreen", "firebrick1", "deepskyblue", "#FFFF33", "gray70", "#F781BF", "wheat", "cyan", "#FF7F00")
    n_mem <- length(members_to_add)
    members_colors <- member_palette[seq_len(n_mem)]
    pedigree_plot_path_full <- plot_pedigree(base_ped, 
                                             full_ped, 
                                             fam_id, 
                                             members_to_add, 
                                             members_colors)
    
    # --- Simulation ---
    pedigree_seed <- fam_seeds[[fam_id]]
    
    EP <- future_map(all_ped, function(x) {
      tryCatch({
        missingPersonEP(x, missing = "Missing person", verbose = FALSE)
      }, error = function(e) {
        warning("Error in missingPersonEP for ", fam_id, "\n", e)
        return(NULL)
      })
    }, .options = furrr_options(seed = pedigree_seed)) %>% set_names(labels)
    
    IP <- future_map2(all_ped, seq_along(all_ped),
                      function(x, seed_val) {
                        tryCatch({
                          missingPersonIP(x,
                                          missing          = "Missing person",
                                          nsim             = nsims,
                                          threshold        = threshold,
                                          disableMutations = FALSE,
                                          seed             = pedigree_seed,
                                          verbose          = FALSE)
                        }, error = function(e) {
                          warning("Error in missingPersonIP for ", fam_id, "\n", e)
                          return(NULL)
                        })
                      }, .options = furrr_options(seed = pedigree_seed)) %>%
      set_names(labels)
    
    # --- Save rdata ---
    save(EP, file = paste0(path_info$EP_path, fam_id, "_EP.rdata"))
    save(IP, file = paste0(path_info$IP_path, fam_id, "_IP.rdata"))
    
    # --- Simulation results ---
    simulation_plot_path_full <- plot_simulation_results(EP, 
                                                         IP, 
                                                         labels, 
                                                         members_colors, 
                                                         fam_id, 
                                                         pedigree_seed, 
                                                         nsims)
    
    # --- Combine final images ---
    tryCatch({
      ped_img <- image_read(pedigree_plot_path_full)
      sim_img <- image_read(simulation_plot_path_full)
      final_img <- image_append(c(ped_img, sim_img), stack = FALSE)
      image_write(final_img, paste0(path_info$output_path, "/", fam_id, ".jpeg"))
    }, error = function(e) {
      warning("Error combining images for ", fam_id, "\n", e)
    })
  }
}