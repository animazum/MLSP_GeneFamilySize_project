Systematic_LOO_Analysis <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it, n_bootstrap = NULL) {
  phenotypes <- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  log_print("#########################", hide_notes = T)
  log_print("##### Systematic LOO PGLS ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  
  # Get all species names
  all_species <- traits$Spp.Name
  
  # If n_bootstrap is not provided, use all species
  if (is.null(n_bootstrap)) {
    n_bootstrap <- length(all_species)
  } else {
    # If n_bootstrap is provided but larger than number of species, cap it
    n_bootstrap <- min(n_bootstrap, length(all_species))
  }
  
  # Storage for results
  loo_results <- data.frame(
    Iteration = integer(),
    Removed_Species = character(),
    Success = logical(),
    Error_Message = character(),
    stringsAsFactors = FALSE
  )
  
  # Systematically go through each species
  for (i in 1:n_bootstrap) {
    # Take species sequentially instead of randomly
    removed_spp <- all_species[i]
    cat("\nRunning LOO iteration", i, "Removing", removed_spp, "\n")
    
    # Subset data
    traits_sub <- traits[traits$Spp.Name != removed_spp, ]
    gene_numbers_sub <- gene.numbers[, colnames(gene.numbers) %in% traits_sub$Spp.Name]
    tree_sub <- drop.tip(tree, removed_spp)
    
    # Ensure dataset consistency
    tree_sub <- drop.tip(tree_sub, setdiff(tree_sub$tip.label, traits_sub$Spp.Name))
    
    if (nrow(traits_sub) < 5) {
      cat("Too few species left after removing", removed_spp, "- Skipping this iteration.\n")
      next
    }
    
    # Try running PGLS
    model_success <- FALSE
    error_msg <- NA
    tryCatch({
      PGLS.GF.WithID(
        traits = traits_sub,
        pheno.col.nums = pheno.col.nums,
        tree = tree_sub,
        spp.col.num = spp.col.num,
        gene.numbers = gene_numbers_sub,
        No.variables = No.variables,
        where.Save.it = where.Save.it, 
        IDForfile = paste(removed_spp, i, sep = ".")
      )
      model_success <- TRUE
    }, error = function(e) {
      error_msg <<- conditionMessage(e)
    })
    
    # Store results
    loo_results <- rbind(loo_results, data.frame(
      Iteration = i,
      Removed_Species = removed_spp,
      Success = model_success,
      Error_Message = error_msg,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save results
  saveRDS(loo_results, file = paste(where.Save.it, "Systematic_LOO_PGLS_Results.rds", sep = "/"))
  write.csv(loo_results, paste(where.Save.it, "Systematic_LOO_PGLS_Results.csv", sep = "/"), row.names = FALSE)
  
  log_print("#### Systematic Leave-One-Out Analysis Completed ####", hide_notes = T)
  return(loo_results)
}

# Usage example:
pathSystematic<-"/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_LOO_systematic/"

results_Systematic_LOO <- Systematic_LOO_PGLS(
  traits = traits.filtered, 
  pheno.col.nums = c(12,13),
  tree = tree, 
  spp.col.num = 2, 
  gene.numbers = gene.numbers.filtered, 
  No.variables = 2, 
  where.Save.it = pathSystematic
  # Note: n_bootstrap is not specified, so it will run for all species
)
