Systematic_LOO_Analysis <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, 
                                    No.variables, where.Save.it, Filt.Near.Uni = 90, 
                                    ssp.treshold = 90, n_bootstrap = NULL) {
  
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
  
  # Create directory if it doesn't exist
  dir.create(where.Save.it, showWarnings = FALSE, recursive = TRUE)
  
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
      log_print(paste("Too few species left after removing", removed_spp, "- Skipping this iteration."))
      next
    }
    
    # Add genome completion filtering
    tryCatch({
      # Generate a unique ID for this iteration
      iter_id <- paste0("iter_", i, "_", removed_spp)
      
      # Run genome completion filtering with updated return values
      genome_results <- Genome.completion.filtering(
        FamData = gene_numbers_sub,
        Filt.Near.Uni = Filt.Near.Uni,
        ssp.treshold = ssp.treshold,
        IDnameGC = iter_id
      )
      
      # Run trait orthogroup filtering with updated return values
      trait_results <- Trait.orthoGr.filtering(
        traits = traits_sub,
        pheno.col.nums = pheno.col.nums,
        spp.col.num = spp.col.num,
        gene.numbers = genome_results$gene_numbers,
        IDNameTOF = iter_id
      )
      
      # Use the returned filtered results
      gene_numbers_pgls <- trait_results$gene_numbers_filtered
      traits_pgls <- trait_results$traits_filtered
      
      # Ensure tree matches filtered data
      tree_pgls <- drop.tip(tree_sub, setdiff(tree_sub$tip.label, traits_pgls$Spp.Name))
      
      ## sort gene numbers by traits.filtered
      gene_numbers_pgls <- gene_numbers_pgls[, match(traits_pgls$Spp.Name, colnames(gene_numbers_pgls))]
      
      # Try running PGLS
      model_success <- FALSE
      error_msg <- NA
      
      if (nrow(traits_pgls) >= 5) {  # Recheck sample size after filtering
        PGLS.GF.WithID(
          traits = traits_pgls,
          pheno.col.nums = pheno.col.nums,
          tree = tree_pgls,
          spp.col.num = spp.col.num,
          gene.numbers = gene_numbers_pgls,
          No.variables = No.variables,
          where.Save.it = where.Save.it, 
          IDForfile = paste(i, removed_spp,sep = ".")
        )
        model_success <- TRUE
      } else {
        error_msg <- "Too few species after filtering"
      }
      
    }, error = function(e) {
      error_msg <<- conditionMessage(e)
      model_success <<- FALSE
      log_print(paste("Error in iteration", i, ":", error_msg))
    })
    
    # Store results
    loo_results <- rbind(loo_results, data.frame(
      Iteration = i,
      Removed_Species = removed_spp,
      Success = model_success,
      Error_Message = if(is.na(error_msg)) "Success" else error_msg,
      stringsAsFactors = FALSE
    ))
    
    # Save intermediate results after each iteration
    write.csv(loo_results, 
              file = paste0(where.Save.it, "/Systematic_LOO_PGLS_Results_interim.csv"), 
              row.names = FALSE)
  }
  
  # Save final results
  saveRDS(loo_results, file = paste0(where.Save.it, "/Systematic_LOO_PGLS_Results.rds"))
  write.csv(loo_results, 
            file = paste0(where.Save.it, "/Systematic_LOO_PGLS_Results_final.csv"), 
            row.names = FALSE)
  
  log_print("#### Systematic Leave-One-Out Analysis Completed ####", hide_notes = T)
  return(loo_results)
}


results <- Systematic_LOO_Analysis(
  traits = traits.filtered,
  pheno.col.nums = c(12,13),
  tree = tree,
  spp.col.num = 2,
  gene.numbers = gene.numbers.filtered,
  No.variables = 2,
  where.Save.it = "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_LOO_sys_errorhandding/PGLS_github/",
  Filt.Near.Uni = 90,               # Custom genome completion threshold (95%)
  ssp.treshold = 90,                # Custom species threshold (85%)
  #n_bootstrap = 3                   # Limit to first 10 species
)
