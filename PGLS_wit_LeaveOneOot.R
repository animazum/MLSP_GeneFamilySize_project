Bootstrap_LOO_PGLS_Parallel <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it, n_bootstrap = 100) {
  phenotypes <- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  log_print("#########################", hide_notes = T)
  log_print("##### Parallel Bootstrap LOO PGLS ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  
  # Function to run a single iteration of the bootstrap
  run_iteration <- function(iteration, traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it) {
    set.seed(iteration) # Ensure reproducibility for each iteration
    removed_spp <- sample(traits$Spp.Name, 1)
    cat("\nRunning bootstrap iteration", iteration, "Removing", removed_spp, "\n")
    
    # Subset data
    traits_sub <- traits[traits$Spp.Name != removed_spp, ]
    gene_numbers_sub <- gene.numbers[, colnames(gene.numbers) %in% traits_sub$Spp.Name]
    tree_sub <- drop.tip(tree, removed_spp)
    
    # Ensure dataset consistency
    tree_sub <- drop.tip(tree_sub, setdiff(tree_sub$tip.label, traits_sub$Spp.Name))
    
    if (nrow(traits_sub) < 5) {
      return(data.frame(
        Iteration = iteration,
        Removed_Species = removed_spp,
        Success = FALSE,
        Error_Message = "Too few species left after removal",
        stringsAsFactors = FALSE
      ))
    }
    
    # Try running PGLS
    model_success <- FALSE
    error_msg <- NA
    
    tryCatch({
      PGLS.GF.size(
        traits = traits_sub,
        pheno.col.nums = pheno.col.nums,
        tree = tree_sub,
        spp.col.num = spp.col.num,
        gene.numbers = gene_numbers_sub,
        No.variables = No.variables,
        where.Save.it = where.Save.it
      )
      model_success <- TRUE
    }, error = function(e) {
      error_msg <<- conditionMessage(e)
    })
    
    # Return results for this iteration
    return(data.frame(
      Iteration = iteration,
      Removed_Species = removed_spp,
      Success = model_success,
      Error_Message = ifelse(is.na(error_msg), NA, error_msg),
      stringsAsFactors = FALSE
    ))
  }
  
  # Number of cores for parallel processing
  num_cores <- detectCores() - 1 # Leave one core free
  cl <- makeCluster(num_cores) # Create cluster
  
  # Ensure cluster is stopped even if an error occurs
  on.exit(stopCluster(cl), add = TRUE)
  
  # Export necessary variables and libraries to the cluster
  clusterExport(cl, c("traits", "pheno.col.nums", "tree", "spp.col.num", "gene.numbers", "No.variables", "where.Save.it", "PGLS.GF.size"), envir = environment())
  clusterEvalQ(cl, library(ape)) # Ensure required libraries are loaded on workers
  
  # Run the bootstrap iterations in parallel
  loo_results <- parLapply(
    cl, 
    1:n_bootstrap, 
    function(iteration) {
      run_iteration(iteration, traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it)
    }
  )
  
  # Combine results
  loo_results <- do.call(rbind, loo_results)
  
  # Save results
  saveRDS(loo_results, file = paste(where.Save.it, "Bootstrap_LOO_PGLS_Results_Parallel.rds", sep = "/"))
  write.csv(loo_results, paste(where.Save.it, "Bootstrap_LOO_PGLS_Results_Parallel.csv", sep = "/"), row.names = FALSE)
  
  log_print("#### Parallel Bootstrap Leave-One-Out Analysis Completed ####", hide_notes = T)
  return(loo_results)
}