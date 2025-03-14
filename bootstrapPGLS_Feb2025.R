bootstrap_PGLS <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, 
                           No.variables, where.Save.it, n_bootstrap = 100) {
  
  clusterExport(cl = cl , varlist = c("gene.numbers", "traits", "gls", "tree", "corBrownian", "phenotypes", "n_bootstrap"), envir = environment())
  
  trait <- "log10MLSP_RelBrainMass"
  tmp <- file.path(getwd(), paste("Bootstrap-GenefamilySize.",trait,".associated.genes.",Correcting,".correction.log", sep = ""))
  lf <- log_open(tmp, traceback = F)
  results_list <- list()
  
  for (i in 1:n_bootstrap) {
    cat("Running bootstrap iteration", i, "\n")
    
    # Randomly remove 5% of species
    drop_spp <- sample(traits$Spp.Name, size = 1)
    print(drop_spp)
    traits_boot <- traits[!traits$Spp.Name %in% drop_spp, ]
    gene_numbers_boot <- gene.numbers[, colnames(gene.numbers) %in% traits_boot$Spp.Name]
    tree_boot <- drop.tip(tree, drop_spp)
    
    # Run PGLS
    result <- try(PGLS.GF.size(traits = traits_boot, pheno.col.nums = pheno.col.nums, tree = tree_boot,
                               spp.col.num = spp.col.num, gene.numbers = gene_numbers_boot, 
                               No.variables = No.variables, where.Save.it = where.Save.it))
    
    if (!inherits(result, "try-error")) {
      results_list[[i]] <- result
    }
  }
  
  return(results_list)
}

bootstrap_results <- bootstrap_PGLS(
  traits = traits.filtered,
  pheno.col.nums = c(12,13),
  tree = tree,
  spp.col.num = 2,
  gene.numbers = gene.numbers.filtered,
  No.variables = 2,
  where.Save.it = path1,
  n_bootstrap = 10
)
