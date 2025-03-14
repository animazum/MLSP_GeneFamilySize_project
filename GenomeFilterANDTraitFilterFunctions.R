

##### This script has the Genome.completion.filtering and Trait.orthoGr.filtering functions

###############################################
#** The Genome.completion.filtering function filters and *
#** processes gene family data to identify high-quality *
#** single-copy genes and assess genome completeness for *
#** phylogenetic analysis. *
###############################################


Genome.completion.filtering <- function(FamData, Filt.Near.Uni, ssp.treshold, IDnameGC) {
  # Input validation
  if (!is.data.frame(FamData)) stop("FamData must be a data frame")
  if (!is.numeric(Filt.Near.Uni) || Filt.Near.Uni < 0 || Filt.Near.Uni > 100) 
    stop("Filt.Near.Uni must be between 0 and 100")
  if (!is.numeric(ssp.treshold) || ssp.treshold < 0 || ssp.treshold > 100) 
    stop("ssp.treshold must be between 0 and 100")
  
  want.percent <- Filt.Near.Uni/100
  NumCols1 <- length(FamData)
  
  # Filter ortho counts for single-copy orthologs
  list1 <- FamData[1:NumCols1] %>% filter_all(all_vars(.< 2))
  list1$Total <- rowSums(list1)
  
  # Calculate nearly universal genes
  filtering.90.1 <- filter(list1, Total >= round(NumCols1 * want.percent))
  
  # Calculate threshold for species
  treshold <- round(nrow(filtering.90.1) * (ssp.treshold/100))
  
  log_print(paste(round(nrow(filtering.90.1))))
  log_print(paste("All species must have at least", treshold, "genes"))
  
  # Filter species based on threshold
  filtering.90.1 <- filtering.90.1[, colSums(filtering.90.1) >= treshold]
  
  removed_species <- setdiff(colnames(FamData), colnames(filtering.90.1))
  log_print(paste("Removed species:", paste(removed_species, collapse = ", ")))
  
  # Create final filtered dataset
  gene_numbers <- FamData[colnames(FamData) %in% colnames(filtering.90.1)]
  
  # Output files
  output_file1 <- paste0("nearly.universal.", Filt.Near.Uni, ".", IDnameGC, 
                         ".percent.", ssp.treshold, "ssp.treshold.csv")
  output_file2 <- paste0("Total.counts.", IDnameGC, ".", ncol(gene_numbers), "spp.csv")
  
  write.csv(filtering.90.1, output_file1)
  write.csv(gene_numbers, output_file2)
  
  log_print(paste("Created files:", output_file1, "and", output_file2))
  
  return(list(
    gene_numbers = gene_numbers,
    filtering_results = filtering.90.1
  ))
}

##################################################
#** The Trait.orthoGr.filtering function filters and aligns trait and gene number datasets *
#** for downstream comparative analyses. It removes species with missing traits, ensures *
#** species alignment between datasets, and filters out gene families with excessive missing values, *
#** low variance, or presence in only a single species.*
###################################################

Trait.orthoGr.filtering <- function(traits, pheno.col.nums, spp.col.num, gene.numbers, IDNameTOF) {
  # Input validation
  if (!is.data.frame(traits)) stop("traits must be a data frame")
  if (!is.numeric(pheno.col.nums)) stop("pheno.col.nums must be numeric")
  if (!is.numeric(spp.col.num)) stop("spp.col.num must be numeric")
  
  phenotypes <- colnames(traits)[pheno.col.nums]
  
  # Filter traits and align datasets
  traits.filtered <- traits
  for (phenotype in phenotypes) {
    traits.filtered <- traits.filtered[!is.na(traits.filtered[, phenotype]), ]
  }
  
  spp.name <- colnames(traits)[spp.col.num]
  traits.filtered <- traits.filtered[traits.filtered[, spp.name] %in% colnames(gene.numbers), ]
  gene.numbers.filtered <- gene.numbers[, colnames(gene.numbers) %in% traits.filtered[, spp.name]]
  
  # Log missing species information
  missing_in_traits <- setdiff(colnames(gene.numbers), traits.filtered[, spp.name])
  missing_in_genes <- setdiff(traits[[spp.name]], traits.filtered[[spp.name]])
  
  if (length(missing_in_traits) > 0) {
    log_print("Species missing in traits:", paste(missing_in_traits, collapse = ", "))
  }
  if (length(missing_in_genes) > 0) {
    log_print("Species missing in gene numbers:", paste(missing_in_genes, collapse = ", "))
  }
  
  # Filter gene families
  zeros_threshold <- ncol(gene.numbers.filtered) * 0.2
  gene.numbers.filtered <- gene.numbers.filtered[rowSums(gene.numbers.filtered == 0) <= zeros_threshold, ]
  gene.numbers.filtered <- gene.numbers.filtered[apply(gene.numbers.filtered, 1, var) > 0, ]
  gene.numbers.filtered <- gene.numbers.filtered[apply(gene.numbers.filtered, 1, max) > 2, ]
  
  # Save results
  output_file <- paste0("Filtered.GeneNum.", IDNameTOF, ".", paste(phenotypes, collapse = "_"), ".csv")
  write.csv(gene.numbers.filtered, output_file)
  log_print(paste("Created file:", output_file))
  
  return(list(
    traits_filtered = traits.filtered,
    gene_numbers_filtered = gene.numbers.filtered
  ))
}


# genome_results <- Genome.completion.filtering(
#   FamData = gene.numbers.filtered,
#   Filt.Near.Uni = 90,
#   ssp.treshold = 90,
#   IDnameGC = "genomeResults"
# )
# # 
# trait_results <- Trait.orthoGr.filtering(
#   traits = traits.filtered,
#   pheno.col.nums = c(12,13),        # Adjust based on your data
#   spp.col.num = 2,                # Adjust based on your data
#   gene.numbers = genome_results$gene_numbers,
#   IDNameTOF = "traitsResults"
# )

