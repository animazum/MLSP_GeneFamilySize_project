# Diagnostic function
check_data_patterns <- function(traits2, pheno1, pheno2) {
  # Check for duplicate patterns
  duplicates <- duplicated(traits2[, c(pheno2, pheno1)])
  if(any(duplicates)) {
    cat("Duplicate patterns found in rows:", which(duplicates), "\n")
  }
  
  # Check for extreme values
  for(col in c(pheno2, pheno1)) {
    outliers <- boxplot.stats(traits2[[col]])$out
    if(length(outliers) > 0) {
      cat("Outliers found in", col, ":", outliers, "\n")
    }
  }
  
  # Check correlations
  cor_matrix <- cor(traits2[, c(pheno2, pheno1)])
  if(any(abs(cor_matrix[upper.tri(cor_matrix)]) > 0.9)) {
    cat("High correlations detected between variables\n")
    print(cor_matrix)
  }
}

PGLS.GF.WithID <- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it, IDForfile) {
  phenotypes <- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  print("#########################", hide_notes = T)
  print("##### PGLS PGLS PGLS ####", hide_notes = T)
  print("#########################", hide_notes = T)
  
  gene_families <- rownames(gene.numbers)
  n_species <- nrow(traits)
  
  if (No.variables == 2) {
    # Create results matrix directly instead of using apply
    results <- matrix(NA, nrow = nrow(gene.numbers), ncol = 18)
    
    for(i in 1:nrow(gene.numbers)) {
      x <- gene.numbers[i,]
      gene_name <- gene_families[i]
      
      pheno1 <- phenotypes[3]
      pheno2 <- phenotypes[2]
      traits2 <- cbind(traits[,c("Spp.Name", pheno2, pheno1)], as.numeric(x))
      names(traits2)[dim(traits2)[2]] <- "GFS"
      
      var_gfs <- var(traits2$GFS)
      if(var_gfs < 1e-10) {
        cat("Gene Family:", gene_name, "\n")
        cat("GFS values:", toString(unique(traits2$GFS)), "\n")
        cat("GFS variance:", var_gfs, "\n\n")
      }
      
      # Remove species with NAs
      tree <- drop.tip(tree, names(x)[is.na(x)])
      
      # Create formulas
      formilin <- as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
      formilin2 <- as.formula(paste("~", phenotypes[1], sep = ""))
      
      # Try running the model with additional checks
      model_results <- tryCatch({
        # Scale predictors
        traits2_scaled <- traits2
        for(col in c(pheno2, pheno1, "GFS")) {
          var_val <- var(traits2[[col]])
          if(var_val < 1e-10) {
            warning(paste("Very low variance in", col, ":", var_val))
            return(rep(NA, 18))
          }
          traits2_scaled[[col]] <- scale(traits2[[col]], center = TRUE, scale = TRUE)
        }
        
        # Check correlations
        cors <- cor(traits2_scaled[,c(pheno2, pheno1, "GFS")])
        if(any(abs(cors[upper.tri(cors)]) > 0.99)) {
          warning("Perfect correlation detected between variables")
          return(rep(NA, 18))
        }
        
        # Fit model
        pglsModel <- gls(model = formilin,
                         correlation = corBrownian(phy = tree, form = formilin2),
                         method = "ML",
                         data = traits2_scaled,
                         control = glsControl(tolerance = 1e-6,
                                              maxIter = 1000,
                                              opt = "nlminb"))
        
        # Extract statistics
        f_values <- anova(pglsModel)$`F-value`
        p_values <- anova(pglsModel)$`p-value`
        t_table <- summary(pglsModel)$tTable
        
        c(f_values[1], p_values[1], 
          f_values[2], p_values[2],
          f_values[3], p_values[3], 
          t_table[1,1], t_table[2,1], t_table[3,1],
          t_table[1,2], t_table[2,2], t_table[3,2],
          t_table[1,3], t_table[2,3], t_table[3,3],
          t_table[1,4], t_table[2,4], t_table[3,4])
        
      }, error = function(e) {
        warning(paste("Model fitting error:", conditionMessage(e)))
        return(rep(NA, 18))
      })
      
      results[i,] <- model_results
    }
    
    # Convert to data frame
    out <- as.data.frame(results)
    rownames(out) <- gene_families
    
    # Set column names
    colnames(out) <- c("F_value Intercept", "p_value Intercept", 
                       paste("F_value", phenotypes[3], sep=" "), paste("p_value", phenotypes[3], sep=" "),
                       paste("F_value", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), 
                       "Coefficient Intercept", paste("Coefficient", phenotypes[3], sep=" "),
                       paste("Coefficient", phenotypes[2], sep=" "), "Coef SE Intercept", 
                       paste("Coef SE", phenotypes[3], sep=" "), paste("Coef SE", phenotypes[2], sep=" "),
                       "Coef tval Intercept", paste("Coef tval", phenotypes[3], sep=" "), 
                       paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept",
                       paste("Coef pval", phenotypes[3], sep=" "), paste("Coef pval", phenotypes[2], sep=" "))
    
    # Calculate R²
    for (i in 1:No.variables) {
      print("Calc. Rs", hide_notes = T)
      DeFr <- n_species - (1 + No.variables)
      Coef.tval1 <- toString(paste("Coef tval", phenotypes[i+1], sep=" "))
      R1 <- toString(paste("R t value", phenotypes[i+1], sep=" "))
      
      out[[R1]] <- sapply(1:nrow(out), function(row) {
        tval <- as.numeric(out[row, Coef.tval1])
        if(!is.na(tval)) {
          return(tval / sqrt(tval^2 + DeFr))
        } else {
          return(NA)
        }
      })
    }
    
    # Calculate benjamini correction
    for(i in 1:No.variables) {
      p_val_col <- paste("p_value", phenotypes[i+1], sep=" ")
      adj_p_col <- paste("p.adjusted GFS vs.", phenotypes[i+1], sep="")
      out[[adj_p_col]] <- p.adjust(out[[p_val_col]], method = "fdr")
    }
    
  } else {
    stop("At the moment the function only works with 2 variables")
  }
  
  # Save results
  nameout <- paste("PGLS_LOO", IDForfile, 
                   paste(phenotypes[2:length(phenotypes)], collapse = "_"),
                   n_species, "Spp", dim(out)[1], "GeneFams",
                   "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  
  print(paste("Location and name of your file:", 
              paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  print(paste("object created:", "out"), hide_notes = T)
  
  write.csv(out, paste(where.Save.it, nameout, sep = "/"))
  
  return(out)
}


sppout<-c("Bos_taurus")
traits.filtered2out <- traits.filtered[!traits.filtered$Spp.Name %in% sppout, ]
tree2out<- drop.tip(tree, tip = sppout)
gene.numbers2out<- gene.numbers.filtered[!colnames(gene.numbers.filtered) %in% sppout]

pathSystematic<-"/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_LOO_sys_errorhandding/"

PGLS.GF.WithID(
  traits = traits.filtered,
  pheno.col.nums = c(12,13),
  tree = tree,
  spp.col.num = 2,
  gene.numbers = gene.numbers.filtered,
  No.variables = 2,
  where.Save.it = pathSystematic, 
  IDForfile = "errorHandling"
)
