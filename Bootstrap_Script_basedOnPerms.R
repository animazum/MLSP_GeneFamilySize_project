###################################
#### PERMS PERMS PERMS PERMS ######
###################################

PGLS.GF.size.perms<-function(traits.cl, pheno.col.nums, cores,tree.cl, spp.col.num, gene.numbers.cl, where.Save.it, No.perms, percentTotalSample){
  log_print("###########################", hide_notes = T)
  log_print("##### Perms perms PGLS ####", hide_notes = T)
  log_print("###########################", hide_notes = T)
  log_print("calc. cores", hide_notes = T)
  nucleos<-detectCores()
  nucleos<-nucleos + cores - nucleos
  cl <- makeCluster(nucleos) 
  log_print(paste(nucleos + cores - nucleos," cores detected and ready to roll"), hide_notes = T)
  
  # bootstrap dataset creation
  # size of the data frame 
  # Sizesample<-as.numeric(round(dim(traits.filtered)[1]*(99/100)))
  Sizesample<-as.numeric(round(dim(traits.cl)[1]*(percentTotalSample/100)))
  # Resample species with replacement
  # resampled_species <- sample(traits.filtered[,1], replace = F, size = 55)
  resampled_species <- sample(traits.cl[, spp.col.num], replace = F, size = Sizesample)

  # traits.cl <- traits.filtered[traits.filtered[, 1] %in% resampled_species, , drop = FALSE]
  traits.cl <- traits.cl[traits.cl[, spp.col.num] %in% resampled_species, , drop = FALSE]
  
  # Drop tips from the tree that are not in the resampled dataset
  # tree.cl <- drop.tip(tree, setdiff(tree$tip.label, traits.filtered[, 1]))
  tree.cl <- drop.tip(tree.cl, setdiff(tree.cl$tip.label, traits.cl[, 1]))
  
  # Drop species in gene numbers
  # gene.numbers.cl <- gene.numbers.filtered[,colnames(gene.numbers.filtered) %in% resampled_species , drop = FALSE]
  gene.numbers.cl <- gene.numbers.cl[,colnames(gene.numbers.cl) %in% resampled_species , drop = FALSE]
  
  #gene.numbers.cl<- gene.numbers.cl
  #traits.cl<- traits.cl
  #tree.cl<- tree.cl
  # phenotypes<<- colnames(x = traits.cl)[c(1, 10)]
  phenotypes<<- colnames(x = traits.cl)[c(spp.col.num, pheno.col.nums)]
  # traits.cl<<-traits.cl
  #No.perms <- No.perms
  
  ##### this way  clusterExports stops using .GlobalEnv and uses the eviropnment inside the function. 
  clusterExport(cl = cl , varlist = c("gene.numbers.cl", "traits.cl", "gls", "tree.cl", "corBrownian", "phenotypes", "No.perms"), envir = environment())
  
  log_print("Starting perms PGLS", hide_notes = T)
  
  if (length(pheno.col.nums) == 2) {
    perms<<- parLapply(cl, 1:No.perms, function(i,...){
      traits.cl<-traits.cl
      y<-sample(c(1:length(gene.numbers.cl)))
      gene_numbers2 <- gene.numbers.cl[,y]
      outperm1<-apply(gene_numbers2, MARGIN = 1,function(x){
        traits2<<-cbind(traits.cl[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c("Spp.Name" , "Av.Body.mass" , "SSD")],x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c(cat(paste0(dQuote(phenotypes[1:length(phenotypes)], F), collapse=" , ")))],x)) #add phenotypes as needed
        # The dQuote adds double quotes, the paste0 inserts the commas and cat shows the result without escaping special characters.
        #list2env(traits2, envir = .GlobalEnv)## Assign them to the global environment
        names(traits2)[dim(traits2)[2]]<-"GFS"
        
        # formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
        # formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        # formilin2<<-as.formula(paste("~",colnames(traits.cl)[spp.col.num], sep = ""))
        formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
        
        pglsModel<-try(gls(formilin, correlation = corBrownian(phy = tree.cl, form= formilin2), method = "ML", data = traits2))  #add phenotypes as needed
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], 
                 anova(pglsModel)$`p-value`[2], anova(pglsModel)$`F-value`[3], anova(pglsModel)$`p-value`[3], 
                 summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[3,1], 
                 summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], summary(pglsModel)$tTable[3,2], 
                 summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[3,3], 
                 summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4], summary(pglsModel)$tTable[3,4]))
      })
      outperm1<-as.data.frame(t(outperm1))
      colnames(outperm1)<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "), paste("p_value", phenotypes[2], sep=" "), 
                            paste("F_value", phenotypes[3], sep=" "), paste("p_value", phenotypes[3], sep=" "), "Coefficient Intercept", 
                            paste("Coefficient", phenotypes[2], sep=" "), paste("Coefficient", phenotypes[3], sep=" "),  "Coef SE Intercept", 
                            paste("Coef SE", phenotypes[2], sep=" "), paste("Coef SE", phenotypes[3], sep=" "), "Coef tval Intercept", 
                            paste("Coef tval", phenotypes[2], sep=" "), paste("Coef tval", phenotypes[3], sep=" "), "Coef pval Intercept",
                            paste("Coef pval", phenotypes[2], sep=" "), paste("Coef pval", phenotypes[3], sep=" "))
      
      pglsList<- list(outperm1)
      return(outperm1)
      
    })
  } else if (length(pheno.col.nums) == 1) {
    
    perms<<- parLapply(cl, 1:No.perms, function(i,...){
      traits.cl<-traits.cl
      y<-sample(c(1:length(gene.numbers.cl)))
      gene_numbers2 <- gene.numbers.cl[,y]
      outperm1<-apply(gene_numbers2, MARGIN = 1,function(x){
        traits2<<-cbind(traits.cl[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c("Spp.Name" , "SSD" , "Av.Body.mass")],x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c(cat(paste0(dQuote(phenotypes[1:length(phenotypes)], F), collapse=" , ")))],x)) #add phenotypes as needed
        # The dQuote adds double quotes, the paste0 inserts the commas and cat shows the result without escaping special characters.
        #list2env(traits2, envir = .GlobalEnv)## Assign them to the global environment
        names(traits2)[dim(traits2)[2]]<-"GFS"
        
        formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        # formilin2<<-as.formula(paste("~",colnames(traits.cl)[spp.col.num], sep = ""))
        formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
        
        
        # pglsModel<-pgls(formilin, correlation = corBrownian(1, phy = tree.cl, form = formilin2), method = "ML", data = traits2) #add phenotypes as needed
        # pglsModel<-gls(formilin, correlation = corBrownian(1, phy = tree.cl, form = formilin2), method = "ML", data = traits2, ) #add phenotypes as needed
        
        pglsModel<-try(gls(formilin, correlation = corBrownian(phy = tree.cl, form= formilin2), method = "ML", data = traits2))  #add phenotypes as needed
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2],
                 anova(pglsModel)$`p-value`[2], summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1],
                 summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], summary(pglsModel)$tTable[1,3],
                 summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4]))
      })
      outperm1<-as.data.frame(t(outperm1))
      colnames(outperm1)<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "),
                            paste("p_value", phenotypes[2], sep=" "), "Coefficient Intercept",
                            paste("Coefficient", phenotypes[2], sep=" "), "Coef SE Intercept",
                            paste("Coef SE", phenotypes[2], sep=" "), "Coef tval Intercept",
                            paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept",
                            paste("Coef pval", phenotypes[2], sep=" "))
      
      pglsList<- list(outperm1)
      return(outperm1)
      
    })
  }
  
  log_print("Perms PGLS done", hide_notes = T)
  
  perms<<- as.data.frame(do.call(rbind, perms), rownames=TRUE)
  log_print(colnames(perms))
  
  for (i in 1:length(pheno.col.nums)) {
    
    log_print("creating variables for stats", hide_notes = T)
    
    #traits2<-as.data.frame(cbind(traits.cl[,c(phenotypes[1:length(phenotypes)])],x)) #add phenotypes as needed
    
    # perms<- bind_rows(perms,.id = NULL)
    sps<-nrow(traits.cl)
    DeFr<- nrow(traits.cl) - (1 + (length(phenotypes)-1)) ## the -1 removes the col.name Spp.Name that its not a phenotype.
    
    head(traits.cl)
    #### the 62 are the numbers of spp minus 2 = degrees of freedom
    log_print("tval and R values", hide_notes = T)
    
    Coef.tval1<<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    
    head(perms)
    colnames(perms[[1]])
    rownames(perms)
    R_t_value <<- ((perms[,Coef.tval1]))/sqrt((((perms[,Coef.tval1])^2) + DeFr )) # square root of (tvalue^2 / (tvalue^2 + DF))
    log_print(head(R_t_value), hide_notes = T)
    
    R1 <- paste("R t value", phenotypes[i+1], sep=" ")
    perms[R1]<<- R_t_value
    log_print(" benjamini correction", hide_notes = T)
    
    ## benjamini correction
    perms[paste("p.adjusted Benjamini GFS vs ", phenotypes[i+1], sep = "")] <<- p.adjust(perms[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
    #permsRs<<-permsRs
  }
  log_print("Creating output file", hide_notes = T)
  nameout<- paste("PGLS_perms", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits.cl)[1], "Spp", dim(perms)[1], "Perm_GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  write.csv(perms, paste(where.Save.it, nameout, sep = "/"))
  perms<<- perms
  
  closeAllConnections()
}

# works with 2 phenos only at the moment
PGLS.GF.size.perms(traits.cl = traits.filtered, cores = 5, pheno.col.nums = c(10), tree.cl = tree, spp.col.num = 1, 
                   gene.numbers.cl = gene.numbers.filtered, where.Save.it = dirSave, No.perms = 10, percentTotalSample = 99)
colnames(perms)
