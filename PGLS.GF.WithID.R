
PGLS.GF.WithID<- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it, IDForfile){
  phenotypes<<- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  print("#########################", hide_notes = T)
  print("##### PGLS PGLS PGLS ####", hide_notes = T)
  print("#########################", hide_notes = T)
  if (No.variables == 1){
    ###################################
    print("#### PGLS ONE Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,1,FUN = function(x){
      
      traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      tree<- drop.tip(tree,names(x)[is.na(x)]) 
      x
      ## creating the formula for the model for 1 variable
      formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin2<-as.formula(paste("~",phenotypes[1], sep = ""))
      
      pglsModel<<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(pglsModel, "try-error")) 
        return(c(NA,NA))
      else
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2], 
                 summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], 
                 summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4]))
    })))
    colnames(out)<<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "), paste("p_value", phenotypes[2], sep=" "), 
                      "Coefficient Intercept", paste("Coefficient", phenotypes[2], sep=" "), "Coef SE Intercept", paste("Coef SE", phenotypes[2], sep=" "), 
                      "Coef tval Intercept", paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept", paste("Coef pval", phenotypes[2], sep=" "))
    
  } 
  else if (No.variables == 2) {
    ###################################
    print("#### PGLS TWO Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,MARGIN = 1,FUN = function(x){
      pheno1<<-phenotypes[3]
      pheno2<<-phenotypes[2]
      traits2<<-cbind(traits[,c("Spp.Name", pheno2, pheno1)],as.numeric(x)) #add phenotypes as needed
      # traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      tree<- drop.tip(tree,names(x)[is.na(x)])
      ## creating the formula for the model for 2 variables
      # formilin<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
      #print((paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ ")))
      # formilin2<-as.formula(paste("~",colnames(traits)[spp.col.num], sep = ""))
      formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
      #print((paste("~",colnames(traits)[spp.col.num], sep = "")))
      # pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      pglsModel<<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(pglsModel, "try-error"))
        return(c(NA,NA))
      else
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2],
                 anova(pglsModel)$`F-value`[3], anova(pglsModel)$`p-value`[3], summary(pglsModel)$tTable[1,1],
                 summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[3,1], summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2],
                 summary(pglsModel)$tTable[3,2], summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[3,3],
                 summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4], summary(pglsModel)$tTable[3,4]))
    })))
    colnames(out)<<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[3] , sep=" "), paste("p_value", phenotypes[3], sep=" "),
                      paste("F_value", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), "Coefficient Intercept", paste("Coefficient", phenotypes[3], sep=" "),
                      paste("Coefficient", phenotypes[2], sep=" "),  "Coef SE Intercept", paste("Coef SE", phenotypes[3], sep=" "), paste("Coef SE", phenotypes[2], sep=" "),
                      "Coef tval Intercept", paste("Coef tval", phenotypes[3], sep=" "), paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept",
                      paste("Coef pval", phenotypes[3], sep=" "), paste("Coef pval", phenotypes[2], sep=" "))
  }
  else if (No.variables > 2){
    ###################################
    ## PGLS more than two Variables ###
    ###################################
    print("At the moment the function only works with upto 2 variables, sorry :P", hide_notes = T)
  }
  for (i in 1:No.variables) {
    ###########################
    ### R^2 calculation ###
    ###########################
    print("Calc. Rs", hide_notes = T)
    sps<-nrow(traits2)
    DeFr <- sps - (1 + No.variables)
    Coef.tval1<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    R1 <- toString(paste("R t value", phenotypes[i+1], sep=" "))
    R_t_value1 <- out[[Coef.tval1]] / (sqrt(((out[[Coef.tval1]])^2) + DeFr)) # tvalue / square root of(tvalue^2 + DF))
    out[[R1]] <<- R_t_value1
    ### benjamini correction
    print("Calc. benjamini correction", hide_notes = T)
    out[[paste("p.adjusted GFS vs.", phenotypes[i+1], sep = "")]] <<- p.adjust(out[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
  }
  nameout<- paste("PGLS_LOO", IDForfile, paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits2)[1], "Spp", dim(out)[1], "GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  print(paste("object created:", paste( "out", sep = "")), hide_notes = T)
  write.csv(out, paste(where.Save.it, nameout, sep = "/"))
}

# PGLS.GF.WithID(traits = trait_results$traits_filtered, 
#                pheno.col.nums = c(12,13), 
#                No.variables = 2, 
#              tree = tree, 
#              spp.col.num = 2, 
#              gene.numbers = trait_results$gene_numbers_filtered, 
#              where.Save.it = "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_LOO_sys_errorhandding/PGLS/", 
#              IDForfile = "000000hola")
