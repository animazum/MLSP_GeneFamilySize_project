
############################################
######## PGLS     PGLS      PGLS ###########
############################################

# testgf<-head(x = gene.numbers.filtered, n = dim(gene.numbers.filtered)[2])  ### test gene number file with 100 Orthogroups.

# write.tree(tree2, paste("/Users/mona/Dropbox/SSD/mammalian.phylo.",tree2$Nnode,".ssp.SSDchapter.tre",sep = ""))
################
##### PGLS #####
### function ###
################

## for pre pgls analyses
library("logr")
library("dplyr")
library("ggpubr")

## for pgls
library("ape")
library("nlme")
library("parallel")

## for the plot 
library("reshape")
library("ggplot2")

print("change always the ulimit of the computer $ulimit -s 21000 !!!!!!!!!!!!")
rm(list=ls())

path1<- "~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/"

SSDbias <- ""
trait<- "Github.log10MLSP_RelBrainMass"
Correcting<- "benjamini"

tmp <- file.path(getwd(), paste("GenefamilySize.",SSDbias,"longevity.",trait,".associated.genes.",Correcting,".correction.log", sep = ""))
lf <- log_open(tmp, traceback = F)

Nombre.file<-c(paste("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/Datos/phenotype_data_for_46_species.csv", sep = ""))
traits.filtered<-read.csv(Nombre.file, stringsAsFactors = F)
traits.filtered$log10MLSP<- log10(traits.filtered$Maximum.longevity..yrs.)

tree <-read.tree(paste("~/Dropbox/kilili_shared_project/GFSandMLSP_92_Mammals-First_PhD_Project/R Analysis Folders/05-Correlogram/species_tree_edited.nwk", sep = ""))

# Keep only the tips that are in traits.filtered$Genus_species
tips_to_remove <- setdiff(tree$tip.label, traits.filtered$Genus_species)
tree <- drop.tip(phy = tree, tip = tips_to_remove)

#*************   IT IS IMPORTANT TO HAVE THE SAME SPP ORDER BETWEEN TRAITS AND GENE NUMBERSS!!!!!!!   ********###
gene.numbers.filtered <- read.csv(paste("~/Dropbox/kilili_shared_project/GFSandMLSP_92_Mammals-First_PhD_Project/R Analysis Folders/02-Ensembl_GFS/Gene_family_size_per_specie_92spp_r99.csv", sep = ""),header = T, row.names = 1)
dim(traits.filtered)
gene.numbers.filtered<- gene.numbers.filtered[colnames(gene.numbers.filtered) %in% traits.filtered$Ensembl_Species]

# remove families with zeros
FamsWithOs <- gene.numbers.filtered[rowSums(gene.numbers.filtered) < 1, ]
gene.numbers.filtered<-gene.numbers.filtered[!rownames(gene.numbers.filtered) %in% rownames(FamsWithOs),]

# Create a mapping vector with old names as keys and new names as values
name_mapping <- setNames(traits.filtered$Genus_species, traits.filtered$Ensembl_Species)

# Replace column names only where they match, keeping unmatched names unchanged
colnames(gene.numbers.filtered) <- ifelse(colnames(gene.numbers.filtered) %in% traits.filtered$Ensembl_Species, name_mapping[colnames(gene.numbers.filtered)], colnames(gene.numbers.filtered))

colnames(traits.filtered)[2] <- "Spp.Name"
# total of gene counts per spp
sumgenes<-rowSums(t(gene.numbers.filtered))

shapiro.test((traits.filtered$log10MLSP))
shapiro.test((traits.filtered$cEi))

# Remove Homo_sapiens column
# gene.numbers.filtered <- gene.numbers.filtered[, !colnames(gene.numbers.filtered) %in% "Homo_sapiens"]
# traits.filtered <- traits.filtered[ !traits.filtered$Spp.Name %in% "Homo_sapiens",]
# tree<- drop.tip(tree, "Homo_sapiens")


###############################################
#** The Genome.completion.filtering function filters and *
#** processes gene family data to identify high-quality *
#** single-copy genes and assess genome completeness for *
#** phylogenetic analysis. *
###############################################

Genome.completion.filtering<-function(FamData,Filt.Near.Uni, ssp.treshold){
  want.percent<- Filt.Near.Uni/100
  #percetout<<-Filt.Near.Uni
  #ssp.treshold.out<<- ssp.treshold
  ### Single-copy genes, REMOVE!!! the last column that counts genes per orthogroup ###
  NumCols1<-length(FamData)
  
  # Filter the ortho counts by the single-copy orthologs. We are going to remove the orthogroups that have 0 or 1 number of genes per family 
  list1<- FamData[1:NumCols1] %>% filter_all(all_vars(.< 2)) #  all orthogroups have 0 or 1 gene
  list1$Total<- rowSums(list1) # summatory of all the sigle-copy genes in new col. 
  ## Nearly universal genes calculation (filtering to get our core sigle-copy gene set)
  # NumCols1+1 (is the column with the total gene counts per orthogroup)
  
  # Why do I need to filter the orthogroups to keep the ones that have as many genes as the 90% or 95%  of the number of species?
  filtering.90.1<<- filter(list1, Total >= round(NumCols1*want.percent))    
  filtering.90.1<<- filter(list1, Total >= round(NumCols1*want.percent))    
  #filtering.90.1$Total<- NULL 
  # Notes:
  ## So we have 1339 orthogroups with at least a count of 117 genes (90% of the number of the spp).
  
  # keep the spp that have 90% of these single-copy genes orthogroups.
  # Find the treshold for the spp
  # All the spp most have at least the number of genes designated by the treshold. 
  treshold <- round(nrow(filtering.90.1)*(ssp.treshold/100))
  
  log_print(paste(round(nrow(filtering.90.1))))
  log_print(paste("All the spp most have at least the number of genes designated by the treshold:", treshold, "genes", sep = " "), hide_notes = T)
  
  # Count the number of single-copy genes per spp and remove those that does not have at least the number of genes dictated by the treshold.
  filtering.90.1<<- filtering.90.1[,colSums(filtering.90.1) >= treshold]
  filtering.90.1<- filtering.90.1[,colSums(filtering.90.1) >= treshold]
  # Note: apparently the 130 spp have at least 1205 single-copy genes. 
  
  log_print(paste( "Which spp or columns were removed? ", setdiff(colnames(FamData),colnames(filtering.90.1))), hide_notes = T)
  
  ### Jackpot! the final table that will be used to run the PGLS with all the spp that we need. 
  gene_numbers<<-FamData[colnames(FamData) %in% colnames(filtering.90.1)]
  
  log_print(paste("You have", NumCols1, "spp that passed the genome completness filter", sep = " "),hide_notes = T )
  log_print("#### Jackpot! Now we have finished the filtering to know if the genomes used are complete ####", hide_notes = T)
  log_print("The dataframe gene_numbers was created", hide_notes = T)
  log_print("The filtering.90.1 object was created", hide_notes = T)
  write.csv(filtering.90.1, paste("nearly.universal.",Filt.Near.Uni,".percent.",ssp.treshold,"ssp.treshold",".csv", sep = ""))
  log_print(paste("File"," nearly.universal.",Filt.Near.Uni,".percent.",ssp.treshold,"ssp.treshold",".csv", " was created created", sep = ""))
  list.40_orthogroups_counts<-FamData[colnames(FamData) %in% colnames(filtering.90.1)]
  sppnumero<-dim(list.40_orthogroups_counts)[2]
  write.csv(list.40_orthogroups_counts, paste("Total.counts.",sppnumero,"spp.csv", sep = ""))
  log_print(paste("The file Total.counts.",sppnumero,"spp.csv was created", sep = ""))
}

### function to filter orthogroups with gene counts
Genome.completion.filtering(FamData = gene.numbers.filtered, Filt.Near.Uni = 90, ssp.treshold = 90)

##################################################
#** The Trait.orthoGr.filtering function filters and aligns trait and gene number datasets *
#** for downstream comparative analyses. It removes species with missing traits, ensures *
#** species alignment between datasets, and filters out gene families with excessive missing values, *
#** low variance, or presence in only a single species.*
###################################################

Trait.orthoGr.filtering<-function(traits, pheno.col.nums, spp.col.num, gene.numbers){
  phenotypes<- colnames(x = traits)[pheno.col.nums]
  for (i in 1:length(phenotypes)) {
    traits.filtered<<-traits[!is.na(traits[,phenotypes[i]]),]
  }
  spp.name<-colnames(x = traits)[spp.col.num]
  traits.filtered<<-traits.filtered[traits.filtered[,spp.name] %in% colnames(gene.numbers),]
  gene.numbers.filtered<-gene_numbers[,colnames(gene_numbers) %in% traits.filtered[,spp.name],]
  #gene.numbers.filtered<<-gene_numbers[colnames(gene_numbers) %in% traits.filtered$spp.name,]
  log_print("***********************************************************", hide_notes = T)
  log_print("****Spp missing in traits but present in gene.numbers:****", hide_notes = T)
  log_print("***********************************************************", hide_notes = T)
  log_print(setdiff(colnames(gene.numbers),traits.filtered[,spp.name]))
  log_print("***********************************************************", hide_notes = T)
  log_print("****Spp missing in gene.numbers but present in traits:****", hide_notes = T)
  log_print("***********************************************************", hide_notes = T)
  log_print(setdiff(traits$Spp.Name, traits.filtered$Spp.Name))
  if (dim(traits.filtered)[1] == 0) {
    log_print("The col. names from gene.number might be different to traits spp names", hide_notes = T)
  }
  log_print("****************************************************************************", hide_notes = T)
  log_print("Removing the orthogroups (gene families) that have 0 genes in 20% of the spp", hide_notes = T)
  log_print("****************************************************************************", hide_notes = T)
  # Keep only genes with a minimum number of missing values - in this case 20
  # As we have 40 species we want to keep rows with 20 or less zeros ??? ask araxi about this...
  ### change 20 to maximum zeros allowed e.g. 20/40 or 30/50 the diference should be 20...
  #Filter gene families to keep only those which have at least one gene in at least 6 species.
  
  #[huki] Remove any fams with more than 20% of zeros.
  zeros20perc <- ncol(gene.numbers.filtered)*0.2
  gene.numbers.filtered<-gene.numbers.filtered[rowSums(gene.numbers.filtered == 0) <= zeros20perc, ]
  
  log_print("************************************************************", hide_notes = T)
  log_print("Removing any orthogroupos (gene families) with variance of 0", hide_notes = T)
  log_print("************************************************************", hide_notes = T)
  #[huki]Remove any fams with variance of zeros.
  gene.numbers.filtered<-gene.numbers.filtered[apply(gene.numbers.filtered, 1, var) > 0,]
  
  #[huki]Remove any fams where the max number of genes in a row is 1
  # to avoid gene families present in only one sp. 
  log_print("****************************************************************************", hide_notes = T)
  log_print("Removing orthogroups (gene families) of 1 gene or less present in only 1 spp", hide_notes = T)
  log_print("*****************************************************************************", hide_notes = T)
  gene.numbers.filtered<<-gene.numbers.filtered[apply(gene.numbers.filtered, 1, max) > 2,]
  log_print(paste("The file ", paste("Filtered.GeneNum.",phenotypes,".csv"," was created", sep = "")))
  
  
  write.csv(gene.numbers.filtered, paste("Filtered.GeneNum.",phenotypes,".csv", sep = ""))
}

Trait.orthoGr.filtering(traits = traits.filtered, pheno.col.nums = 13, spp.col.num = 2, gene.numbers = gene_numbers)

#################
#** PGLS FUNCTION *
#################

PGLS.GF.size<- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it){
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
  
  nameout<- paste("PGLS", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits2)[1], "Spp", dim(out)[1], "GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  print(paste("object created:", paste( "out", sep = "")), hide_notes = T)
  write.csv(out, paste(where.Save.it, nameout, sep = "/"))
}

### remove two species

#sppout<-c("Dasypus_novemcinctus", "Ailuropoda_melanoleuca")
#traits.filtered2out <- traits.filtered[!traits.filtered$Spp.Name %in% sppout, ]
#tree2out<- drop.tip(tree, tip = sppout)
#gene.numbers2out<- gene.numbers.filtered[!colnames(gene.numbers.filtered) %in% sppout]

### PGLS function
dirSave<-path1

print(paste("You have ",dim(gene.numbers.filtered)[1], " gene families to input the PGLS", sep = ""), hide_notes = T)

## sort gene numbers by traits.filtered
gene.numbers.filtered <- gene.numbers.filtered[, match(traits.filtered$Spp.Name, colnames(gene.numbers.filtered))]

# write.csv(traits.filtered, "traits.filtered_SSDLog10testessize_paper.csv")
PGLS.GF.size(traits = traits.filtered, pheno.col.nums = c(12,13), No.variables = 2, 
             tree = tree, spp.col.num = 2, gene.numbers = gene.numbers.filtered, 
             where.Save.it = dirSave)

colnames(out)
# Calculate counts and percentages
summary_data <- out %>%
  filter(`p.adjusted GFS vs.log10MLSP` < 0.05) %>% ### CHANGE THIS NAME TO THE TRAIT YOU ARE INTERESTED IN
  summarise(
    Positive = sum(`R t value log10MLSP` > 0), ### ALSO CHANGE R VALUES
    Negative = sum(`R t value log10MLSP` < 0)  ### ALSO CHANGE R VALUES
  ) %>%
  tidyr::pivot_longer(everything(), 
                      names_to = "Direction", 
                      values_to = "Count") %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Create horizontal 100% stacked barplot
Nombeplot<- paste("Gene Families Sig.Associated (p < 0.05)_", trait, "_correctedByLog10MLSP_PGLS.pdf", sep = "")
pdf(Nombeplot, width = 8, height = 6, useDingbats = FALSE)

ggplot(summary_data, aes(y = "Significant Genes", x = Percentage, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Negative" = "indianred", "Positive" = "steelblue")) +
  geom_text(aes(label = sprintf("%.1f%% \n(n=%d)", Percentage, Count)), 
            position = position_stack(vjust = 0.5),
            size = 4) +
  labs(title = paste("Proportion of Gene Families with p < 0.05", trait, "PGLS"),
       x = "Percentage",
       y = "") +
  theme_minimal() +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  theme(legend.position = "top",
        legend.title = element_blank())
dev.off()

### output object called 'out'
#traits.filtered$Spp.Name
#out<- read.csv("PGLS.Av.Body.mass.125.Spp.13753.GeneFams.Rs_benjamini.2022-08-03.csv")
#perms<- read.csv("PGLS.Av.Body.mass.125.Spp.13753.GeneFams.Rs_benjamini.2022-08-03.csv")


