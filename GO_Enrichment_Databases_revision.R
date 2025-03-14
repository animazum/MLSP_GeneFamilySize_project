

########################################################################################
#*** We are going to create the script to run the GO enrichment analysis for the MLSP **
# 
# Script to measure if a list of Families is particularly enriched in genes pertaining to a Gene Onthology Term
########################################################################
##*** ONLY VARIABLES THAT NEEDS TO BE CHANGED :                     *###
#### infile , GOanott,  phenotype, filenameout 
#### Those variables reduced to one so we also reduce the redundancy...
#########################################################################

###########################################
##### Script for GO enrichment relative testes size
###########################################

rm(list=ls())

setwd("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/Background and gene list/filtered/")
setwd("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/Background and gene list/GO_enrichment_GeneFamilies/")
PATH<- getwd()
setwd(PATH)

library(parallel)
library(colorout)


## one or two phenotype SSD PGLS

# List all CSV files in the directory
csv_files <- list.files(path = PATH, pattern = "*.csv", full.names = TRUE)


# Read each CSV file into a list, keeping the file names as identifiers
#BG_genelists <- setNames(lapply(csv_files, read.csv, stringsAsFactors = FALSE), 
                         #tools::file_path_sans_ext(basename(csv_files)))

BG_genelists <- readRDS("BG_genelists_processed.rds")

###************************************###
###*** CHOSE a database to work with **###
print(names(BG_genelists))

  #rm(list=ls())
  #setwd("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/Background and gene list/filtered/")
  #PATH<- getwd()
  setwd(PATH)
  
  # List all CSV files in the directory
  #csv_files <- list.files(path = PATH, pattern = "*.csv", full.names = TRUE)
  
  # Read each CSV file into a list, keeping the file names as identifiers
  #BG_genelists <- setNames(lapply(csv_files, read.csv, stringsAsFactors = FALSE), 
   #                        tools::file_path_sans_ext(basename(csv_files)))
  
databaseName<-BG_genelists[1]
infile<-BG_genelists[[1]]
head(infile, 3)
colnames(infile)
infile$X <- NULL
infile<- unique(infile)
colnames(infile)

GOenrichment<-function(filename1, BPFilt, PHENO, DATABASE, METODO, focaltrait,FAMIDCol, focusGenes, nameFile) { 
  
  # cl <- makeCluster(detectCores()-1) 
  ###########*************
  infile<-filename1
  head(infile,3)
  #infile <-  subset(infile, select = -c(X.1)) #para quitar la columna extra de rownames numerados en la Perm
  
  ########################                                                ####################
  ######################## Paramenters section. PLEASE change as you need ####################
  ########################                                                ####################
  
  BPFilt= BPFilt  ## Index to filter of minimun genes per gene family 
  PHENO<- PHENO
  DATABASE <- DATABASE
  METODO<- METODO ###  select "fdr" or "bonferroni" methods for the analisi, 
  
  if (DATABASE == "GENES") {
    ## Ensembl GO terms
    # GO.TERMS <-"GO.terms.table.Ensembl.bioproc.human.csv" ### GO terms file to be used
    GO.TERMS <-"~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/GO_terms/GOids_to_EnsemblID_v.ensembl113_human_short.csv" ### GO terms file to be used
    print("You have selected Ensembl GOs", quote=F)
  } else if (DATABASE == "FAMILIES") {
    ## NCBI GO terms
    GO.TERMS <-"~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/GO_terms/GOids_Ensemblid_GeneFams.ensembl113-CLEAN.csv" ### GO terms file to be used
    print("you have selected NCBI but please go with Ensembl", quote=F)
  } else { 
    print(" please write NCBI or ENsembl in DATABASE")
  }
  
  VALS<-paste(gsub(x = nameFile, pattern = "csv", replacement = ""), PHENO, "Phenotype", sep = "") ### Notes to be included in the outfile's name
  
  
  print("Please change the 'FAMIDCol','Pvalcol' & 'Rvalcol' with the number of your PGLS's file respective columns", quote=F)
  colnames(infile)
  head(infile,3)
  FAMIDCol = FAMIDCol
  focusGenes <- focusGenes
  print(paste(colnames(infile[c(FAMIDCol,focusGenes)])))
  
  dir.create(file.path(gsub(nameFile, pattern = ".csv", replacement = "")))
  setwd(file.path(gsub(nameFile, pattern = ".csv", replacement = "")))
  
  phenotype<-toString(paste(VALS, ".",METODO,sep = ""))# Put Phenotype of interest###########
  print(phenotype)
  
  ###################################################################################################
  ###################################################################################################
  
  rm(list=setdiff(ls(), c("infile", "focaltrait","DATABASE", "GOanott", "phenotype", "filenameout", "PATH","BPFilt", "VALS", "METODO","GO.TERMS", "FAMIDCol","focusGenes", "nameFile")))
  #setwd("/Users/Administrator/Documents/PhD/Proyecto/GliaN_Paper_Obs/Pearson/PIC")
  library(parallel)
  
  cl <- makeCluster(detectCores()-1) 
  
  filterbp<-toString(paste(METODO,"BP",BPFilt, sep = "."))#Put filter and BP characteristics for name###########*************
  
  outfile <- toString(paste("GO_", focaltrait, "_", phenotype, "_", filterbp, sep=""))
  outfile1 <<- toString(paste(outfile, ".csv", sep=""))
  print(outfile1)
  numreps<-10000
  
  
  ##################GOanott
  AssocGF<-function(infile, Genes){
    print("Here we are going to select the associated genefamilies")
    ssd_assoc <<- infile[infile[,Genes] == 1 ,]
  }
  
  AssocGF(infile = infile, Genes = focusGenes)
  
  dim(ssd_assoc)
  head(ssd_assoc)
  
  print( "Selecting rows shared between the PGLS file and the associated Rs" , quote=F)
  head(infile,3)
  Familyset <- as.data.frame(infile[infile[,FAMIDCol] %in% ssd_assoc$X, FAMIDCol]) ### change number top where your family ID are
  dim(infile)
  print("Check if 'Familyset' selected the gene families", quote=F)
  
  names(Familyset)<- c("family")
  dim(Familyset)
  head(Familyset)
  
  #Population of all gene families studied (Background)
  dat<<-as.data.frame(infile[,FAMIDCol]) ### change number top where your family ID are
  names(dat)<- c("family")
  
  ### read the GO terms that you want to use for the analysis
  print("#Gene onthology annotations per family*******", quote=F)
  GOanott<-read.csv(paste(GO.TERMS, sep = "/"),header=TRUE,stringsAsFactors=FALSE) #Gene onthology annotations per family###########*************
  colnames(GOanott)<-c("family", "GO_accession", "GO_name", "GO_domain")
  lovedCols<-c("family", "GO_accession", "GO_name", "GO_domain")
  
  if (colnames(GOanott)[1] == lovedCols[1] |colnames(GOanott)[2] == lovedCols[2]|colnames(GOanott)[3] == lovedCols[3]|colnames(GOanott)[4] == lovedCols[4]) {
    print("all good, the col names of you your GO annotations are perfect, lets hope that the data is good aswell")  
  } else {
    print(paste("change your GO annotations colnames to be:", sep = ""), quote=F)
    print(lovedCols)
    print(" your col names look like this:", quote=F)
    print(colnames(GOanott))
  }
  
  ##### setting up the file to be used by the enrichment program
  print("Change the column names as in 'lovedCols'", quote=F)
  print("Change the column names as in 'lovedCols'", quote=F)
  print("Change the column names as in 'lovedCols'", quote=F)
  print("Change the column names as in 'lovedCols'", quote=F)
  print("Change the column names as in 'lovedCols'", quote=F)
  head(GOanott)
  
  ##################Get only GO terms that include at least n families
  print(paste("Getting only GO terms that include at least",BPFilt,"families"), quote=F)
  GOsize <-table(GOanott$GO_accession)
  head(GOsize, 10)
  
  bigGOs2<-names(GOsize[which(GOsize>=BPFilt)])
  head(bigGOs2, 10)
  print(paste("We have", length(bigGOs2), "out of", length(GOsize), "GO terms including", BPFilt, "families"), quote=F)
  out5<-GOanott
  head(out5)
  dim(out5)
  
  out5$GO_accession[!(out5$GO_accession %in% bigGOs2)]<-"GO:000STBP"
  out5$GO_name[!(out5$GO_accession %in% bigGOs2)]<-"small biological process GO terms"
  head(out5)
  dim(out5[(out5$GO_accession == "GO:000STBP"),])
  
  out5<-out5[!(out5$GO_accession == "GO:000STBP"),]
  head(out5)
  print(paste("We have", length(bigGOs2), "out of", length(GOsize), "GO terms including", BPFilt, "families"), quote=F)
  print(paste("You have",dim(out5)[1],"GF terms that include at least", BPFilt, "families.",sep=" "), quote=F)
  
  print(paste("output file name:", paste("GO.Terms.BP",BPFilt,phenotype,"csv", sep=".")), quote=F)
  write.csv(out5,paste("GO.Terms.BP",BPFilt,phenotype,"csv", sep="."),row.names=F)
  GOanott<- out5
  ################3 END
  
  ### eliminates the stbp witch are small gene fams. 
  dim(GOanott)
  GOanott<<-GOanott[!is.na(GOanott$GO_accession),]
  dim(GOanott)
  
  GOsBPn<<-unique(GOanott[,-1])
  dim(GOsBPn)
  GOsBP<<-GOsBPn[,1]
  dim(GOsBP)
  
  GOanott2<<-matrix(NA,dim(GOanott)[1],2,dimnames = list(c(1:dim(GOanott)[1]),c("GO_accession", "family")))
  GOanott2[,1]<<-GOanott[,2]
  GOanott2[,2]<<-GOanott[,1]
  GOanott2<<-GOanott2[!duplicated(GOanott2),]
  
  ### creating the family Id, GO ID dataframe
  counts<-vector("numeric")
  annotation<-vector("numeric")
  join<-merge(Familyset,GOanott2,by="family")
  head(join)
  print("Check that your families and GO ids are in the 'join' object", quote = F)
  head(join) 
  
  head(Familyset)
  w <<- dim(Familyset)[1]
  w
  
  #countssmpl<-matrix(NA,length(GOsBP),numreps)
  zscore<-matrix(NA,length(GOsBP),1,dimnames=list(GOsBP,"Z-score"))
  #annotationsmpl<-matrix(NA,w,numreps)
  countssmpl2<-matrix(NA,length(GOsBP),numreps)
  zscore2<-matrix(NA,length(GOsBP),1,dimnames=list(GOsBP,"Z-score"))
  
  
  counts<-sapply(GOsBP, function(x) #Counts in Familyset per GO
    sum(join[,2]==x))
  
  annotation<-apply(Familyset,1, function(x) #Annotation density of each family in Familyset
    sum(join[,1]==x))
  
  counts2<-counts/mean(annotation) #Counts in Familyset per GO corrected against annotation density
  
  ### start the Go enrichment analysis???
  clusterExport(cl, c("dat", "GOanott2", "GOsBP", "w"))
  run3 <- function(c) {
    smpl <- as.data.frame(sample(dat[,1],w),stringsAsFactors=FALSE)
    names(smpl)<- c("family")
    joinsmpl<-merge(smpl,GOanott2,by="family")
    
    countssmpl<-sapply(GOsBP, function(x) #Cuentas por GO Term de genes en la smpl
      sum(joinsmpl[,2]==x))
    
    annotationsmpl<-apply(smpl,1, function(x) #Densidad de anotacion en reales
      sum(joinsmpl[,1]==x))
    
    
    return(list(countssmpl=countssmpl, annotationsmpl=annotationsmpl))
    print("GO enrichment done??", quote=F)
  }
  
  MCsmpl<-parSapply(cl,rep(1,numreps),run3)
  stopCluster(cl)
  
  countssmpl<-do.call(cbind,MCsmpl[1,])
  annotationsmpl<-do.call(cbind,MCsmpl[2,])
  
  means<-rowMeans(countssmpl)
  meansannotation<-colMeans(annotationsmpl)
  for(p in c(1:numreps))
    countssmpl2[,p]<-countssmpl[,p]/meansannotation[p]
  means2<-rowMeans(countssmpl2)
  SEMs<-apply(countssmpl,1,sd)
  SEMs2<-apply(countssmpl2,1,sd)
  zscore[,1]<-(counts-means)/SEMs
  zscore2[,1]<-(counts2-means2)/SEMs2
  pval<-pnorm(-abs(zscore))
  adjpval<-p.adjust(pval, METODO)
  #adjpval<-p.adjust(pval,"bonferroni")
  pval2<-pnorm(-abs(zscore2))
  adjpval2<-p.adjust(pval2,METODO)
  #adjpval2<-p.adjust(pval2,"bonferroni")
  numpval<-rowSums(countssmpl>=counts)/numreps
  numpval2<-rowSums(countssmpl2>=counts2)/numreps
  adjnumpval<-p.adjust(numpval, METODO)
  #adjnumpval<-p.adjust(numpval,"bonferroni")
  adjnumpval2<-p.adjust(numpval2, METODO)
  #adjnumpval2<-p.adjust(numpval2,"bonferroni")
  
  out1<-as.data.frame(cbind(GOsBPn[,1],GOsBPn[,2],counts,means,SEMs,zscore,pval, adjpval, numpval, adjnumpval),stringsAsFactors=FALSE)
  #out2<-as.data.frame(cbind(GOsBPn[,1],GOsBPn[,2],counts2,means2,SEMs2,zscore2,pval2, adjpval2, numpval, adjnumpval),stringsAsFactors=FALSE)
  
  names(out1)<-c("GOID","GO Term","Genes in the GO","Expected genes per GO","Standard Deviation of the samples","Z-score","One tail p-value","FDR", "Numeric p-value", "Numeric FDR")
  #names(out2)<-c("GOID","GO Term","Gene Families in the GO","Expected families per GO","Standard Deviation of the samples","Z-score","One tail p-value","FDR", "Numeric p-value", "Numeric FDR")
  head(out1,4)
  
  #write.csv(out1, outfile, row.names=F)
  print("######################################")
  print(phenotype, quote=F)
  print(paste("We have", length(bigGOs2), "out of", length(GOsize), "GO terms including", BPFilt, "families"), quote=F)
  
  print(paste("You have",dim(out5)[1],"GF terms that include at least", BPFilt, "families.",sep=" "), quote=F)
  
  print(paste("output file name:", paste("GO.Terms.BP",BPFilt,phenotype,"csv", sep=".")), quote=F)
  
  print(paste("Your file is called: ",outfile1, sep = ""), quote =F)
  
  print(paste("You have ", sum(as.numeric(out1$FDR) < 0.05), " significant FDR values", sep = ""))
  phenotype
  
  write.csv(out1, outfile1, row.names=F)
  print(paste("Check the FDR to be less than 0.05"), quote =F)
  print("########################################")
  alarm()
  closeAllConnections()
}
print(paste("Please check your FAMIDCol, Pvalcol and Rvalcol:"))
colnames(infile)

###****************Change biological process filtering, PHENO= the number of phenotypes, focaltrait = name of the main trait of interest ***************
###****************DATABASE = the name of the database used, METODO = the statistical method to be used, FAMIDCol= gene family ID column number ***************
###**************** Pvalcol = the p value column number Rcalcol = the number of the R value column  ***************
GOenrichment(filename1 = infile, BPFilt = 50, PHENO = 2, focusGenes = 2,
             DATABASE = "FAMILIES", METODO = "fdr", FAMIDCol = 1,focaltrait = "MLSP_Ei", nameFile = names(databaseName))


print("###########################")
print("####### plot plot plot#####")
print("######## plot heat map ####")
print("###########################")
BPFilt <- 50
METODO = "fdr"

archivos <- list.files(pattern = paste("^",outfile1,"$", sep = ""))
archivos 
#VAL <- "negative"
nameout <- paste("All_Groups")####put name of file TodosM4top20
filterchar<-toString(paste("BP",BPFilt,"_", METODO, sep = ""))#Put filter and BP characteristics for name###########*************
outfilename <- toString(paste("GO-HM.",gsub(archivos, pattern = ".csv", replacement = ""), nameout, "_", filterchar, ".", sep=""))
outfilename
#*****####If one of the "archives" is modified before, then quit manually single or double quotes (' or ") from the ID names so it doesn't complain when reading the table.*******##
##fill the ones to omit with NAs###

HEATMAPNAME <- paste("1.1heatmap.",archivos,sep = "")
HEATMAPNAME

noms<-sapply(strsplit(archivos,"_"),"[",2) 

tb1<-read.table(archivos[1],sep=",",header=TRUE) #### change number according to the file that you want to use 

dim(tb1)
head(tb1)

print(paste("You have ", sum(as.numeric(tb1$FDR) < 0.05), " significant FDR values", sep = ""))

TABLAGOenrich<-tb1[(tb1$FDR) < 0.05,]
tb1<-tb1[!is.na(tb1$GO.Term),]
dim(tb1)

ps<-matrix(NA,dim(tb1)[1],length(archivos),dimnames=list(tb1$GO.Term,noms))
zs<-matrix(NA,dim(tb1)[1],length(archivos),dimnames=list(tb1$GO.Term,noms))
obs<-matrix(NA,dim(tb1)[1],length(archivos),dimnames=list(tb1$GO.Term,noms))
for (i in c(1:length(archivos))){
  tb1<-read.table(archivos[i],sep=",",header=TRUE)
  #tb1<-tb1[order(tb1$GOID),]
  ps[,i]<-tb1$FDR
  zs[,i]<-tb1$Z.score
  obs[,i]<-tb1$Genes.in.the.GO
}
head(ps)
dim(ps)
### positive z scores
GOn<-apply(zs,2,function(x) length(which(x>0))) #valores de zscore positivos (por lo tanto solo categorias overrepresented), tanto significativos como no

GOn2<-apply(ps,1,function(x) length(which(x<0.05)))
ind <- which(GOn2==0)

newnames<-paste(colnames(ps),GOn)
colnames(ps)<-newnames
ps[which(ps>0.05)]<-NA
ps[which(zs<0)]<-NA
ps[which(obs<3)]<-NA #filter that more than 1 Gene Family is observed

head(ps)
pscol<-colnames(ps)
PS1<-ps
ps <- ps[-ind,]
ps<- sort(ps, decreasing = F)

ps<-as.data.frame(ps)
colnames(ps)<- pscol

head(ps)

tablename<-toString(paste(outfilename, "_Table", sep=""))
print(paste("Name of the GO enrichment table: ", tablename, ".csv",sep = ""), quote =F)
write.csv(ps, paste(tablename, ".csv", sep=""))


ps<-read.csv(toString(paste(tablename, ".csv", sep="")), row.names = 1)
# the heatmap with all SSD alone and SSD + BM
# ps<-read.csv("/Users/lisalisa/Library/CloudStorage/Dropbox/SSD/GO-HM.GO_SSD_ALLPositiveNegative.GO.terms.table.Process.mammals.ALLPhenotype.fdr_fdr.BP.50All_Groups_BP50_fdr._TableOCT2023.csv", row.names = 1, )
# the heatmap with SD + BM
# ps<-read.csv("~/Dropbox/TestesSize_SSD/Datos/PGLS.SSD.77.Spp.4251.GeneFams.Rs_benjamini.2024-11-06/GO-HM.GO_SSD_Negative.GO.terms.table.Process.mammals.1Phenotype.fdr_fdr.BP.50All_Groups_BP50_fdr._Table.csv", row.names = 1, )

#Heatmap plot #
library(reshape)
library(lattice)
library(RColorBrewer)
# colores<-c(rev(brewer.pal(9,"YlOrRd")[2:9]))
colores<-c(rev(brewer.pal(9,"YlGnBu")[2:9]))

paleta<-colorRampPalette(colores, space = "Lab")

colores<-paleta(50)
interv<-seq(0,0.05,by=.001)

HEATMAPNAME<-"heatmap_MLSP.Ei_GOenrichment.csv"
pdfname <- toString(paste(HEATMAPNAME,".pdf", sep=""))
pdf(pdfname, paper="a4r",width =8, height=8) 
levelplot(t(ps), main="", xlab="", ylab="", col.regions=colores, cuts=100, at=interv, scales=list(x=list(rot=90,cex=.7), y=list(rot=0,cex=0.7), draw=TRUE) , colorkey=list(TRUE,space="right",contour=FALSE,pretty=TRUE),   par.settings = list(axis.line = list(col = 0)),
          panel = function(...) {
            panel.fill(col = "gray96")
            panel.levelplot(...)
            #panel.text(x,y, a, cex=.6)
          })
dev.off()

png(paste(HEATMAPNAME,".png",sep = "") ,units="px", width=1920, height=1920, res=290)
levelplot(t(ps), main="", xlab="", ylab="", col.regions=colores, cuts=100, at=interv, scales=list(x=list(rot=90,cex=.7), y=list(rot=0,cex=0.7), draw=TRUE) , colorkey=list(TRUE,space="right",contour=FALSE,pretty=TRUE),   par.settings = list(axis.line = list(col = 0)),
          panel = function(...) { 
            panel.fill(col = "gray96")
            panel.levelplot(...)
            #panel.text(x,y, a, cex=.6)
          })
dev.off()

print(paste("Name of the GO enrichment table: ", tablename, ".csv",sep = ""), quote =F)
print(paste("Name of the plot in PDF: ",HEATMAPNAME, ".pdf",sep = ""), quote=F)
getwd()

