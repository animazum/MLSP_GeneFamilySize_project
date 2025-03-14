
rm(list=ls())

# Set the working directory (if not already set)
PATH <- getwd()

# List all CSV files in the directory
csv_files <- list.files(path = PATH, pattern = "*.csv", full.names = TRUE)

# Read each CSV file into a list, keeping the file names as identifiers
BG_genelists <- setNames(lapply(csv_files, read.csv, stringsAsFactors = FALSE), 
                         tools::file_path_sans_ext(basename(csv_files)))

# Merge all data frames in the list by the "X" column
merged_df <- Reduce(function(x, y) merge(x, y, by = "X", all = TRUE), BG_genelists)
# colnames(merged_df)<-c("X","MLSP_Ei_1","MLSP_Ei_2","MLSP_Ei_3", "MLSP_Ei_4", "MLSP_Ei_5", 
#                    "MLSP_Ei_6", "MLSP_Ei_7", "MLSP_Ei_8", "MLSP_Ei_9", "MLSP_Ei_10", 
#                    "MLSP_Ei_11", "MLSP_Ei_12", "MLSP_Ei_13", "MLSP_Ei_14", "MLSP_Ei_15")

colnames(merged_df) <- c("X", "AgeCellExp_Neg", "AgeCellExp_Pos", "Turnover_Neg", "Turnover_Pos", 
                         "AgeExp_Neg", "AgeExp_Pos", "Apoptosis", "Senescence", 
                         "DNA_Repair", "Diet_Resist", "Ei_Assoc", "Longevity_Var", 
                         "LifeExt_Drugs", "MLSP_Assoc", "PM_Longevity")
# View the merged data frame
head(merged_df)


write.csv(merged_df, "00.GO-HM.GO_MLSP_Ei_ALL_GOenrichments.csv", row.names = F)

ps<-read.csv("00.GO-HM.GO_MLSP_Ei_ALL_GOenrichments.csv", row.names = 1)
ps <- as.matrix(ps)  # Ensure data is in matrix format
rownames(ps) <- rownames(ps)

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
levelplot(t(ps), main="", xlab="", ylab="", col.regions=colores, cuts=100, at=interv, scales=list(x=list(rot=90,cex=.5), y=list(rot=0,cex=0.5), draw=TRUE) , colorkey=list(TRUE,space="right",contour=FALSE,pretty=TRUE),   par.settings = list(axis.line = list(col = 0)),
          panel = function(...) {
            panel.fill(col = "gray96")
            panel.levelplot(...)
            #panel.text(x,y, a, cex=.6)
          })
dev.off()

png(paste(HEATMAPNAME,".png",sep = "") ,units="px", width=1920, height=1920, res=290)
levelplot(t(ps), main="", xlab="", ylab="", col.regions=colores, cuts=100, at=interv, scales=list(x=list(rot=90,cex=.5), y=list(rot=0,cex=0.5), draw=TRUE) , colorkey=list(TRUE,space="right",contour=FALSE,pretty=TRUE),   par.settings = list(axis.line = list(col = 0)),
          panel = function(...) { 
            panel.fill(col = "gray96")
            panel.levelplot(...)
            #panel.text(x,y, a, cex=.6)
          })
dev.off()

print(paste("Name of the GO enrichment table: ", tablename, ".csv",sep = ""), quote =F)
print(paste("Name of the plot in PDF: ",HEATMAPNAME, ".pdf",sep = ""), quote=F)
getwd()

