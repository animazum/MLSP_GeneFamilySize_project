
library(parallel)
library(colorout)

rm(list=ls())

setwd("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/Background and gene list/GO_enrichment_GeneFamilies/")
PATH<- getwd()
setwd(PATH)
## one or two phenotype SSD PGLS

# List all CSV files in the directory
csv_files <- list.files(path = PATH, pattern = "*.csv", full.names = TRUE)


# Read each CSV file into a list, keeping the file names as identifiers
BG_genelists <- setNames(lapply(csv_files, read.csv, stringsAsFactors = FALSE), 
                         tools::file_path_sans_ext(basename(csv_files)))


GO.TERMS <-"~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/GO_terms/GOids_Ensemblid_GeneFams.ensembl113.csv" ### GO terms file to be used
GOanott<-read.csv(paste(GO.TERMS, sep = "/"),header=TRUE,stringsAsFactors=FALSE) #Gene onthology annotations per family###########*************
GOanott<- GOanott[c(1,2)]
GOanott$Gene.stable.ID<- NULL
GOanott <- unique(GOanott)
GOanott <- GOanott[GOanott$PANTHER.ID != "", ]
write.csv(GOanott, "~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/GO_terms/GOids_Ensemblid_GeneFams.ensembl113.csv")

# Function to process each dataframe
process_dataframe <- function(df, GOanott) {
  # Step 1: Merge with GOanott
  merged_df <- merge(x = df, y = GOanott, by.x = "X", by.y = "Gene.stable.ID", all.x = TRUE, all.y = FALSE)
  
  # Step 2: Remove rows with NA values in any column
  merged_df <- na.omit(merged_df)
  
  # Step 3: Remove rows with NA in specific columns (e.g., PANTHER.ID)
  merged_df <- merged_df[!is.na(merged_df$PANTHER.ID), ]
  
  # Step 4: Remove rows with blank cells in specific columns (e.g., PANTHER.ID)
  merged_df <- merged_df[merged_df$PANTHER.ID != "", ]
  
  # Step 5: Change order of columns
  merged_df <- merged_df[c(1,4,2,3)]
  
  # Step 6: Remove duplicate rows
  merged_df <- unique(merged_df)
  
  return(merged_df)
}

# Apply the function to each dataframe in the list using lapply
BG_genelists_processed <- lapply(BG_genelists, process_dataframe, GOanott = GOanott)

saveRDS(BG_genelists_processed, file = "BG_genelists_processed.rds")


