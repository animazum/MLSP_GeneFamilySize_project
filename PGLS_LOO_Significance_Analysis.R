

#** PGLS_LOO_Significance_Analysis.R *
# 
# Description:
# This script processes and analyzes Phylogenetic Generalized Least Squares (PGLS) results 
# from a Leave-One-Out (LOO) analysis. It combines multiple CSV files, cleans and structures 
# the data, and performs statistical tests to assess the significance of PGLS results.
#
# Key steps include:
# - Loading and combining PGLS results from multiple CSV files.
# - Cleaning and formatting the dataset.
# - Performing Shapiro-Wilk normality tests.
# - Conducting Wilcoxon tests to compare significant PGLS results with the original dataset.
# - Calculating effect sizes using Cohen's d.
# - Generating summary statistics and frequency tables.
#
# Dependencies:
# - dplyr (for data manipulation)
# - effsize (for effect size calculations)
#
# Author: [BEnjamin]
# Date: [Feb-2025]
#

# filename <- file.choose()
# Canteen_clean <- readRDS(filename)
rm(list=ls())

# setwd("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_1000reps_LOO")  # Adjust based on your actual folder
setwd("~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_LOO_sys_errorhandding/PGLS")  # Adjust based on your actual folder

files <- list.files(pattern = "\\.csv$")
library(dplyr)  # Useful for data manipulation

data_list <- lapply(files, read.csv, stringsAsFactors = FALSE)
names(data_list) <- files  # Assign file names to the list elements

combined_data <- bind_rows(data_list, .id = "source")

str(combined_data)
summary(combined_data)
head(combined_data)
colSums(is.na(combined_data))

cleaned_data <- combined_data %>%
  mutate(source = gsub("^PGLS_LOO\\.\\d+\\.|\\.cEi_log10MLSP.*$", "", source))

sumatoriaSPPS <- data.frame(Source = names(table(cleaned_data$source)), 
                            Freq = as.numeric(table(cleaned_data$source)))
mean(sumatoriaSPPS$Freq)
shapiro.test(sumatoriaSPPS$Freq)

originalPGLS<-read.csv("/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/PGLS_results/PGLS.cEi_log10MLSP.46.Spp.4136.GeneFams.Rs_benjamini.2025-02-05.csv")
dim(originalPGLS)[1]
mean(originalPGLS$R.t.value.log10MLSP)
originalmlspPval<-originalPGLS[originalPGLS$p_value.log10MLSP < 0.05,]
originalPGLSmlspPval<-originalPGLS[originalPGLS$p.adjusted.GFS.vs.log10MLSP < 0.05,]
originalPGLSeiPval<-originalPGLS[originalPGLS$p.adjusted.GFS.vs.cEi < 0.05,]

##** is the sample normal? *
shapiro.test(originalPGLSmlspPval$R.t.value.log10MLSP)

# Split the dataframe based on the "source" column
split_data <- split(cleaned_data, cleaned_data$source)

# Check how many species (sources) are present
length(split_data)  # Number of unique sources
names(split_data)   # Names of all sources

###** Wilcox of all PGLS results from LOO analysis. Corrected by MLSP *
###** Values of 0 - 0.2 is negligible so the groups are very similar *
###*
library(effsize)

wilcoxR2 <- lapply(names(split_data), function(x) {
  cat("\nWilcox-test of significant PGLS results from LOO analysis for", x, "\n")
  # Filter significant values within each species' dataframe
  SignificantPvalsR2mlsp <- subset(split_data[[x]], p_value.log10MLSP < 0.05)
  SignificantR2mlsp <- subset(split_data[[x]], p_value.log10MLSP < 0.05)
  cohentest<-cohen.d(SignificantPvalsR2mlsp$Coef.tval.log10MLSP, originalmlspPval$Coef.tval.log10MLSP)
  cat("Cohen test for effect size","Estimate", cohentest$estimate, "| Conf.int:", cohentest$conf.int, "| Magnitude:",cohentest$magnitude,"\n")
  
  # Ensure there are enough significant values to run the t-test
  if (nrow(SignificantPvalsR2mlsp) > 1) {
    wilcox_result <- wilcox.test(SignificantPvalsR2mlsp$R.t.value.log10MLSP, originalmlspPval$R.t.value.log10MLSP)
    
    # Return a named list containing both species name and the t-test result
    return(list(species = x, wilcox = wilcox_result))
  } else {
    return(list(species = x, wilcox = NA))  # Return NA if not enough values
  }
})

###** Check if pvalues are significat for each species *
for (i in 1:length(wilcoxR2)) {
  # Check if wilcox test result is not NULL
  if (!is.null(wilcoxR2[[i]]$wilcox)) {
    
    # Check if wilcox test result is a list or atomic vector
    if (is.list(wilcoxR2[[i]]$wilcox) && !is.null(wilcoxR2[[i]]$wilcox$p.value)) {
      p_value <- wilcoxR2[[i]]$wilcox$p.value[1]  # Extract safely
    } else if (is.atomic(wilcoxR2[[i]]$wilcox)) {
      p_value <- wilcoxR2[[i]]$wilcox  # Directly use value if it's an atomic vector
    } else {
      p_value <- NA  # Handle unexpected cases
    }
    
    # Check if p_value is valid
    if (!is.na(p_value)) {
      if (p_value < 0.05) {
        cat("\n", wilcoxR2[[i]]$species, "has a significant p-value of", p_value, "\n")
      } else {
        cat("\n", wilcoxR2[[i]]$species, "has a NOT significant p-value of", p_value, "\n")
      }
    } else {
      cat("\n", wilcoxR2[[i]]$species, "has no valid test result (NA or missing values). \n")
    }
    
  } else {
    cat("\n", wilcoxR2[[i]]$species, "has no valid test result (NULL value). \n")
  }
}

# Extract selected dataframes and combine them
do.call(rbind, split_data[c(16)])
SignificantPvalsR2 <- do.call(rbind, split_data[c(14,15,35,45)])
SignificantPvalsR2mlsp<- SignificantPvalsR2[SignificantPvalsR2$p.adjusted.GFS.vs.log10MLSP < 0.05,]
SignificantPvalsR2ei<- SignificantPvalsR2[SignificantPvalsR2$p.adjusted.GFS.vs.cEi < 0.05,]

# Create a frequency table for the "source" column
table(SignificantPvalsR2$source)
