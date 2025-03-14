

#### lets create the backgrounds for the GO enrichment
library(tools)

rm(list=ls())

# Set the directory path
directory_path <- "~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/Background and gene list/"  # Change this to your actual directory

# List all CSV files in the directory
csv_files <- list.files(path = directory_path, pattern = "*.csv", full.names = TRUE)

# Read each CSV file into a list, keeping the file names as identifiers
BG_genelists <- setNames(lapply(csv_files, read.csv, stringsAsFactors = FALSE), 
                         tools::file_path_sans_ext(basename(csv_files)))

# Initialize a list to store the resulting data frames
BG_dfs <- list()

# Loop over each data frame in BG_genelists
for (file_name in names(BG_genelists)) {
  df <- BG_genelists[[file_name]]  # Get the data frame
  
  # Remove duplicate rows based on the first column (gene names)
  df <- df[!duplicated(df[, 1]), ]  
  
  # Extract the first column as row names
  row_names <- df[, 1]
  firstCol <- df[, 2]  # Background column (binary)
  
  # Create individual data frames for each column (except the first one)
  for (col_name in colnames(df)[-c(1,2)]) {
    new_df <- data.frame(focusGene = df[[col_name]])  # New data frame
    rownames(new_df) <- row_names  # Assign row names
    new_df$background <- firstCol  # Add background column
    
    # Store the new data frame in the list with a unique name
    BG_dfs[[paste(file_name, col_name, sep = "_")]] <- new_df
  }
}

# Apply filtering: Keep only rows where 'background' == 1
BG_dfs_filtered <- lapply(BG_dfs, function(df) df[df$background == 1, ])


# Define output directory (change to your preferred location)
output_directory <- "~/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/results/GO_enrichment_databases/Background and gene list/"  # Change this to your actual directory
dir.create(output_directory, showWarnings = FALSE)  # Create directory if it doesn't exist

# Save each filtered dataframe as a CSV file
for (df_name in names(BG_dfs_filtered)) {
  file_path <- file.path(output_directory, paste0(df_name, ".csv"))
  write.csv(BG_dfs_filtered[[df_name]], file_path, row.names = TRUE)  # Save with row names
}

# Print confirmation message
print(paste("Saved", length(BG_dfs_filtered), "filtered CSV files in", output_directory))


# Print names of the created and filtered data frames
print(names(BG_dfs_filtered))

# Return the filtered data frames
BG_dfs_filtered
