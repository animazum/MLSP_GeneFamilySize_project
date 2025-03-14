"""
Phylogenetic and Statistical Analysis of MLSLP Sexual Size Dimorphism (SSD) in Mammalian Species

Key Features and Steps:
1. **Data Setup**:
   - Sets the working directory and global variables for key parameters in the analysis.
   - Loads and cleans data on gene orthogroups and traits related to body mass and testes mass.
   
2. **Phylogenetic Tree Processing**:
   - Loads a phylogenetic tree in Newick format, cleans tip labels, and prunes the tree to retain only 
     relevant species based on the trait data.
   - Provides a visual representation of both the original and pruned trees and outputs the final number of 
     species.

3. **Trait Data Processing**:
   - Computes the predicted testes size based on body mass using a specified formula, then calculates 
     relative testes size by comparing actual to predicted values.
   - Removes identified outliers based on species names.

4. **Statistical Analysis**:
   - Tests for normality in relative testes size and log-transformed body mass using the Shapiro-Wilk test.
   - Adds log-transformed columns for key variables to support log-based analyses.
   - Performs correlation analyses between:
     - Log-transformed male body mass and log-transformed testes size.
     - Log-transformed testes size and SSD.
     - Relative testes size and SSD.
   - Outputs correlation coefficients, R² values, and p-values.

5. **Data Visualization**:
   - Creates and saves scatter plots with regression lines and annotated statistical results for each 
     correlation analysis:
       - Correlation between log-transformed body mass and log-transformed testes size.
       - Correlation between log-transformed testes size and SSD.
       - Correlation between relative testes size and SSD.

Dependencies:
- `pandas` for data handling.
- `numpy` for numerical calculations.
- `scipy.stats` for statistical tests.
- `matplotlib.pyplot` and `seaborn` for plotting.
- `Bio.Phylo` (from Biopython) for phylogenetic tree processing.

Usage:
Run this script in a Python environment with the necessary libraries installed. Ensure all input data files 
are accessible and paths are correctly set.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import os


# for name in dir():
#     if not name.startswith('_'):
#         del globals()[name]

# Set working directory
# os.chdir("~/Dropbox/TestisSize_SSD/Datos/")
os.chdir("/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/Datos")
# Global variables
#SSD_bias = "MALE"
trait = "MLSP"
correcting = "benjamini"

path1= os.getcwd()

# Read the orthogroups count data
orthogroups_counts = pd.read_csv("~/Dropbox/SSD/Orthogroups.GeneCount.tsv", sep="\t", index_col=0)
print(f"Species in the list: {orthogroups_counts.shape[1]}")

# Clean column names
orthogroups_counts.columns = [col.replace(".LS", "") for col in orthogroups_counts.columns]
orthogroups_counts = orthogroups_counts.drop("Total", axis=1)

# Read the testis data
Nombre_file = os.path.join(path1, "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/Datos/anage/anage_data.txt")

# Define full file path in one step
Nombre_file = os.path.join(
    "/Users/lisalisa/Library/CloudStorage/Dropbox/kilili_shared_project/MLSP_GeneFamilySize_project/Datos/anage",
    "anage_data.txt")
traits = pd.read_table(Nombre_file, encoding='ISO-8859-1')

# Merge columns with '_'
traits['species_name'] = traits['Genus'].astype(str) + '_' + traits['Species'].astype(str)
# Remove NaNs only from 'col1'
traits = traits.dropna(subset=["Maximum longevity (yrs)"])

# Subset DataFrame to keep only specific columns
traits = traits[['species_name', 'Maximum longevity (yrs)']]

# Define full file path in one step
Nombre_file2 = os.path.join(
    "~/Dropbox/SSD/",
    "SSD_mammalian_124sppWithOrtogroup_data_150721.csv")
traitsSSD = pd.read_csv(Nombre_file2, encoding='ISO-8859-1')

# # Subset DataFrame to keep only specific columns
traitsSSD = traitsSSD[['Spp.Name', 'SSD']]

# Merge using different column names
traits_MLSPSSD = pd.merge(traits, traitsSSD, left_on='species_name', right_on='Spp.Name', how='inner')
# Drop the duplicate column ('spp_name') if needed
traits_MLSPSSD = traits_MLSPSSD.drop(columns=['Spp.Name'])



# Read and process phylogenetic tree
tree = Phylo.read("SpeciesTree_rooted_node_labels.txt", "newick")
# Remove the ".LS" suffix from each tip label
for tip in tree.get_terminals():
    tip.name = tip.name.replace(".LS", "")
print("Original tree:")
Phylo.draw_ascii(tree)

# Load the traits data (assuming traits is a DataFrame that has been loaded already)
print(traits_MLSPSSD.columns)
spp_names = set(traits_MLSPSSD['species_name'])  # List of species names to keep

# Iterate over tree tips and prune those not in spp_names
tree2 = tree  # Make a copy to preserve the original tree
for tip in tree2.get_terminals():
    if tip.name not in spp_names:
        tree2.prune(tip)  # Removes the tip from tree2
print("Pruned tree:")        
# Display the pruned tree (optional)
Phylo.draw_ascii(tree2)
num_species = len(tree2.get_terminals())
print(f"Number of species in the tree: {num_species}")

# Remove outliers
#outliers = ["Pteropus_alecto", "Phyllostomus_discolor"]
#traits = traits[~traits['Spp.Name'].isin(outliers)]

# Add log-transformed columns
#traits_MLSPSSD.loc[:, 'Log10.SSD'] = np.log10(traits_MLSPSSD['SSD'])
#traits_MLSPSSD.loc[:, 'Log10.MLSP'] = np.log10(traits_MLSPSSD['Maximum longevity (yrs)'])

# Test for normality
#print("Shapiro-Wilk test for SSD:")
#print(stats.shapiro(np.log10(traits_MLSPSSD['SSD'])))
#print("\nShapiro-Wilk test for MLSP:")
#print(stats.shapiro(np.log10(traits_MLSPSSD['Maximum longevity (yrs)'])))

# Add log-transformed columns
#traits.loc[:, 'Log10.Rel.testes.size'] = np.log10(traits['Rel.testes.size'])
#traits.loc[:, 'Log10Bodymass_male.gr'] = np.log10(traits['Bodymass_male..g.'])
#traits.loc[:, 'Log10testes.size'] = np.log10(traits['Combined.Testes.Mass..g.'])

# Correlation analyses
import scipy.stats as stats

def perform_correlation(x, y, data):
    """
    Computes Pearson correlation between two columns in a Pandas DataFrame.

    Parameters:
    - x (str): Name of the first column.
    - y (str): Name of the second column.
    - data (pd.DataFrame): The DataFrame containing the columns.

    Returns:
    - dict: Pearson correlation coefficient (r), r-squared value, and p-value.
    """
    r_value, p_value = stats.pearsonr(data[x], data[y])  # Compute correlation
    r_squared = r_value ** 2  # Compute R² value

    return {"r": r_value, "r_squared": r_squared, "p_value": p_value}
  
traits_MLSPSSD.to_csv("traits_MLSPSSD.csv", index=False)


# Correlation 1: body mass vs testes size
corr1 = perform_correlation('SSD', 'Maximum longevity (yrs)', traits_MLSPSSD)
print("\nCorrelation between log male body mass and log Testes Size:")
print(print(f"R: {corr1['r']:.3f}, R²: {corr1['r_squared']:.3f}, p-value: {corr1['p_value']:.3e}"))


# Plotting function
def create_correlation_plot(data, x, y, title, xlabel, ylabel, r_squared, p_value, filename):
    plt.figure(figsize=(10, 6))
    
    # Scatter plot
    sns.scatterplot(data=data, x=x, y=y, alpha=0.7)  # Added transparency for better visibility
    
    # Regression line
    sns.regplot(data=data, x=x, y=y, scatter=False, color='blue')

    # Labels and title
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Add statistics annotation
    stats_text = f'R² = {r_squared:.3f}\np-value = {p_value:.3e}'
    plt.text(0.95, 0.90, stats_text, 
             transform=plt.gca().transAxes,
             ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# Create plot
create_correlation_plot(
    traits_MLSPSSD, 'SSD', 'Maximum longevity (yrs)',
    'Correlation between SSD and MLSP',
    'SSD', 'Maximum longevity (yrs)',
    corr1['r_squared'], corr1['p_value'],  # Ensure correct dictionary keys
    'Correlation_pearson_SSDvsMLSP.pdf'
)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr


# Scatterplot with LOESS smoothing curve
def create_nonlinear_plot(data, x, y, title, xlabel, ylabel, filename):
    plt.figure(figsize=(10, 6))

    # Scatterplot
    sns.scatterplot(data=data, x=x, y=y, alpha=0.7)

    # LOESS smoothing (nonparametric trend)
    sns.regplot(data=data, x=x, y=y, scatter=False, lowess=True, color='red', line_kws={"lw": 2})

    # Polynomial trend line (degree=2 for quadratic fit)
    sns.regplot(data=data, x=x, y=y, scatter=False, order=2, color='blue', line_kws={"lw": 2})

    # Labels and title
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# Generate plot
create_nonlinear_plot(
    traits_MLSPSSD, 'SSD', 'Maximum longevity (yrs)',
    'Nonlinear Relationship: SSD vs. MLSP',
    'SSD', 'Maximum longevity (yrs)',
    'Nonlinear_Correlation_SSDvsMLSP.pdf'
)  
