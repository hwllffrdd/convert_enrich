# Proteomics Pathway Enrichment Analysis

This R script performs pathway enrichment analysis on proteomics data using the ReactomePA package. It converts UniProt accessions to Entrez IDs and performs enrichment analysis for both up- and down-regulated proteins.

## Prerequisites

The script requires R (>= 4.0.0) and the following packages:

### CRAN packages:
- dplyr

### Bioconductor packages:
- BiocManager
- biomaRt
- ReactomePA
- enrichplot
- org.Hs.eg.db

## Installation

1. First, ensure you have R installed on your system.

2. Clone this repository:
```bash
git clone https://github.com/hwllffrdd/convert_enrich.git
cd convert_enrich
```

3. The script will automatically install required packages if they're not present. However, you can manually install them:

```R
# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("biomaRt", "ReactomePA", "enrichplot", "org.Hs.eg.db"))

# Install CRAN packages
install.packages("dplyr")
```

## Input Data Format

The script expects a CSV file with semicolon (;) as separator. Your CSV file should contain the following columns (column names can be customized):

Required columns:
- UniProt accessions (e.g., "uniprot_id", "accession", etc.)
- Log2 fold change values (e.g., "log2_fold_change", "l2fc", etc.)
- P-values (e.g., "p_value", "pval", etc.)

Example CSV format:
```csv
uniprot_id;log2_fold_change;p_value;other_columns...
P12345;1.5;0.001;...
P67890;-2.3;0.003;...
```

Notes:
- Fold changes should be in log2 scale
- By default:
  - Proteins with log2FC ≥ 1 (FC ≥ 2) and p-value ≤ 0.05 are considered upregulated
  - Proteins with log2FC ≤ -1 (FC ≤ 0.5) and p-value ≤ 0.05 are considered downregulated
- These thresholds can be customized in the function call

## Usage

Basic usage with default column names:
```R
source("proteomics_enrichment.R")
results <- main(
  data_file = "your_data.csv",
  fc_up_col = "log2_fold_change",
  p_value_col = "p_value",
  uniprot_col = "uniprot_id",
  output_prefix = "your_analysis"
)
```

Advanced usage with custom settings:
```R
results <- main(
  data_file = "your_data.csv",
  fc_up_col = "L2FC",                # your fold change column name
  p_value_col = "pval",              # your p-value column name
  uniprot_col = "Accession",         # your UniProt column name
  fc_up_cutoff = 1.5,                # custom fold change threshold
  fc_down_cutoff = -1.5,             # custom fold change threshold
  p_value_cutoff = 0.01,             # custom p-value threshold
  output_prefix = "your_analysis"
)
```

## Output

The script generates:
1. CSV files with enrichment results for up- and down-regulated proteins
2. Dotplots visualizing the enrichment results
3. Returns a list containing the full enrichment results and plots
