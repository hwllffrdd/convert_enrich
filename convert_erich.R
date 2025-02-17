#!/usr/bin/env Rscript

#' Proteomics Pathway Enrichment Analysis
#' 
#' This script performs pathway enrichment analysis on proteomics data using ReactomePA.
#' It converts UniProt accessions to Entrez IDs and performs enrichment analysis
#' for both up- and down-regulated proteins.
#'
#' Required packages:
#' - dplyr
#' - biomaRt
#' - ReactomePA
#' - BiocManager

# Package management ----

#' Install and load required packages
#' @param bioc_packages Vector of Bioconductor package names
#' @param cran_packages Vector of CRAN package names
install_and_load_packages <- function(bioc_packages, cran_packages) {
  # Install BiocManager if not present
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install and load CRAN packages
  for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # Install and load Bioconductor packages
  for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# Define required packages
bioc_packages <- c(
  "biomaRt",
  "ReactomePA",
  "enrichplot",
  "org.Hs.eg.db"
)
cran_packages <- c("dplyr")

# Install and load all required packages
install_and_load_packages(bioc_packages, cran_packages)

#' Convert UniProt accessions to Entrez IDs
#'
#' @param uniprot_acc Vector of UniProt accession numbers
#' @return Vector of corresponding Entrez IDs
convert_uniprot_to_entrez <- function(uniprot_acc) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  entrez_ids <- getBM(
    attributes = c('uniprotswissprot', 'entrezgene_id'),
    filters = 'uniprotswissprot',
    values = uniprot_acc,
    mart = ensembl
  )
  
  return(as.character(entrez_ids[["entrezgene_id"]]))
}

#' Perform pathway enrichment analysis
#'
#' @param gene_list Vector of Entrez IDs
#' @param background Vector of background Entrez IDs
#' @param p_cutoff P-value cutoff for enrichment
#' @return ReactomePA enrichment results
perform_enrichment <- function(gene_list, background, p_cutoff = 0.05) {
  enrichment <- enrichPathway(
    gene_list,
    pvalueCutoff = p_cutoff,
    universe = background
  )
  return(enrichment)
}

#' Main analysis workflow
#'
#' @param data_file Path to input CSV file
#' @param fc_column Column number for fold change
#' @param p_value_column Column number for p-value
#' @param uniprot_column Column number for UniProt accessions
#' @param output_prefix Prefix for output files
main <- function(data_file, 
                fc_up_col = "log2_fold_change",
                fc_down_col = "log2_fold_change",
                p_value_col = "p_value",
                uniprot_col = "uniprot_id",
                fc_up_cutoff = 1,    # log2(2)
                fc_down_cutoff = -1,  # log2(0.5)
                p_value_cutoff = 0.05,
                output_prefix = "enrichment") {
  
  # Read data
  data <- read.csv(data_file, sep = ";")
  
  # Check if required columns exist
  required_cols <- c(fc_up_col, p_value_col, uniprot_col)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), 
         "\nAvailable columns: ", paste(colnames(data), collapse = ", "))
  }
  
  # Extract regulated proteins
  uniprot_acc_up <- data %>%
    filter(.[[fc_up_col]] >= fc_up_cutoff, 
           .[[p_value_col]] <= p_value_cutoff) %>%
    pull(uniprot_col)
  
  uniprot_acc_down <- data %>%
    filter(.[[fc_down_col]] <= fc_down_cutoff, 
           .[[p_value_col]] <= p_value_cutoff) %>%
    pull(uniprot_col)
  
  # Generate background from all identified proteins
  background_uniprot <- unique(data[, uniprot_column])
  
  # Convert all accessions to Entrez IDs
  list_up <- convert_uniprot_to_entrez(uniprot_acc_up)
  list_down <- convert_uniprot_to_entrez(uniprot_acc_down)
  backgroundlist <- convert_uniprot_to_entrez(background_uniprot)
  
  # Perform enrichment analysis
  enrich_up <- perform_enrichment(list_up, backgroundlist)
  enrich_down <- perform_enrichment(list_down, backgroundlist)
  
  # Generate plots
  dotplot_up <- dotplot(enrich_up, showCategory = 20, font.size = 8)
  dotplot_down <- dotplot(enrich_down, showCategory = 20, font.size = 8)
  
  # Save results
  summary_up <- as.data.frame(summary(enrich_up))
  summary_down <- as.data.frame(summary(enrich_down))
  
  write.csv(summary_up, 
            file = paste0(output_prefix, "_upregulated_pathways.csv"))
  write.csv(summary_down, 
            file = paste0(output_prefix, "_downregulated_pathways.csv"))
  
  # Return results as a list
  return(list(
    enrichment_up = enrich_up,
    enrichment_down = enrich_down,
    plots = list(up = dotplot_up, down = dotplot_down)
  ))
}

# Example usage:
if (sys.nframe() == 0) {
  results <- main(
    data_file = "data2.csv",
    output_prefix = "enrichment_results"
  )
}
