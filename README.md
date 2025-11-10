Welcome to the Pre-processing scRNA-Seq Data repository! This repository contains a step-by-step R workflow for cleaning and preparing single-cell RNA sequencing data (scRNA-seq) prior to downstream analyses such as clustering, dimensionality reduction, and differential expression.

# Overview
Single-cell RNA-seq (scRNA-seq) provides powerful, high-resolution views of gene expression at the level of individual cells. However, it also introduces unique challenges in data preparation and quality control. This workflow demonstrates how to:

1. Load Raw Count Matrices and associated sample annotations.
2. Create a SingleCellExperiment object that encapsulates count data, annotations, and metadata in a convenient structure.
3. Perform Quality Control (QC) by removing low-quality cells and uninformative genes.
4. Normalise and Transform the raw counts (e.g., by applying log transformations).
5. Visualise and Examine QC Metrics, including outliers and potential batch effects.
6. Identify Mitochondrial and Spike-In Reads for additional technical QC.
7. Save the Cleaned/Annotated Object for downstream analyses (e.g., clustering and differential expression).

# Repository Structure
- **scRNA_preprocessing.Rmd**  
  *RMarkdown document with the full walkthrough*

- **data/**  
  - `molecules.txt`  
    *Example count matrix (genes x cells)*  
  - `annotation.txt`  
    *Cell-level annotation/metadata*  
  - `umi.rds`  
    *Saved SingleCellExperiment object after QC*

# Quick Start
1. Clone or Download this repository.
2. Install R and R packages required for scRNA-seq data processing. Major dependencies:
   - SingleCellExperiment
   - scater
   - AnnotationDbi, org.Hs.eg.db, EnsDb.Hsapiens.v86
   - tidyverse, ggplot2, dplyr
   - readr, readxl
3. Open `scRNA_preprocessing.Rmd` in RStudio (or your preferred environment).
4. Update file paths in the code if necessary (e.g., to point to your own local directory where the data are stored).
5. Run the code chunk by chunk or knit the R Markdown to produce an HTML report of the process.

# Data Description
The example dataset used in the workflow is from Tung et al. 2017, where induced pluripotent stem cells (iPSCs) were generated from three different individuals.
The data were produced on the Fluidigm C1 platform and include:
- Raw count (UMI) matrix (`molecules.txt`), rows = genes, columns = cells.
- Cell annotation (`annotation.txt`), which includes metadata such as cell IDs, individual, replicate, etc.
- ERCC spike-ins and mitochondrial genes for QC.

# References
Tung, P.-Y., Blischak, J. D., Hsiao, C. J., Knowles, D. A., Burnett, J. E., Pritchard, J. K., & Gilad, Y. (2017). Batch effects and the effective design of single-cell gene expression studies. *Scientific Reports*, 7, 39921.
