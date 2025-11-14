# Single-Cell RNA-seq Analysis: iPSC Heterogeneity

> **A comprehensive end-to-end scRNA-seq analysis pipeline demonstrating best practices in quality control, batch correction, cell type annotation, and differential expression analysis**


## Overview

This repository provides a **complete, production-ready workflow** for analysing single-cell RNA-sequencing (scRNA-seq) data, from raw counts through biological interpretation. Using iPSC data from three individuals (Tung et al., 2017), this pipeline demonstrates:

- **Rigorous quality control** with pre-normalisation doublet detection
- **Multiple normalization strategies** with method comparison
- **Batch effect correction** while preserving biological variation
- **Unsupervised clustering** and cell type annotation
- **Cell cycle analysis** with canonical markers
- **Dual differential expression strategy** for robust testing
- **Pathway enrichment analysis** with biological interpretation
- **Publication-quality visualizations** throughout

### Why This Pipeline?

Unlike basic tutorials, this pipeline:
- Goes beyond preprocessing to complete downstream analysis
- Implements quantitative validation of batch correction
- Uses pseudo-bulk aggregation for robust differential expression
- Includes both cluster-level and cell-type-level comparisons
- Provides detailed biological interpretation at every step
- Follows current best practices (Lun et al., Bioconductor OSCA)


## Dataset

**Source**: [Tung et al. (2017)](https://doi.org/10.1038/srep39921)

**Experimental Design**:
- **Biological Question**: How much heterogeneity exists within iPSC populations?
- **Samples**: 3 Yoruba individuals (NA19098, NA19101, NA19239)
- **Cell Type**: Induced pluripotent stem cells (iPSCs)
- **Replicates**: 3 technical replicates per individual
- **Cells**: 864 cells = 650 high-quality cells after QC
- **Genes**: 14,677 genes = 13,839 analyzed

**Key Feature**: Nested experimental design (replicates within individuals) demonstrates handling of confounded batch effects — a common challenge in real-world scRNA-seq studies.


## Analysis Pipeline

### Section 1: Quality Control & Preprocessing
- Initial cell/gene filtering based on QC metrics
- **Pre-normalisation doublet detection** using scDblFinder
- Removal of low-quality cells and doublets (~3% doublet rate)
- Visualisation of QC distributions

**Key Innovation**: Doublets removed *before* normalisation to prevent artificial clusters

### Section 2: Normalisation
- Comparison of three methods:
  - Raw log transformation (baseline)
  - CPM normalisation  
  - **Scran deconvolution** (selected for downstream analysis)
- Side-by-side PCA comparison showing method effectiveness
- Scran uses pooling-based size factors to handle composition biases


### Section 3: Batch Correction
- Highly Variable Gene (HVG) selection: Top 2000 genes
- **fastMNN** correction for technical replicates
- Quantitative validation using R² metrics:
  - Replicate R2 reduced: 0.69 → 0.07 (10-fold reduction)
  - Individual R2 preserved: 0.56 → 0.52 (biological signal retained)
- UMAP visualisation before/after correction

**Critical Design Decision**: Correct for `replicate` (technical) but preserve `individual` (biological) variation

### Section 4: Clustering & Cell Type Annotation
- Graph-based clustering using Shared Nearest Neighbors (SNN)
- Walktrap community detection algorithm
- 3 main clusters identified, each corresponding to one individual
- **SingleR** automatic cell type annotation
  - Reference: Human Primary Cell Atlas
  - Result: All cells confirmed as iPSC/ESC (no contamination)
- Marker gene discovery with AUC-based ranking
- Heatmap of top cluster-specific markers

**Finding**: Clusters represent individual genetic variation, not cell types

### Section 5: Cell Cycle Analysis
- Gene sets: Canonical S/G2M markers from Macosko et al. (2015)
  - S phase: 42 genes (PCNA, MCM5, CDC6, etc.)
  - G2M phase: 54 genes (MKI67, TOP2A, CDK1, etc.)
- Per-cell scoring and phase classification (G1, S, G2M)
- Distribution analysis: 49% G1, 28% S, 23% G2M
- UMAP overlays showing cell cycle scores

**Key Result**: Cell cycle heterogeneity exists but doesn't drive main clustering structure

### Section 6: Differential Expression Analysis

**Two Complementary Strategies**:

#### Strategy 1: Cluster Comparisons (Donor-Level)
- 3 pairwise comparisons between individuals
- Pseudo-bulk aggregation by donor × replicate × cluster
- EdgeR quasi-likelihood framework
- Blocking for technical replicates in model design

**Results**:
- NA19101 vs NA19098: 161 DE genes
- NA19101 vs NA19239: 437 DE genes (most divergent)
- NA19098 vs NA19239: 144 DE genes

#### Strategy 2: Cell-Type Comparisons (Within-Type)
- Tests donor differences within epithelial cell population
- Controls for cell type composition
- Same statistical framework as Strategy 1

**Statistical Features**:
- FDR < 0.05 significance threshold
- Model: `~ 0 + group + replicate`
- Variance-stabilised normalisation
- Robust to dropouts via pseudo-bulk approach

### Section 7: Pathway Analysis

**Gene Set Enrichment Analysis (GSEA)**:
- Database: Hallmark gene sets (MSigDB)
- Method: fgsea with 10,000 permutations
- Ranking: Signed -log10(p-value) × sign(logFC)
- Output: Enriched pathways per comparison (FDR < 0.1)

**Per-Cell Pathway Scoring**:
- Method: singscore for single-cell resolution
- Pathways: Selected from GSEA-significant results
- Visualisation: UMAP overlays of pathway activity
- Interpretation: Pathway differences across individuals

**Common Enriched Pathways**:
- Cell cycle checkpoints (G2M, E2F targets)
- Metabolic pathways
- DNA repair mechanisms
- Stem cell maintenance signatures


## Repository Structure

```
scrnaseq-ipsc-analysis/
├── README.md                              # This file
├── .gitignore                             # Git exclusions
│
├── Complete_scRNAseq_Analysis.Rmd         # Main analysis pipeline (49KB)
├── Pre-processing-scRNA-Seq-data.Rmd      # Preprocessing of scRNA seq pipeline
│
├── data/                                  # Input data folder
│   ├── molecules.txt                      # Raw count matrix (genes × cells)
│   ├── annotation.txt                     # Cell metadata
│   └── umi.rds                            # Preprocessed SingleCellExperiment
│
└── output/                                # Generated results
    ├── doublet_detection.png
    ├── normalisation_comparison.png
    ├── batch_correction_umap.png
    ├── clustering_overview.png
    ├── cell_cycle_analysis.png
    ├── markers/
    │   ├── markers_cluster_1.csv
    │   ├── markers_cluster_2.csv
    │   └── markers_cluster_3.csv
    └── DE_results/
        ├── DE_clusters_Clusters.csv
        ├── volcano_Clusters.png
        

## Key Results

### Cell Population Structure

| Cluster | Cells | Individual | Cell Cycle Distribution |
|---------|-------|------------|-------------------------|
| 1 | 222 | NA19101 | 48% G1, 29% S, 23% G2M |
| 2 | 167 | NA19098 | 51% G1, 27% S, 22% G2M |
| 3 | 261 | NA19239 | 49% G1, 27% S, 24% G2M |

**Key Findings**:
- Each cluster = one individual (strong genetic signal)
- Cell cycle similar across donors (consistent proliferation)
- No differentiated cell contamination (all iPSCs confirmed)
- Batch effects successfully removed without over-correction

### Differential Expression Summary

**Inter-Individual Variation**:
- Total unique DE genes: 742 across all comparisons
- Most divergent pair: NA19101 vs NA19239 (437 genes)
- Most similar pair: NA19098 vs NA19239 (144 genes)

**Top Differentially Expressed Gene**: ENSG00000106153 (logFC > 10 in multiple comparisons)

### Biological Interpretation

1. **Batch Correction Validation** 
   - Technical replicates successfully merged (R2 = 0.07)
   - Biological variation preserved (R2 = 0.52)
   - Quantitative metrics confirm appropriate correction

2. **Individual Genetic Variation**
   - Strong donor-specific expression signatures
   - Reproducible across technical replicates
   - Affects hundreds of genes

3. **Cell Cycle Heterogeneity**
   - Asynchronous proliferation (expected for iPSCs)
   - Not driving main clustering structure
   - Similar distribution across individuals

4. **Pluripotency Confirmed**
   - All cells annotated as iPSC/ESC
   - No evidence of differentiation
   - Consistent stem cell signatures


## Methodological Highlights

### What Makes This Analysis Robust?

1. **Pre-normalisation doublet detection**
   - Prevents doublets from forming artificial clusters
   - Identified 20 doublets (~3%) before normalisation
   - Uses scDblFinder expectation-maximisation algorithm

2. **Quantitative batch correction validation**
   - Not just visual inspection
   - R2 metrics quantify technical vs biological variation
   - Confirms appropriate correction level

3. **Dual differential expression strategy**
   - Strategy 1: Broad donor comparisons
   - Strategy 2: Refined within-cell-type comparisons
   - Complementary insights from both approaches

4. **Pseudo-bulk aggregation**
   - Robust to dropout and zero-inflation
   - More statistical power than single-cell DE
   - Follows current best practices (Squair et al. 2021)

5. **Proper handling of confounded design**
   - Individual × replicate nested structure
   - Correct for technical batch, preserve biology
   - Documented decision-making process



## References

### Primary Dataset
Tung, P.-Y. et al. (2017). Batch effects and the effective design of single-cell gene expression studies. *Scientific Reports*, 7, 39921. [https://doi.org/10.1038/srep39921](https://doi.org/10.1038/srep39921)

### Key Methodological Papers

**Quality Control & Normalization**
- Lun, A. T. L. et al. (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data. *F1000Research*, 5, 2122.
- Lun, A. T. L. et al. (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. *Genome Biology*, 17, 75.

**Batch Correction**
- Haghverdi, L. et al. (2018). Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nature Biotechnology*, 36, 421–427.

**Doublet Detection**
- Germain, P.-L. et al. (2022). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research*, 10, 979.

**Differential Expression**
- Robinson, M. D. et al. (2010). edgeR: a Bioconductor package for differential expression analysis. *Bioinformatics*, 26(1), 139-140.
- Squair, J. W. et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12, 5692.

**General scRNA-seq**
- Amezquita, R. A. et al. (2020). Orchestrating single-cell analysis with Bioconductor. *Nature Methods*, 17, 137–145.

