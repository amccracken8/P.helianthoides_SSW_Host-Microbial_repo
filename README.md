# Precursors of Sea Star Wasting: Immune and Microbial Disruption During Initial Disease Outbreak in Southeast Alaska

This repository contains code and data processing scripts for the study:

**"Precursors of Sea Star Wasting: Immune and Microbial Disruption During Initial Disease Outbreak in Southeast Alaska"**

DOI: TBA

The analysis compares wild-sampled individuals from outbreak-free sites (“Naive”) to those from sites where wasting disease is present but individuals appear visually healthy (“Exposed”). Scripts include:

- Data cleaning
- Statistical analyses
- Figure generation

---

## Related Repositories and Data

- **Dryad DOI:** TBA
- **Manuscript DOI:** TBA

## Folders and Code Overview

### **`Cleaning_Mapping_Reads/`**

*Cleaning and mapping paired-end reads and transcript quantification*

- **`star_loop.sh`** – Bash script to run gene alignment and transcript quantification on paired FASTA files to a reference.
- **`star_counts_matrix.R`** – R script to combine outputs from STAR gene mapping into a single counts matrix.

---

### **`DEG/`**

*Differential gene expression analysis and plotting*

- **`Annotations/`** – Annotation of mapped transcripts:
    - **`diamond.sh`** – Bash script for running DIAMOND blastX.
    - **`subset_fasta.py`** – Python script to subset FASTA files for differentially expressed genes.
- **`DESeq2_pycnoRNA_STAR_NE.R`** – Differential expression analysis of transcript counts.
- **`heatmaps_by_catagory/`** – Heatmaps of differential gene expression.

---

### **`16S_Counts/`**

*Microbial classification and quantification*

- **`Qiime2_pipeline_Silva.txt`** – QIIME2 pipeline for microbial classification.
- **`pycno_samples_NE.txt`** – Sample metadata.
- **`level-7.csv`** – Counts matrix of classified microbial taxa per sample.
- **`NE_ANCOM-BC-Lev7.R`** – Analysis of differentially abundant microbes.
- **`deseqPCA_microbes.R`** – PCA analysis of microbial composition.
- **`RelAbundancePlot.R`** – Microbial abundance visualizations.

---

### **`GO/`**

*Gene ontology enrichment analysis*

- **`TOP_GO.R`** – Gene ontology analysis of differentially expressed genes.
- **`topgo_wgcna_mods.R`** – Gene ontology analysis of genes clustered in WGCNA modules “red” and “purple”.

---

### **`Gene_Microbe_WGCNA/`**

*Weighted gene-microbe correlation network analysis and plotting*

- **`WGCNA.R`** – WGCNA analysis of gene expression and microbial abundance counts.
- **`WGCNA_expression_abund_plots_nf.R`** – Abundance and expression plots of module membership.
- **`new.netolot.R`** – Network visualization of WGCNA modules.
