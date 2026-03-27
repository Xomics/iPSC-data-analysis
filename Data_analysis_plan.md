# Data analysis plan (DAP)

**Table of contents**
- [Data analysis plan (DAP)](#data-analysis-plan-dap)
  - [Data stewardship](#data-stewardship)
  - [Outline](#outline)
  - [Data types and sources](#data-types-and-sources)
  - [Data Pre-processing](#data-pre-processing)
  - [Quality Assurance](#quality-assurance)
    - [Genomics](#genomics)
    - [Transcriptomics](#transcriptomics)
    - [Epigenomics methylation](#epigenomics-methylation)
    - [Untargeted proteomics](#untargeted-proteomics)
    - [Untargeted Metabolomics (Nijmegen)](#untargeted-metabolomics-nijmegen)
    - [Untargeted Metabolomics (Leiden)](#untargeted-metabolomics-leiden)
    - [Untargeted Lipidomics (Leiden)](#untargeted-lipidomics-leiden)
    - [Targeted Lipidomics (Leiden)](#targeted-lipidomics-leiden)
  - [Data integration and analysis](#data-integration-and-analysis)
    - [Single Omics Analysis](#single-omics-analysis)
    - [Multi-Omics Analysis](#multi-omics-analysis)
  - [FAIR Workflow](#fair-workflow)

## Data stewardship
- **Raw Data Storage:** DRE X-omics will be stored in blob storage, with methylation data occupying a significant space (>10 TB). 
- **Metadata Management:** Metadata will be organized using ISA templates through the FAIR Data Station. Data providers will fill in ISA files, and pre-processing steps will be clearly described in the Assay files. 
- **Phenotypic Information:** Assessment of available phenotypic information to determine the extent of data provided. 
- **Data Organization:** Multi-omics data will be structured into one R's MultiAssayExperiment object.

## Outline

The general outline of the data analysis plan includes:

1. Data FAIRification
2. Pre-processing of the different omics
3. QC per omic
4. Single omic analysis
5. Multi-omics analysis

A brief overview of the data analysis plan for step 4 and 5 is shown below.

![alt text](DAP.png)

## Data types and sources

There are 7 seven different data types across 8 different data sources. For each data type the data provider is indicated.

| Data type | Details | Provider (department, institute, city) |
|-|-|-|
|Genomics | Low coverage WGS | Stem Cell Facility, Radboudumc Technology Center for Stem Cells, Nijmegen|
|Transcriptomics| RNA-seq |Utrecht Sequencing Facility (USEQ), University Medical Center Utrecht, Utrecht|
|Epigenomics |methylation-seq | SNP&SEQ Technology Platform, National Genomics Infrastructure (NGI), Science for Life Laboratory, Uppsala University, Uppsala|
|Proteomics |Untargeted |Translational Metabolic Laboratory, Radboudumc Technology Center for Mass Spectrometry, Nijmegen|
|Metabolomics |Untargeted |Translational Metabolic Laboratory, Radboudumc Technology Center for Mass Spectrometry, Nijmegen|
|Metabolomics |Untargeted|Center for Proteomics and Metabolomics, Leiden University Medical Center, Leiden|
|Lipidomics |Untargeted|Center for Proteomics and Metabolomics, Leiden University Medical Center, Leiden|
|Lipidomics |Targeted|Center for Proteomics and Metabolomics, Leiden University Medical Center, Leiden|

## Data Pre-processing 

- **Confirm data formats** per omic based on information from the provider
- **Batch Adjustment and Transformation**: Data will undergo batch adjustment and transformation to enhance comparability, based on the delivery format. 
- **Noise Reduction**: Consider implementing the Variational Autoencoder (VAE) approach for noise reduction. 
- **Quality Assurance**: Check for possible sample mix-ups by comparing for example RNA and WES for individual samples. 

## Quality Assurance
The general QA steps include:
- Descriptive statistics
- Missing data reports
- Outlier detection
- Data distribution analysis
- Data normalization
- Batch effect analysis and correction
- Check for sample mislabelling and incorrect samples

Below the additional QC steps per omic are described.

### Genomics
- Sex matching
- Expected genetic defect observed
- Sequencing depth
- read mapping rate

### Transcriptomics
- Sex matching
- Sample matching RNA-DNA
- Expected genetic defect observed
- Sequencing depth
- Read mapping rate

### Epigenomics methylation
- Sex matching
- Sample matching RNA-DNA
- Sequencing depth

### Untargeted proteomics
- RT time variability
- Mass accuracy
- 

### Untargeted Metabolomics (Nijmegen)

### Untargeted Metabolomics (Leiden)

### Untargeted Lipidomics (Leiden)

### Targeted Lipidomics (Leiden)

## Data integration and analysis

### Single Omics Analysis 

- **Unsupervised Analyses**:  Use Principal Component Analysis (PCA) to find main variation within omics layers, identify technical variations, and associate with phenotypic characteristics. 

- **Hierarchical clustering heatmap**: (Euclidallucian/ k-means) 

- **Variability Assessment**: Identify omics features with variable expression across different cell lines, within the same cell lines, and evaluate stability using Coefficient of Variation. 

- **Supervised Analyses**: Investigate linear associations of single features with phenotypic variables and disease states. Identify potential biomarkers for individual diseases and disease severity, as well as shared biomarkers across diseases. 

- **ANOVA and globalANCOVA**: Employ Analysis of Variance (ANOVA) and globalANCOVA to assess contributions of factors like clone, patient, and disease state to omics readouts, and to discern disease-specific associations. 

 

### Multi-Omics Analysis 

- **Multi-Omics Factor Analysis (MOFA)**: Apply MOFA to discover latent dimensions in multi-omics layers and identify biomarkers for different diseases. 

- **Pathway Enrichment Analysis**: Investigate overlap in disease pathways across disorders using pathway enrichment analyses.
  - Moreover, the MOFA factor weights can be employed for pathway enrichment analysis, potentially identifying more biological background to the MOFA factors. 

- **Multilayer network analyses**, including multilayered diffusion propagation methods (random walk, heat diffusion) that make combined use of established knowledge networks (for example containing the interactions between genes, proteins and metabolites) and the multi-omics data obtained       

- **GO-Term enrichment** to identify the most functional terms
  - With **DOMINO active module identification** to determine the most relevant modules per disease background

- **DIABLO**: Employ DIABLO to identify predictive latent variables in multi-omics data for categorical outcome variables (disease), aiding in the discovery of new multi-omics biomarkers. 

## FAIR Workflow 
Accesibility, eproducability, and reusability: Ensure reproducibility by sharing the workflow on Workflowhub.eu along with processed data in R/Python object formats. 

 