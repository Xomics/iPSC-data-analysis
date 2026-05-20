# iPSC Data Analysis

## Description
This project is part of the Dutch X-omics initiative and focuses on the analysis of multi-omics data from several iPSC lines. The data included genomics, transcriptomics, proteomics, and several different types of metabolomics.

## Goal
The goal is to identify disease-specific traits of the patient-derived iPSC compared to various control lines.

## Pipeline overview


## File outline
|- [README.md](README.md) &emsp;&emsp;&emsp;&emsp;# This readme\
|- [Data analysis plan.md](Data_analysis_plan.md) # Description of the data analysis plan (DAP)\
|- [DAP.png](DAP.png) &emsp;&emsp;&emsp;&emsp;&emsp;&ensp; # Overview image of step 4 and 5 of the DAP

## Software and packages used
- R version 4.3.2 in the <a href= "file:///Z:/software/R"> software folder </a>
(Z:\software\R)


## Data sources
- <a href= "file:///Z:/data/iPSC_data/iPSC_CNV-WES/README.md"> CNV WES </a> (Z:\data\iPSC_data\iPSC_CNV-WES\README.md)
- <a href= "file:///Z:/data/iPSC_data/iPSC_RNAseq/README.md"> RNA-seq </a> (Z:\data\iPSC_data\iPSC_RNAseq\README.md)
- <a href= "file:///Z:/data/iPSC_data/iPSC_proteomics/README.md"> Proteomics </a> (Z:\data\iPSC_data\iPSC_proteomics\README.md)
- <a href= "file:///Z:/data/iPSC_data/iPSC_metabolomics_TML/README.md"> Metabolomics (RUMC, Nijmegen) </a> (Z:\data\iPSC_data\iPSC_metabolomics_TML\README.md)
- <a href= "file:///Z:/data/iPSC_data/iPSC_metabolomics_lipidomics_LUMC/README.md"> Metabolomics and Lipidomics (LUMC, Leiden) </a> (Z:\data\iPSC_data\iPSC_metabolomics_lipidomics_LUMC\README.md)


## Digital Reserarch Environment:

#### Support MyDRE:
* https://support.mydre.org/

#### Gitlab-DRE connection errors
If Gitlab gives a timeout warning when used from within the DRE:

1) Confirm if https://gitlab.cmbi.umcn.nl/x-omics_ipscs/ipsc-data-analysis.git/ can be accessed through the browser
2) Correct the proxy configuration with the following git command:

    ````git config --global http.proxy http://proxy.mydre.org:3128````

#### CRAN/Bioconductor
* Configuration instructions https://support.mydre.org/portal/en/kb/articles/proxy-configurations-rstudio#Introduction

Use a Biconductor mirror:
* Add `options(BioC_mirror = "https://ftp.gwdg.de/pub/misc/bioconductor/")` to the `C:\Program Files\R\R-[your R version]\etc\Rprofile.site`(requires DRE admin rights)

## Nextflow Workflow

Run the analysis workflow using Nextflow.

### Usage
```bash
nextflow run main.nf \
    -profile singularity \        # Use `-profile docker` for Docker
    --mae_object <path/to/h5MAE/> \
    --output <output/dir> \
    --r_config <bin/Analysis/config.yml> \
    --isa_file <path/to/ISA_file>
```

**Optional parameters:**

| Parameter | Description |
|---|---|
| `--markers_dir` | Path to markers files directory |
| `--alt_omics_dir` | Path to alternative omics directory |

---

### Running in a Digital Research Environment (DRE)

Use the `dre` profile and provide a `--container_dir` pointing to your Singularity (`.sif`) files.
```bash
nextflow run main.nf \
    -profile dre \
    --container_dir /mnt/data/software/singularity/ \
    --mae_object /mnt/data/data/iPSC_data/iPSC_MAE/h5MAE/ \
    --output /output/dir \
    --r_config bin/Analysis/config.yml \
    --isa_file /mnt/data/data/iPSC_data/iPSC_MAE/20240822_Xomics_iPSC_Neurons_new_READ-ONLY.xlsx
```

**Optional parameters:**

| Parameter | Description |
|---|---|
| `--markers_dir` | `/mnt/data/data/iPSC_data/Markers_files/` |
| `--meth_dir` | `/mnt/data/data/iPSC_data/Meth_NF_files/` |
| `--alt_omics_dir` | `/mnt/data/data/iPSC_data/iPSC_MAE/` |
## Data files on Zenodo
To run the Nextflow pipeline, input files can be downloaded on Zenodo:
[Link_to_Zenodo](https://doi.org/10.5281/zenodo.19129856)

## Licenses

