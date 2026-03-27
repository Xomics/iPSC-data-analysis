<<<<<<< HEAD
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

## Nextflow workflow
Run the analysis workflow using Nextflow


Example run command:
```
nextflow run main.nf 
	--mae_object path/to/h5MAE/ 
	--output output/dir 
	--r_config bin/Analysis/config.yml 
	--isa_file path/to/20240822_Xomics_iPSC_Neurons_new_READ-ONLY.xlsx 
	--markers_dir path/to/Markers_files/
	--alt_omics_dir path/to/iPSC_MAE/
```



Example run command in Digital Research Environment (DRE):
* Need to specify dre.config and container_dir wit .sif files
```
nextflow run main.nf 
	--mae_object /mnt/data/data/iPSC_data/iPSC_MAE/h5MAE/ 
	-c dre.config 
	--container_dir /mnt/data/software/singularity/ 
	--output /output/dir 
	--r_config bin/Analysis/config.yml 
	--isa_file /mnt/data/data/iPSC_data/iPSC_MAE/20240822_Xomics_iPSC_Neurons_new_READ-ONLY.xlsx 
	--markers_dir /mnt/data/data/iPSC_data/Markers_files/
	--meth_dir /mnt/data/data/iPSC_data/Meth_NF_files/
	--alt_omics_dir /mnt/data/data/iPSC_data/iPSC_MAE/
```

## Licenses

=======
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

## Nextflow workflow
Run the analysis workflow using Nextflow


Example run command:
```
nextflow run main.nf 
	--mae_object path/to/h5MAE/ 
	--output output/dir 
	--r_config bin/Analysis/config.yml 
	--isa_file path/to/20240822_Xomics_iPSC_Neurons_new_READ-ONLY.xlsx 
	--markers_dir path/to/Markers_files/
	--alt_omics_dir path/to/iPSC_MAE/
```



Example run command in Digital Research Environment (DRE):
* Need to specify dre.config and container_dir wit .sif files
```
nextflow run main.nf 
	--mae_object /mnt/data/data/iPSC_data/iPSC_MAE/h5MAE/ 
	-c dre.config 
	--container_dir /mnt/data/software/singularity/ 
	--output /output/dir 
	--r_config bin/Analysis/config.yml 
	--isa_file /mnt/data/data/iPSC_data/iPSC_MAE/20240822_Xomics_iPSC_Neurons_new_READ-ONLY.xlsx 
	--markers_dir /mnt/data/data/iPSC_data/Markers_files/
	--meth_dir /mnt/data/data/iPSC_data/Meth_NF_files/
	--alt_omics_dir /mnt/data/data/iPSC_data/iPSC_MAE/
```

## Licenses

>>>>>>> 10b2857eac83a420e0723015000057fafaa163fb
