#### Params replacements of markdown paths ####
params <- c()
params$mrna_deseq_data_tab <- "Z:/data/iPSC_data/iPSC_RNAseq/Data_integration/mrna_dseq_filtered.csv"
params$metadata <- "metadata: C:/Users/Cenna.Doornbos/Documents/ipsc-data-analysis/20240514_Xomics iPSC Neurons_new_READ-ONLY.xlsx"


#### Libraries ####
library(readxl)
library(tidyverse)
library(ggrepel)
library(plotly)
library(ggbiplot)

#### params ####
print(Sys.Date())
print(params)

#### Read data ####
# read in data
data_frame_rna <- read.csv(params$mrna_deseq_data_tab) %>%
  mutate(sample_name = gsub("\\.", "-", X), #stratify names
         sample_name = gsub("X409", "409", sample_name), #stratify names
         sample_name = gsub("CDH2", "CHD2", sample_name), #correct CHD2 mislabeling
         .before = X) %>%
  select(-X)

# read in metadata
metadata_rna <- read_excel(params$metadata, sheet = "Assay - transcriptomics") %>%
  left_join(read_excel(params$metadata, sheet = "Sample - Induced pluripoten "), by = "sample identifier", suffix = c("_rna", "_sam")) %>%
  left_join(read_excel(params$metadata, sheet = "ObservationUnit"), by = "observation unit identifier", suffix = c("", "_OU")) %>%
  mutate("sample_name" = gsub("_Transcriptomics", "", `assay identifier`), .after = "assay identifier")

#### Sex matching ####
#ensemble labels xist (female) and uty (male)
# xist <- "ENSG00000229807"
# uty <- "ENSG00000183878"

# adapt data to generate figure
df <- data_frame_rna %>%
  select(sample_name, ENSG00000229807, ENSG00000183878) %>%
  rename("XIST" = ENSG00000229807,
         "UTY" = ENSG00000183878) %>%
  mutate("measured sex" = case_when(UTY <= 0.1 ~ "female",
                                    UTY > 1.5 ~ "male",
                                    .default = "NA")) %>%
  left_join(metadata_rna, by = "sample_name") %>%
  dplyr::rename("sex" = gender) %>%
  mutate("label" = ifelse(`measured sex` != sex, sample_name, ""))


#generate sex match figure
p <- ggplot(df, aes(XIST, UTY, colour = sex, label = label, text = sample_name, shape = `measured sex`)) +
  # geom_point(size = 4, alpha = 0.4, stroke = NA) +
  geom_jitter(size = 4, alpha = 0.4, stroke = NA) +
  labs(title = "RNA determined sex") +
  theme_minimal()

#save figure with labels
# ggsave(paste0(Sys.Date(), "_RNA determined sex.pdf"), path = file.path(dirname(getwd()), "results"),
#        plot = p + geom_text_repel(max.overlaps = 100, force = 50), width = 7, height = 7)
ggsave(paste0(Sys.Date(), "_RNA determined sex.pdf"), path = file.path(dirname(getwd()), "Documents", "ipsc-data-analysis", "results"),
       plot = p + geom_text_repel(max.overlaps = 100, force = 50), width = 7, height = 7)

#display figure
ggplotly(p, tooltip = c("text", "colour", "shape", "UTY", "XIST"))

#### Additional count matrices ####
# params
params$mrna_deseq <- "Z:/data/iPSC_data/iPSC_RNAseq/iPSC_USEQ_processed/KUL9276_processed/featureCounts/KUL9276_processed_featureCounts_deseq2.txt"
params$mrna_cpm <- "Z:/data/iPSC_data/iPSC_RNAseq/iPSC_USEQ_processed/KUL9276_processed/featureCounts/KUL9276_processed_featureCounts_CPM.txt"
params$mrna_rpkm <- "Z:/data/iPSC_data/iPSC_RNAseq/iPSC_USEQ_processed/KUL9276_processed/featureCounts/KUL9276_processed_featureCounts_RPKM.txt"

# functions
ReadIn <- function(path){
  read.table(path) %>%
    t %>% as.data.frame() %>%
    rownames_to_column(var = "sample_name") %>%
    mutate(sample_name = gsub("\\.", "-", sample_name), #stratify names
           sample_name = gsub("X409", "409", sample_name), #stratify names
           sample_name = gsub("CDH2", "CHD2", sample_name)) #correct CHD2 mislabeling
}

SubSet <- function(df){
  df_deseq_sub <- df %>%  
    select(sample_name, ENSG00000229807, ENSG00000183878) %>%
    rename("XIST" = ENSG00000229807,
           "UTY" = ENSG00000183878)
}

# read in data and select select XIST and UTY subset of the columns for each analysis
df_deseq <- ReadIn(params$mrna_deseq)
df_deseq_sub <- SubSet(df_deseq)

df_cpm <- ReadIn(params$mrna_cpm)
df_cpm_sub <- SubSet(df_cpm)

df_rpkm <- ReadIn(params$mrna_rpkm)
df_rpkm_sub <- SubSet(df_rpkm)

# perform PCA
# source("bin/support_functions/performPCA.R") #doesn't work, see custom function from Casper below

  # function
  performPCA2 <- function(df){
    sampleType <- sub("^(.*?)_.*","\\1", rownames(df))
    colnames(df)[ncol(df)] <- "sampleType"
    
    data.pca <- prcomp(df, center = TRUE, scale.= TRUE)
    g <- ggbiplot(
      data.pca,
      groups = sampleType,
      var.axes = F,
      labels = rownames(df)
    ) + theme_minimal()
    return(g)
  }
  
  #PCA subset
  df_deseq_sub %>%
    column_to_rownames(var = "sample_name") %>%
    performPCA2()
  
  #PCA all
  df_deseq %>%
    column_to_rownames(var = "sample_name") %>%
    select(1:15) %>%
    performPCA2()
