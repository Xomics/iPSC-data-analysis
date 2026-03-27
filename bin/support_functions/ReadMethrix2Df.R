## Function to read Methrix object into measurement data and feature metadata
##
read_meth2df <- function(meth, isa_assay) {
  
  ## Feature metadata
  feature_metadata <- elementMetadata(meth)
  feature_metadata <- as.data.frame(feature_metadata)
  rownames(feature_metadata) <- paste0(feature_metadata$chr, '_', feature_metadata$start)
  
  ## Measurement data
  beta_values_df <- assays(meth)[[1]]
  
  ## Feature names
  rownames(beta_values_df) <- rownames(feature_metadata)
  
  ## Re-format sample names
  new_samples_names <- gsub( "WI-3791-", "", colnames(beta_values_df))
  new_samples_names <- gsub( "Epigenomics", "Proteomics", new_samples_names)
  new_samples_names <- gsub( "_1_val_1_bismark_bt2_pe", "", new_samples_names)
  new_samples_names <- gsub( "-", "_", new_samples_names)
  
  ## Get universal names by mapping to ISA file
  isa_assay <- isa_assay[!duplicated(isa_assay$`assay identifier`),]
  universal_names_meth <- lapply(new_samples_names, function(x) {
    name_uni <- isa_assay[isa_assay$`assay identifier` == x, "sample identifier"]
    bio_rep <- isa_assay[isa_assay$`assay identifier` == x, "biological replicate"]
    name_uni <- paste0(name_uni, "_", bio_rep)
    name_uni
  })
  beta_values_df <- data.frame(beta_values_df)
  colnames(beta_values_df) <- as.character(universal_names_meth)
  colnames(beta_values_df) <- gsub( "-", "_", colnames(beta_values_df))
  colnames(beta_values_df) <- gsub( "409B", "X409B", colnames(beta_values_df))
  
  return(list(beta_values_df, feature_metadata)) 
}



## Function to read Methrix object into measurement data and feature metadata
##
read_meth2df_OLD <- function(meth, isa_assay) {
  
  ## Feature metadata
  feature_metadata <- elementMetadata(meth)
  feature_metadata <- as.data.frame(feature_metadata)
  rownames(feature_metadata) <- paste0(feature_metadata$chr, '_', feature_metadata$start)
  
  ## Measurement data
  beta_values_df <- assays(meth)[[1]]
  
  ## Feature names
  rownames(beta_values_df) <- rownames(feature_metadata)
  
  ## Re-format sample names
  new_samples_names <- gsub( "WI-3791-", "", colnames(beta_values_df))
  new_samples_names <- gsub( "Epigenomics", "Proteomics", new_samples_names)
  new_samples_names <- gsub( "_1_val_1_bismark_bt2_pe", "", new_samples_names)
  new_samples_names <- gsub( "-", "_", new_samples_names)
  
  ## Get universal names by mapping to ISA file
  isa_assay <- isa_assay[!duplicated(isa_assay$`assay identifier2`),]
  universal_names_meth <- lapply(new_samples_names, function(x) {
    name_uni <- isa_assay[isa_assay$`assay identifier2` == x, "sample identifier"]
    bio_rep <- isa_assay[isa_assay$`assay identifier2` == x, "biological replicate"]
    name_uni <- paste0(name_uni, "_", bio_rep)
    name_uni
  })
  beta_values_df <- data.frame(beta_values_df)
  colnames(beta_values_df) <- as.character(universal_names_meth)
  colnames(beta_values_df) <- gsub( "-", "_", colnames(beta_values_df))
  colnames(beta_values_df) <- gsub( "409B", "X409B", colnames(beta_values_df))
  
  return(list(beta_values_df, feature_metadata)) 
}