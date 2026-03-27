library(GenomicRanges)

extract_beta_values_in_DMRs <- function(meth, DMR_list) {
  
  feature_metadata <- elementMetadata(meth)
  feature_metadata <- as.data.frame(feature_metadata)
  rownames(feature_metadata) <- paste0(feature_metadata$chr, '_', feature_metadata$start)
    
  ## Measurement data
  beta_values_df <- assays(meth)[[1]]
  rownames(beta_values_df) <- rownames(feature_metadata)
  
  ## GRanges objects
  regions_gr <- GRanges(
    seqnames = DMR_list$chr,
    ranges = IRanges(start = DMR_list$start, end = DMR_list$end)
  )
  
  coords_gr <- GRanges(
    seqnames = feature_metadata$chr,
    ranges = IRanges(start = feature_metadata$start, end = feature_metadata$start)
  )
  
  overlaps <- findOverlaps(coords_gr, regions_gr)
  
  result <- data.frame(
    coord_idx = queryHits(overlaps),
    region_idx = subjectHits(overlaps),
    chr = as.character(seqnames(coords_gr[queryHits(overlaps)])),
    pos = start(coords_gr[queryHits(overlaps)]),
    region_start = start(regions_gr[subjectHits(overlaps)]),
    region_end = end(regions_gr[subjectHits(overlaps)])
  )
  
  feature_metadata$coord_idx <- 1:nrow(feature_metadata)
  feature_metadata$CpG <- rownames(feature_metadata)
  result <- merge(result, feature_metadata, by='coord_idx')
  
  beta_values_df <- as.data.frame(beta_values_df)
  beta_values_df_sub <- beta_values_df[result$CpG,]
  return(beta_values_df_sub)
  
}

