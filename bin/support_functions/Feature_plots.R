library(ggplot2)
library(ggpubr)

boxplot_per_group <- function(df_input, 
                              feature, 
                              assay_name, 
                              disease_group = 'Disease', 
                              y_lab, 
                              title = NULL,
                              phenotype = 'phenotype', 
                              color_group = NULL,
                              comparisons = NULL, 
                              stat_method = "wilcox.test",
			      width = 25, height=17) {

  # Color is phenotype if not otherwise specified
  if (is.null(color_group)) {
	color_group <- phenotype	
  }

  df2plot <- df_input[, c(feature, phenotype, color_group)]
  colnames(df2plot) <- c("y_feature", 'group', 'color_group')
  
  # Color is phenotype if not otherwise specified
  if (!is.null(title)) {
	  title <- title	
  }
  else {
    title <- stringr::str_wrap(paste0(assay_name, "\t", disease_group), width = 60)
  }
  
  
  ggp <- ggplot(df2plot, aes(x = group, y = y_feature, fill = color_group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.9)) +
    geom_jitter(width = 0.2, size =2, alpha=0.6) +
    ggtitle(title) +
    theme_minimal() +
    labs(y = paste0(y_lab), fill = color_group) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
          axis.text.y = element_text(angle = 90))
  
  # Add stat comparisons if specified
  if (!is.null(comparisons)) {
    ggp <- ggp + 
      stat_compare_means(comparisons = comparisons, 
                         method = stat_method, 
                         label = "p.format")  # or "p.format"
  }

  # Clean and save output
  assay_name_clean <- strsplit(assay_name, "\\|")[[1]][1]
  assay_name_clean <- gsub(" ", "", assay_name_clean)
  full_path <- file.path(outdir, disease_group, assay_name_clean)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
  }
  file_path <- file.path(full_path, paste0(feature, '.pdf'))
  ggsave(file_path, ggp, width = width, height = height, units = "cm", dpi = 300)
  
  return(ggp)
}




## Function to prepare data frame for boxlots of top 10 features
prepare_boxplot_df <- function(assay_name) {
  se <- experiments(mae)[assay_name]
  df <- as.data.frame(t(assays(se)[[1]]))
  df <- cbind(df, pheno_df[rownames(df),])
  return(df)
}



## Function to make boxplots per disease group, per omics data, for each significant feature
make_boxplots_signif_feature <- function(assay_name, ylab = 'Log-scaled intensity') {
  
  rlm_res <- MixedModel_complete[MixedModel_complete$Assay == assay_name,]
  head <- head(rlm_res[order(rlm_res$p.value, decreasing = FALSE),], n = 10)
  rlm_res_sign <- rlm_res[rlm_res$significant ==T,]
  df <- prepare_boxplot_df(assay_name)
  
  CHD2_features <- rlm_res_sign[rlm_res_sign$outcome == 'phenotypeDiseased: CHD2',]$feature
  NANS_features <- rlm_res_sign[rlm_res_sign$outcome == 'phenotypeDiseased: NANS-CDG',]$feature
  DM1_features <- rlm_res_sign[rlm_res_sign$outcome == 'DM1_vs_rescue',]$feature
  Diseases_features <- rlm_res_sign[rlm_res_sign$outcome == 'DiseaseDisease',]$feature
  
  #cat('CHD2 vs controls')
  lapply(CHD2_features, function(x) boxplot_per_group(df, x, assay_name, 'CHD2', ylab))
  #cat('NANS vs controls')
  lapply(NANS_features, function(x) boxplot_per_group(df, x, assay_name, 'NANS-CDG', ylab))
  #cat('DM1 vs controls')
  lapply(DM1_features, function(x) boxplot_per_group(df, x, assay_name, 'DM1', ylab))
  lapply(Diseases_features, function(x) boxplot_per_group(df, x, assay_name, 'All diseases', ylab))
}

## Function to make boxplots per disease group, per omics data, for each significant feature
make_boxplots_molecular_markers <- function(markers, assay_name, ylab = 'Log-scaled intensity') {
  
  df <- prepare_boxplot_df(assay_name)
  lapply(markers, function(x) boxplot_per_group(df, x, assay_name, paste0(x), ylab))
}


## Function to get long table of only the marker expression values
## @param stage select Mature of Early
get_long_table_markers <- function(df, omics_id, stage='Mature'){
  markers_df <- markers_df[markers_df$Stage == stage,]
  markers <- na.omit(markers_df[, omics_id])
  marker_genes <- markers_df[markers_df[, omics_id] %in% markers, "Gene_symbol"]
  
  df2plot <- df[, c(markers, 'sample.description')]
  colnames(df2plot) <- c(marker_genes, 'sample.description')
  df2plot <- reshape2::melt(df2plot, id.vars = 'sample.description')
  colnames(df2plot) <- c('SampleID', 'Gene', 'Expression')
  df2plot$Expression <- as.numeric(df2plot$Expression) 
  return(df2plot)
  
}

## Plot grouped boxplots
boxplot_groups <- function(df2plot, title='Boxplots per group') {
    ggp <- ggplot(df2plot, 
               aes(x = Gene, y = Expression, fill = SampleID)) +
          geom_boxplot(position = position_dodge(width = 0.9))  +
          theme_bw()  +
          ggtitle(title)
    return(ggp)
}



## Heatmap function for neuronal markers
make_heatmap <- function(df, labels, labels_col, markers_df) {
  ggp <- pheatmap(
    df,
    cluster_rows = T,
    cluster_cols = T,
    show_colnames = T,
    show_rownames = T,
    labels_row = labels,
    labels_col = labels_col,
    annnotation_col = markers_df["Stage"],
    scale = 'column',
    fontsize_row = 15,
    fontsize_col = 15
  )
  return(ggp)
}



## Correlation plot for two vector (x and y) equal in size
cor_plot <- function(x, y, x_lab, y_lab, title = "Correlation", method='pearson') {
  
  ## Pearson correlation
  cor_value <- cor(x, y, use="complete.obs", method=method)

  ## ggplot
  ggp <- ggplot(data.frame(x, y), aes(x=x, y=y)) +
    geom_point(alpha =0.6, size=2) +
    geom_smooth(method = "lm", se = T, color = "#2C3E50", fill = "#3498DB") +
    annotate("text",
             x = Inf, y = -Inf,
             label = paste0("r = ", round(cor_value, 3)),
             hjust =1.1, vjust = -0.5,
             size = 4.5,
             fontface = "italic") +
    labs(title = title,
         x = x_lab, y = y_lab) +
    theme_minimal()
  
  return(ggp)

}