#' Convert omics data to long format
wide_to_long_keep_samples <- function(df) {
  df$sample <- rownames(df)
  gather(df, key = "feature", value = "value", -sample, factor_key = TRUE)
}

#' Prepare MultiAssayExperiment for MOFA
#'
#' @param mae A MultiAssayExperiment.
#' @param assay_selection An optional character vector with assay names to be included in analysis.
#' @returns A long data frame.
mae_to_MOFAinput <- function(mae) {

  # get omics data from SummarizedExperiments
  assay_list <- lapply(assay_names, function(assay_name) {
    se <- mae[[assay_name]]
    df_wide <- data.frame(t(assay(se)))
    # standardize data
    df_wide <- as.data.frame(t(scale(t(df_wide), center = TRUE, scale = TRUE)))
    # convert wide to long data frame
    df_long <- wide_to_long_keep_samples(df_wide)
    shorter_name <- ifelse(
      grepl("_metabolite_profiling", assay_name),
      strsplit(assay_name, "_")[[1]][2],
      shoter_names[[assay_name]])
    df_long$view <- shorter_name
    df_long
  }) 
  merged_assay_data <- do.call(rbind, assay_list)

}


## Function to plot scores of two data MOFA factor scores and color by group
Factor_scores_plot <- function(data, Col1, Col2, color_by = 'phenotype', label_col = 'observation.unit.identifier.crispr') {
  ggp <- ggplot(data, aes_string(x= Col1, y= Col2, color = color_by)) + 
    geom_point(size =2.5, alpha =0.85) +  
    theme_classic(base_size=12) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=16), 
          legend.title = element_text(size =13),
          panel.background = element_rect(fill = "white")) +
    geom_text_repel(aes_string(label = label_col)) +
    
  return(ggp)
}


## Function to please feature weights of lipidomics data, color by class
plot_feature_weights <- function(model, view, factor, feature_metadata,
                                 top_n = NULL, color_col = "group", dot_size=2) {
  # 1. extract weights
  w <- get_weights(MOFAobject.trained, views = view, factors = factor)[[1]]
  w <- data.frame(
    feature = rownames(w),
    weight  = as.numeric(w),
    stringsAsFactors = FALSE
  )
  
  # 2. merge with metadata
  lipidomics_feature_metadata$feature <- rownames(lipidomics_feature_metadata)
  df <- w %>%
      left_join(lipidomics_feature_metadata)  
  
  # order by weight
  df <- df %>%
    arrange(desc(weight)) %>%
    mutate(rank = row_number())
	
  if (!is.null(top_n)) {
	df_top <- df%>%
		slice_min(weight, n = top_n) %>%
		bind_rows(
			df %>% slice_max(weight, n = top_n)
		) %>%
		arrange(desc(weight)) %>%
		mutate(rank = row_number())
  } else {
	df_top <- df
  }
  p <- ggplot(df_top, aes(x =rank, y=weight, color = Class)) +
    	geom_point(alpha =0.7, size = dot_size) +
	scale_color_viridis_d(option = "turbo") +
    	theme_minimal() +
      	labs(
		title = paste0("Weights on factor ", factor),
		x = "Feature",
		y  = "Weight",
		color = color_col
      	) 
  
  return(p)
}
