### Function to make pairwise calculation of R-squared, among certain groups
###
###
calc_rsquared_groups <- function(df){
  df_wide <- df %>%
    select(-group) %>%
    pivot_longer(cols = -sample_id, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = sample_id, values_from = value)
  
  r_squared_values <- combn(names(df_wide)[-1], 2, function(cols) {
    group1 <- df$group[df$sample_id == cols[1]]
    group2 <- df$group[df$sample_id == cols[2]]
    
    # Check if groups are the same or different
    if (group1 == group2) {
      cor_value <- cor(df_wide[[cols[1]]], df_wide[[cols[2]]],  use = 'complete.obs')
      r_squared_value <- cor_value^2
      data.frame(pair = paste(cols, collapse = " vs "),
                 r_squared = r_squared_value)
    } else {
      print('No pairwise correlations possible')
    }
  }
  , simplify = FALSE)
  
  bind_rows(r_squared_values)
}


### Function to make dataframe of the R-squares per group
###
###
make_dataframe_rsquared_groups <- function(df, df_metadata, group_column) {
  
  df$group <- df_metadata[[group_column]]
  
  r_squared_by_group <- df %>%
    group_by(group) %>%
    do(calc_rsquared_groups(.))
  
  ## Add column that indicates the group level
  r_squared_by_group$Level <- rep(group_column, nrow(r_squared_by_group))
  
  return(r_squared_by_group)
  
}


### Function to make dataframe of the R-squares per group
###
###
make_dataframe_unrelates_samples <- function(df, df_metadata, group_column) {
  
  # Add column with group info
  df$group <- df_metadata[[group_column]]
  
  # Get the list of numerical columns
  num_cols <- df %>% select(-group) %>% select_if(is.numeric) %>% colnames()
  
  # Initialize an empty list to store correlations
  correlations <- list()
  
  pair_count = 0
  
  # Loop through each pair of rows
  for (i in 1:(nrow(df) - 1)) {
    for (j in (i + 1):nrow(df)) {
      # Check if the groups are different
      if (df$group[i] != df$group[j]) {
        # Calculate correlation for each numerical column
        sample_name1 <- rownames(df)[[i]]
        sample_name2 <- rownames(df)[[j]]
        sample1_values <- as.numeric(df[i, num_cols])
        sample2_values <- as.numeric(df[j, num_cols])
        r_squared <- cor(sample1_values, sample2_values,  use = 'complete.obs')^2
        # Store the correlation values in the list
        correlations[[paste(sample_name1, sample_name2, sep = "-")]] <- r_squared 
        pair_count <- pair_count + 1
        
        if (pair_count > 1000) {
          break
        }
      }
    }
    if (pair_count > 1000) {
      break
    }
  }
  
  # Convert the list to a dataframe for better readability
  cor_df <- do.call(rbind, lapply(names(correlations), function(name) {
    data.frame(pair = name, t(correlations[[name]]))
  }))
  
  # Add columns indicating level and group
  cor_df$Level <-   'unrelated'
  cor_df$group <- 'unrelated'
  colnames(cor_df) <-  c('pair', 'r_squared', 'Level', 'group')
  return(cor_df)
}




### Function to make dataframe of the R-squares per group
###
###
make_dataframe_rsquared_DM1_groups <- function(df, df_metadata, group_column) {
  
  df$group <- df_metadata[[group_column]]
  
  r_squared_by_group <- df %>% 
    group_by(group) %>%
    do(calc_rsquared_groups(.))
  
  ## Add column that indicates the group level
  r_squared_by_group$Level <- rep(group_column, nrow(r_squared_by_group))
  
  return(r_squared_by_group)
  
}


### Function to make boxplots of the R-squares per group
###
###
plot_rsquared_groups_bio <- function(df, isa_assay, assay_name, group_level= "groups") {
  
  ## Add column to color sample identifiers per phenotype
  phenotypes <- lapply(df$group, function(x){
    pheno <-isa_assay[isa_assay$`sample identifier` == x, ncol(isa_assay)]
    pheno[[1]][[1]]
    
  })
  df$phenotype <- as.character(phenotypes)
  
  ## Order sample identifiers based on pheno group
  isa_assay <- isa_assay[order(isa_assay$phenotype),]
  df$group <- factor(df$group, levels = unique(isa_assay$`sample identifier`))
  
  ## Make boxplots
  ggp <- ggplot(df, aes(x=group, y=r_squared, fill = phenotype)) +
    geom_boxplot() +
    ggtitle(paste0(assay_name, '\nPairwise correlations between ', group_level)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    stat_summary(fun.data = function(x){
      return(data.frame(y=mean(x), label = paste0(length(x))))
    }, geom = "text")
  
  return(ggp)
  
}


### Plot pairwise correlations between samples in the same sample group
###
###
plot_rsquared_groups_pheno <- function(df, assay_name, group_level= "groups") {
  
  ## Make boxplots
  ggp <- ggplot(df, aes(x=group, y=r_squared, fill = group)) +
    geom_boxplot() +
    ggtitle(paste0(assay_name, '\nPairwise correlations between ', group_level)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  +
    stat_summary(fun.data = function(x){
      return(data.frame(y=mean(x), label = paste0(length(x))))
    }, geom = "text")
  
  return(ggp)
  
}

## Function to get pairwise correlations on phenotype and bio replicate level
##
## @param name_assay Name of the omics layers in MAE object
## @return Dataframe with pairwise correlation for two levels (phenotype, bio repl.)
get_pairwise_cor_df <- function(name_assay, feature_selection = NULL) {
    
    isa_assay <- isa_mae_list[[name_assay]][[1]]
    df_assay <- isa_mae_list[[name_assay]][[2]]
    
    ## Select features if provided
    if (!is.null(feature_selection)) {
      df_assay <- df_assay[, names(df_assay) %in% feature_selection]
    }
  
    ## Add sample level data from ISA
    df_sample_identifier <- make_dataframe_rsquared_groups(df_assay, isa_assay, 'sample identifier')
    df_phenotype <- make_dataframe_rsquared_groups(df_assay, isa_assay, 'phenotype')
    df_unrelated <- make_dataframe_unrelates_samples(df_assay, isa_assay, 'phenotype')
    
    ## IMPORTANT STEP! REMOVE PAIRS THAT WERE ALSO IN THE OTHER LEVEL
    df_phenotype <- df_phenotype[!df_phenotype$pair %in% df_sample_identifier$pair, ]
    
    ## Combine data frames
    df_input <- rbind(df_unrelated, df_phenotype, df_sample_identifier)
    df_input$Level <- as.factor(df_input$Level)
    df_input$Level <- relevel(df_input$Level, "unrelated")
    return(df_input)
    
}

## Function to get pairwise correlations on phenotype and bio replicate level
##
## @param name_assay Name of the omics layers in MAE object
## @return Dataframe with pairwise correlation for two levels (phenotype, bio repl.)
get_pairwise_cor_df_clones <- function(name_assay, clone_ids, disease) {
  
  isa_assay <- isa_mae_list[[name_assay]][[1]]
  df_assay <- isa_mae_list[[name_assay]][[2]]
  
  ## Select only DM1 patients
  df_assay$group <- isa_assay$phenotype
  isa_assay <- isa_assay[isa_assay$phenotype == disease,]
  df_assay <- df_assay[df_assay$group == disease,]
  df_assay$group <- NULL
  
  df_sample_identifier <- make_dataframe_rsquared_groups(df_assay, isa_assay, 'sample identifier')
  df_phenotype <- make_dataframe_rsquared_groups(df_assay, isa_assay, 'phenotype')
  
  ## IMPORTANT STEP! REMOVE PAIRS THAT WERE ALSO IN THE OTHER LEVEL
  df_phenotype <- df_phenotype[!df_phenotype$pair %in% df_sample_identifier$pair, ]
  
  ## CHANGE DM1 P2 CLONES TO NEW LEVEL
  pairs_matches <- lapply(df_phenotype$pair, function(x) {
    matches <- sapply(clone_ids, grepl, x)
    matches_two <- sum(matches) ==2
    if (matches_two == TRUE){
      x
    }
  })
  pairs_matches[sapply(pairs_matches, is.null)] <- NULL
  df_phenotype <- df_phenotype %>%
    mutate(Level = ifelse(pair %in% pairs_matches, "DM1_P2_clones", Level))
  
  ## Combine data frames
  df_input <- rbind(df_phenotype, df_sample_identifier)
  df_input$Level <- as.factor(df_input$Level)
  df_input$Level <- relevel(df_input$Level, "phenotype")
  return(df_input)
  
}

## Function to calculate pairwise correlation
calc_pairwise_cor <- function(df) {
  
  pairs <- combn(nrow(df), 2, simplify = F)
  map_dfr(pairs, ~{
    row1 <- df[.x[1], ]
    row2 <- df[.x[2], ]
    
    cor_value <- cor(as.numeric(row1), as.numeric(row2))
    r_squared_value <- cor_value^2
    
    data.frame(row1 = .x[1], row2 = .x[2],  r_squared = r_squared_value)
  })
  
}

## Function to calculate pairwise correlations of technical replicates
pairwise_cor_tech <- function(df_input) {
  
  ## Calculate pairwise correlation per group
  result <- df_input %>%
        group_by(group) %>%
        do(calc_pairwise_cor(select(., -group))) %>%
        ungroup()
  ## Add column that indicates the group level
  result$Level <- rep('Technical rep.', nrow(result))
  ## Remove 'pair' column (TODO) to match other data frames
  result$row1 <- NULL
  colnames(result)  <- c("group", "pair", "r_squared", "Level")  
  result$pair <- as.character(result$pair)
  return(result)
}


## Function to make the final plot with all 4 levels
##
## @param name_assay Name of -omics layer in MAE object
## @param df_tech_reps Dataframe with pairwise correlations between technical reps
## @return ggplot
final_combined_plot <- function(df_tech_reps, name_assay, feature_selection = NULL) {
  
  ## Get pairwise correlations for pheno, clone and bio rep levels
  df_input <- get_pairwise_cor_df(name_assay, feature_selection=feature_selection)
  
  ## Combine data frame with technical reps
  df_input <- rbind(df_input, df_tech_reps)
  df_input$Level <- as.factor(df_input$Level)
  df_input$Level <- relevel(df_input$Level, "unrelated")
  ## Change label names
  levels(df_input$Level) <- c("Unrelated samples",  "Same disease group", "Biological repl.", "Tech repl.")
  colors <- c("Unrelated samples" =  "indianred2", "Same disease group" = "forestgreen", "Biological repl." = "skyblue1", "Tech repl." = "purple")
  
  ylim_min <- floor(min(df_input$r_squared * 100)) / 100 - 0.02
    
  ggp <- ggplot(df_input, aes(x=Level, y=r_squared, fill = Level)) +
      geom_boxplot() +
      ggtitle(paste0(name_assay)) +
      theme_bw() +
      ylim(ylim_min, 1) +
      scale_fill_manual(values = colors) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      stat_compare_means() +
      stat_summary(fun.data = function(x){
              return(data.frame(y=mean(x), label = paste0(length(x))))
          }, geom = "text")
  
  return(ggp)
}  