add_nested_cols <- function(df) {
  library(dplyr)
  sample.renames <- lapply(df$sample.name, function(x){
    if (grepl('repeat excised', x)) {
      x <- x
    }
    else if (grepl(' clone', x)) {
      x <- strsplit(x, 'clone')[[1]][[1]]
    }
    else{
      x <- x
    }
    x
  })
  
  ## Add column as factor
  df$sample.renames <- as.factor(as.character(sample.renames))
  
  ## Column indicating biological replicate
  bio_reps <- lapply(df$sample.identifier, function(x) {
    char <- nchar(x)
    rep <- substr(x, char, char)
    rep
  })
  df$bio_reps <- as.factor(as.character(bio_reps))
  
  ## Add column for clone id (only for DM1 patient 2 and CHD2 patient 1)
  clone_ids <- lapply(df$sample.description, function(x) {
    if (grepl('clone ', x)) {
      x <- strsplit(x, "clone ")[[1]][2]
      x <- substr(x, 1, 1)
    }
    else{
      x <- '1'
    }
    x
  })
  df$clone_ids <- as.factor(as.character(clone_ids))
  return(df)
}  


code_phenotype_column <- function(df) {
  ## Add phenotype as factor
  df$phenotype <- as.factor(df$phenotype)
  
  ## Add Healty as first level
  df$phenotype[df$phenotype == 'Diseased: CHD2 CRISPR'] <- 'Diseased: CHD2'
  df$phenotype <- droplevels(df$phenotype)
  df$phenotype <- relevel(as.factor(df$phenotype), ref = 'Healthy')
  
  ## Add separate column for DM1 comparison
  df$phenotype_DM1 <- df$phenotype
  df$phenotype_DM1[df$phenotype_DM1 == 'Rescue: repeat excised CRISPR DM1'] <- 'Healthy'
  
  ## Add Healty as first level
  df$Disease <- relevel(as.factor(df$Disease), ref = 'Healthy')
  
  return(df)
}


run_MixedModel_df <- function(feature, df, predictor, nest_effect = TRUE) {
  #print(feature)
  # run model for one feature and one outcome
  df$y <- df[[feature]]
  # fit mixed model
  if (nest_effect) {
    #str_formula <- paste0("y ~ ", predictor, " + (1 | sample.renames/bio_reps)")
     str_formula <- paste0("y ~ ", predictor, " +  (1 | sample.renames/clone_ids)")
  }
  else {
    str_formula <- paste0("y ~ ", predictor, " + (1 | sample.renames)")
  }
  my_formula <- as.formula(str_formula)
  lm_fit <- try(lmerTest::lmer(formula = my_formula, data=df))
  lm_fit_df <- as.data.frame(summary(lm_fit)$coefficients)
  # Determine random effect on intercepts
  random_effect <- var(lme4::ranef(lm_fit)$sample.renames)
  # Summarize result in data frame
  rlm_result <- lm_fit_df[rownames(lm_fit_df) != '(Intercept)',]
  colnames(rlm_result) <- c('coefficient', 'std.error', 'df', 't.value', 'p.value')
  rlm_result$feature = feature
  rlm_result$outcome = rownames(rlm_result)
  rlm_result$rand_effect <- as.numeric(random_effect)
  if (nest_effect) {
    #nested_effect <- var(lme4::ranef(lm_fit)$`bio_reps:sample.renames`)
    nested_effect <- var(lme4::ranef(lm_fit)$`clone_ids:sample.renames`)
    rlm_result$nest_effect <- as.numeric(nested_effect)
  }
  return(rlm_result)
}


## Function to run the random model on all disease groups (separate or combined)
run_rand_model <- function(assay_name, pred, DM1_ONLY = F) {
  se <- experiments(mae)[assay_name]
  df <- get_df_input_lm(se)
  features <- colnames(df)
  features <- features[features != "ENSG00000129824"] #TODO
  # add biological variables to data frame
  df <- cbind(df, pheno_df[rownames(df), ])
  df <- add_nested_cols(df)
  df <- code_phenotype_column(df)
  # Select only DM1 and control samples if DM1_ONLY = TRUE
  if (DM1_ONLY) {
    df <-  df[df$phenotype_DM1 %in% c('Diseased: DM1', 'Healthy'),]
  }
  # fit linear mixed model
  model_result <- do.call(rbind, lapply(
    features,
    run_MixedModel_df,
    df=df,
    predictor=pred))
  #adjust p-values
  model_result$p.adj <- p.adjust(model_result$p.value, method = "BH")
  model_result$rank <- rank(model_result$p.value)
  ## Remove the CRISPR comparisons
  model_result <- model_result[!grepl('excised', model_result$outcome), ]
  if (!("nest_effect" %in% names(model_result))) {
    model_result$nest_effect <- 0
    model_result <- model_result[,c(1:8, 11,9,10)]
  }
  model_result$Assay <- assay_name
  model_result$Feature_plot <- model_result$feature #If not other feature ID available
  return(model_result)
}


## Function to split data into chunks
split_in_chunks <- function(df, n_chunks) {
  
  chunk_size <- ceiling(ncol(df) / n_chunks)
  split_df <- lapply(seq(1, ncol(df), by = chunk_size), function(i) df[, i:min(i + chunk_size -1, ncol(df))])
  return(split_df)
}

## Add phenotypic columns to each chunk (Add to source file?)
add_pheno_cols <- function(df) {
  
  # add biological variables to data frame
  df <- cbind(df, pheno_df[rownames(df), ])
  df <- add_nested_cols(df)
  df <- code_phenotype_column(df)
  
  return(df)
}


get_df_input_lm <- function(se) {
  # create data frame with omics feature values
  df <- as.data.frame(t(assays(se)[[1]]))
  # remove non-variable features
  df <- df[, apply(df, 2, function(column) {
    var(column, na.rm = T) != 0 })]
  # remove extremely sparse features
  xpercent <- nrow(df)*0.5
  df <- df[, apply(df, 2, function(column) {
    sum(column == 0, na.rm = T) < xpercent })]
  return(df)
  
}


# Extract and prepare feature data
prepare_feature_data <- function(assay_name, mae, exclude_genes = "ENSG00000129824") {
  se <- experiments(mae)[assay_name]
  df <- get_df_input_lm(se)
  features <- colnames(df)
  features <- features[features != exclude_genes]
  
  list(df = df, features = features)
}

# Regress out confounding factors
regress_out_factors <- function(df, features, factor_columns) {
  df_corrected <- as.data.frame(
    lapply(features, function(col) {
      formula_str <- paste(col, "~", paste(factor_columns, collapse = " + "))
      resid(lm(as.formula(formula_str), data = df))
    })
  )
  colnames(df_corrected) <- features
  df_corrected
}

# Add confounding factors to dataframe
add_confounding_factors <- function(df, mofa_scores_df, factor_names = "Factor2") {
  mofa_scores_match <- mofa_scores_df[match(rownames(df), rownames(mofa_scores_df)), , drop = FALSE]
  
  for (factor in factor_names) {
    df[[factor]] <- mofa_scores_match[[factor]]
  }
  
  df
}






# Prepare final dataframe for modeling
prepare_modeling_data <- function(df_corrected, df_original, pheno_df, DM1_ONLY = FALSE) {
  # Add biological variables
  df <- cbind(df_corrected, pheno_df[rownames(df_original), ])
  df <- add_nested_cols(df)
  df <- code_phenotype_column(df)
  
  # Filter for DM1 samples if requested
  if (DM1_ONLY) {
    df <- df[df$phenotype_DM1 %in% c('Diseased: DM1', 'Healthy'), ]
  }
  
  df
}


# Run model for corrected data
run_models_and_save <- function(df, features, predictor, assay_name) {
  # Fit models
  model_result <- do.call(rbind, lapply(
    features,
    run_MixedModel_df,
    df = df,
    predictor = predictor
  ))
  
  # Adjust p-values
  model_result$p.adj <- p.adjust(model_result$p.value, method = "BH")
  model_result$rank <- rank(model_result$p.value)

  model_result$X <- NULL
  ## Remove the CRISPR comparisons
  model_result <- model_result[!grepl('excised', model_result$outcome), ]
  if (!("nest_effect" %in% names(model_result))) {
    model_result$nest_effect <- 0
    model_result <- model_result[,c(1:8, 11,9,10)]
  }
  model_result$p.adj <- p.adjust(model_result$p.value, method = "BH")
  model_result$Assay <- assay_name
  model_result$Feature_plot <- model_result$feature #If not other feature ID available
  model_result
}


# Main function
run_rand_model_CORRECTED <- function(assay_name, pred, DM1_ONLY = FALSE, 
                                     confound_factors = "Factor2") {
  # 1. Prepare data
  data_prep <- prepare_feature_data(assay_name, mae)
  df <- data_prep$df
  features <- data_prep$features
  
  # 2. Add confounding factors
  df <- add_confounding_factors(df, mofa_scores_df, confound_factors)
  
  # 3. Regress out confounding factors
  df_corrected <- regress_out_factors(df, features, confound_factors)
  
  # 4. Prepare final modeling dataframe
  df <- prepare_modeling_data(df_corrected, df, pheno_df, DM1_ONLY)
  
  # 5. Run models and save
  model_result <- run_models_and_save(df, features, pred, assay_name)
  
  return(model_result)
}
