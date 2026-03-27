extract_chains <- function(x) {
  str_extract_all(x, "\\d{2}:\\d")[[1]]
}

parse_chain <- function(chain) {
  data.frame(
    Carbons = as.numeric(str_extract(chain, "^\\d{2}")),
    DBs = as.numeric(str_extract(chain, "(?<=:)\\d$"))
  )
}

sum_chains_df <- function(chains) {
  vals <- bind_rows(lapply(chains, parse_chain))
  data.frame(
    Total_Carbons = sum(vals$Carbons),
    Total_DBs = sum(vals$DBs)
  )
}




# Single chain (CE, PA, FA, LPC, LPE)
parse_single_chain_df <- function(df) {
  df %>%
    mutate(
      Chains = map(Feature, extract_chains),
      parsed = map(Chains, ~ parse_chain(.x[1])),
      Total_Carbons = map_dbl(parsed, "Carbons"),
      Total_DBs = map_dbl(parsed, "DBs"),
      Num_Chains = 1,
      Hydroxyls = NA_real_
    ) %>%
    select(-Chains, -parsed)
}



# Two-chain GPLs (PC, PE, PS, PI, PG)
pparse_two_chain_df <- function(df) {
  df %>%
    mutate(
      Chains = map(Feature, extract_chains),
      parsed = map(Chains, sum_chains_df),
      Total_Carbons = map_dbl(parsed, "Total_Carbons"),
      Total_DBs = map_dbl(parsed, "Total_DBs"),
      Num_Chains = 2,
      Hydroxyls = NA_real_
    ) %>%
    select(-Chains, -parsed)
}




# Sphingolipids (Cer, SM, HexCER, LacCER)
parse_sphingo_df <- function(df) {
  df %>%
    mutate(
      BaseType = str_extract(Feature, "(?<=\\.)[dt](?=\\d{2}:\\d)"),
      Chains = map(Feature, extract_chains),
      parsed = map(Chains, sum_chains_df),
      Total_Carbons = map_dbl(parsed, "Total_Carbons"),
      Total_DBs = map_dbl(parsed, "Total_DBs"),
      Num_Chains = map_int(Chains, length),
      Hydroxyls = case_when(
        BaseType == "d" ~ 2,
        BaseType == "t" ~ 3,
        TRUE ~ NA_real_
      )
    ) %>%
    select(-Chains, -parsed, -BaseType)
}




# TG + DG (total value preferred)
parse_glycerolipid_df <- function(df) {
  df %>%
    mutate(

      Chains = map(Feature, extract_chains),
      parsed = map(Chains, sum_chains_df),
      Total_Carbons = map_dbl(parsed, "Total_Carbons"),
      Total_DBs = map_dbl(parsed, "Total_DBs"),
      Num_Chains = map_int(Chains, length),
      Hydroxyls = NA_real_
    ) %>%
    select(-Chains, -parsed)
}



correlation_plot <- function(df, target_column, filter_column= "Class", filter_value, lipid_property = "Chain", method='spearman') {
  # Subset the dataframe
  subset_df <- df %>%
    filter(.data[[filter_column]] == filter_value)

  # Correlation
  cor_value <- cor(subset_df[[target_column]], y = subset_df[[lipid_property]], method = method, use='complete.obs')

  # Create correlation plot
  ggplot(subset_df, aes(x = .data[[target_column]], y = .data[[lipid_property]])) +
    geom_point(color = "steelblue") +  
    labs(
      title = paste("Correlation between", lipid_property, "and", target_column),
      subtitle = paste("Filtered by", filter_column, "=", filter_value, "|", method, "correlation =", round(cor_value, 2)),
      x = target_column,
      y = lipid_property
    ) +
    theme_minimal()
}

correlation_plot_all_groups <- function(df, target_column, filter_column, fixed_column = "fixed_column_name") {
  
  # Spearman correlation per group
  cor_df <- df %>%
    group_by(.data[[filter_column]]) %>%
    summarise(
      rho = suppressWarnings(
        cor(.data[[target_column]], .data[[fixed_column]], method = "spearman", use = "pairwise.complete.obs")
      )
    ) %>%
    mutate(
      label = glue("{.data[[filter_column]]} (ρ={round(rho, 2)})")
    )
  
  # Merge labels back to main df
  df_plot <- df %>%
    left_join(cor_df, by = filter_column)
  
  ggplot(df_plot, aes(
    x = .data[[target_column]],
    y = .data[[fixed_column]],
    color = label
  )) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    labs(
      title = paste("Relationship between", fixed_column, "and", target_column),
      subtitle = paste("Trend lines per", filter_column),
      x = target_column,
      y = fixed_column,
      color = filter_column
    ) +
    theme_minimal()
}


