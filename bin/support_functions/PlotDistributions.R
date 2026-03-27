##################################################
## Title: PlotDistributions.R
## Description: Function to plot the histograms ("hist"), density plots ("dens") or violinboxplots ("vibo") of all numeric columns.
## Author: Cenna Doornbos
## Email: cenna.doornbos@radboudumc.nl
##################################################

library(ggplot2)
library(tidyverse)

# Create data distributions for all numeric variables in a df in one plot
# use 'type'to specify the type of distribution plot: histograms ("hist"), density plots ("dens") or violinboxplots ("vibo")
PlotDistributions <- function(df, type = "dens"){
  
  # prepare data
  df_long <- df %>%
    pivot_longer(cols = where(is.numeric)) %>%
    drop_na()
  
  # plot data
  if(type == "hist"){
    p <- ggplot(df_long, aes(value)) +
      geom_histogram(bins = 20, color = "grey", fill = "#0073B3") +
      labs(title = "Histograms of Numeric Variables")
  }else if(type == "dens"){
    p <- ggplot(df_long, aes(value)) +
      geom_density(color = "grey", alpha = 0.6, fill = "#0073B3") +
      labs(title = "Density plots of Numeric Variables")
  }else if(type == "vibo"){
    p <- ggplot(df_long, aes(name, value)) +
      geom_violin(fill = "#0073B3", color = "grey", alpha = 0.6) +
      geom_boxplot(alpha = 0.2, width = 0.2) +
      labs(title = "Violinboxplot of Numeric Variables")
  }
  p <- p +
    facet_wrap(~name, scales = "free") +
    labs(x = "Value", y = "Frequency") +
    theme_minimal() +
    theme(axis.text = element_text(size = 6), strip.text = element_text(size = 7))
  return(p)
}
