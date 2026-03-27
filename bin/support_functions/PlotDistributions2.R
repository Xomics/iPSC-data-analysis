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
PlotDistributions <- function(df, type = "dens", main_title = "Distribution of Numeric Variables"){
  
  # prepare data
  df_long <- df %>%
    pivot_longer(cols = where(is.numeric)) %>%
    drop_na()
  
  # plot data
  if(type == "hist"){
    p <- ggplot(df_long, aes(value)) +
      geom_histogram(bins = 20, color = "grey", fill = "#0073B3") 
  }else if(type == "dens"){
    p <- ggplot(df_long, aes(value)) +
      geom_density(color = "grey", alpha = 0.6, fill = "#0073B3") 
  }else if(type == "vibo"){
    p <- ggplot(df_long, aes(name, value)) +
      geom_violin(fill = "#0073B3", color = "grey", alpha = 0.6) +
      geom_boxplot(alpha = 0.2, width = 0.2) 
  }
  return(p + labs(title = main_title))
}
