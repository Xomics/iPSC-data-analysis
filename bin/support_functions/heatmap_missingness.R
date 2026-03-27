##################################################
## Title: heatmap_missingness.R
## Description: Function to plot the a heatmap that shows missing value in a data frame with sample (rows) and features (columns).
## Author: Casper de Visser
## Email: casper.devisser@radboudumc.nl
##################################################

library(ggplot2)

#' Heatmap of missingness data
#' 
#' @param Data Dataframe of which missingness will be represented in heatmap
#' @return Heatmap of the missingness
heatmap_missing <- function(data, title) {
  
  heatmap <- visdat::vis_miss(data, show_perc_col = FALSE, warn_large_data = FALSE) + 
                ggtitle ('Missing values', title) + xlab(paste0('Features: ', ncol(data))) +
                theme(plot.title = element_text(size = 16, hjust = 0.5)) +
                if (ncol(data) > 50) {
                  theme(axis.text.x = element_blank()
                   )} 
  return(heatmap)
}	