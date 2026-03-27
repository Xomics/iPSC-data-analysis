library(ggplot2)
library(reshape2)

## Boxplot to check the distribtuion of omics levels per sample
##
## @param df Dataframe with Sample_ID column and omics features
boxplot_samples <- function(df, title) {
  df_long <- melt(df, id.vars= "Sample_ID")
  
  ggp <- ggplot(df_long, aes(x = Sample_ID, y = value)) +
    geom_boxplot() +
    labs(title = title, x = 'Sample', y = "Value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust =1))
  return(ggp) 
}