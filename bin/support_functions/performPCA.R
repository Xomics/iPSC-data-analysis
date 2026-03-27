##################################################
## Title: perfomrPCA.R
## Description: Function to perform PCA and generate a prcomp object
## Author: Purva Kulkarni
## Email: purva.kulkarni@radboudumc.nl
##################################################

library(ggplot2)
library(ggbiplot)


performPCA <- function(data_table, start_col_name, end_col_name)
{
  #data_table <- read.xlsx(filePath, sheet = 1)
  
  # Prepare data frame based on columns of interest
  start = as.numeric(which(colnames(data_table) == start_col_name))
  end = as.numeric(which(colnames(data_table) == end_col_name))
  
  data_table_sub <- data_table[,start:end]
  
  #transpose matrix
  data_table_sub <- t(data_table_sub)
  # Convert NA to 0
  data_table_sub[is.na(data_table_sub)] <- 0
  
  # Add a new column to indicate sample categories
  # The categories have been extracted from the rownames, as the first word before row names
  # NANS_p3c1_Proteomics_2_S2_A11_1 is NANS
  sampleType <- sub("^(.*?)_.*","\\1", rownames(data_table_sub))
  colnames(data_table_sub)[ncol(data_table_sub)] <- "sampleType"
  
  data.pca <- prcomp(data_table_sub, center = TRUE, scale.= TRUE)
#   
# 
# # Generate scree plot  
   print(ggscreeplot(data.pca))
#   
# # Generate scores plot
  g <- ggbiplot(
    data.pca,
    groups = sampleType,
    var.axes = FALSE
  )
  g
#   
#   # g <-
#   #   ggbiplot(
#   #     data.pca,
#   #     labels = colnames(data_table_sub),
#   #     obs.scale = 1, 
#   #     var.scale = 1,
#   #     ellipse = TRUE, 
#   #     circle = TRUE,
#   #     varname.color = "red",
#   #     varname.size = 4,
#   #     labels.size = 5,
#   #     var.axes = FALSE
#   #   )
  
#  return(data.pca)
}