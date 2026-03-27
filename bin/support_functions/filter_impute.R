##################################################
## Title: Impute_median.R
## Description: Functions to impute missing values by median (column-wise) and filter columns on missigness
## Author: Casper de Visser
## Email: casper.devisser@radboudumc.nl
##################################################


#' Impute by median (column-wise)
#' 
#' @param Data Dataframe 
# '@return Imputed data frame by median (column-wise)
impute_median <- function(data) {
  median_values <- sapply(data, median, na.rm=TRUE)
  for (col in names(data)) {
    data[is.na(data[,col]), col] <- median_values[col]
  }
  return(data)
}

#' Remove columns with 20% or more NA values (column-wise)
#' 
#' @param Data Dataframe with missing values
# '@return Filtered data frame 
remove_column_missingness <- function(data, treshold=0.2) {
  missingness <- colMeans(is.na(data))
  column_remove <- names(missingness)[missingness > treshold]
  data <- data[, !(names(data) %in% column_remove)]
  return(data)
  
}