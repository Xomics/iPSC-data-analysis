###############################################################################
# This file contains some general functions used for the metabolomics analysis 
# pipeline. The pipeline is based on result files from MS-DIAL.
#
#
# Author: R.J.E. Derks
# Center for Proteomics & Metabolomics, LUMC
###############################################################################

#' @title Read MS-DIAL files
#'
#' @description Read the result files from MS-DIAL.
#'
#' @param filename The filename of the MS-DIAL file.
#'
#' @return Returns a tibble
#'
#' @importFrom readr read_delim cols col_double col_character col_integer
#'
#' @author Rico Derks
#'
read_msdial <- function(filename) {
  # determine which version is loaded, read only the column names
  column_names <- colnames(readr::read_delim(file = filename,
                                             delim ="\t",
                                             na = c("", "NA", "null"),
                                             n_max = 1,
                                             skip = 4,
                                             show_col_types = FALSE))
  
  if("Simple dot product" %in% column_names) {
    # version > 5.1
    res <- readr::read_delim(file = filename,
                             delim ="\t",
                             na = c("", "NA", "null"),
                             col_types = readr::cols(`Alignment ID` = readr::col_integer(),
                                                     `Average Rt(min)` = readr::col_double(),
                                                     `Average Mz` = readr::col_double(),
                                                     `Metabolite name` = readr::col_character(),
                                                     `Adduct type` = readr::col_character(),
                                                     `Post curation result` = readr::col_character(),
                                                     `Fill %` = readr::col_double(),
                                                     `MS/MS assigned` = readr::col_character(),
                                                     `Reference RT` = readr::col_character(),
                                                     `Reference m/z` = readr::col_character(),
                                                     Formula = readr::col_character(),
                                                     Ontology = readr::col_character(),
                                                     INCHIKEY = readr::col_character(),
                                                     SMILES = readr::col_character(),
                                                     `Annotation tag (VS1.0)` = readr::col_character(),
                                                     `RT matched` = readr::col_character(),
                                                     `m/z matched` = readr::col_character(),
                                                     `MS/MS matched` = readr::col_character(),
                                                     Comment = readr::col_character(),
                                                     `Manually modified for quantification` = readr::col_character(),
                                                     `Manually modified for annotation` = readr::col_character(),
                                                     `Isotope tracking parent ID` = readr::col_character(),
                                                     `Isotope tracking weight number` = readr::col_character(),
                                                     `Total score` = readr::col_double(),
                                                     `RT similarity` = readr::col_double(),
                                                     `Simple dot product` = readr::col_double(),
                                                     `Reverse dot product` = readr::col_double(),
                                                     `Weighted dot product` = readr::col_double(),
                                                     `Matched peaks percentage` = readr::col_double(),
                                                     `S/N average` = readr::col_double(),
                                                     `Spectrum reference file name` = readr::col_character(),
                                                     `MS1 isotopic spectrum` = readr::col_character(),
                                                     `MS/MS spectrum` = readr::col_character()),
                             skip = 4,
                             show_col_types = FALSE)
  } else {
    # version 4.x
    res <- readr::read_delim(file = filename,
                             delim ="\t",
                             na = c("", "NA", "null"),
                             col_types = readr::cols(`Alignment ID` = readr::col_integer(),
                                                     `Average Rt(min)` = readr::col_double(),
                                                     `Average Mz` = readr::col_double(),
                                                     `Metabolite name` = readr::col_character(),
                                                     `Adduct type` = readr::col_character(),
                                                     `Post curation result` = readr::col_character(),
                                                     `Fill %` = readr::col_double(),
                                                     `MS/MS assigned` = readr::col_character(),
                                                     `Reference RT` = readr::col_character(),
                                                     `Reference m/z` = readr::col_character(),
                                                     Formula = readr::col_character(),
                                                     Ontology = readr::col_character(),
                                                     INCHIKEY = readr::col_character(),
                                                     SMILES = readr::col_character(),
                                                     `Annotation tag (VS1.0)` = readr::col_character(),
                                                     `RT matched` = readr::col_character(),
                                                     `m/z matched` = readr::col_character(),
                                                     `MS/MS matched` = readr::col_character(),
                                                     Comment = readr::col_character(),
                                                     `Manually modified for quantification` = readr::col_character(),
                                                     `Manually modified for annotation` = readr::col_character(),
                                                     `Isotope tracking parent ID` = readr::col_character(),
                                                     `Isotope tracking weight number` = readr::col_character(),
                                                     `Total score` = readr::col_double(),
                                                     `RT similarity` = readr::col_double(),
                                                     `Dot product` = readr::col_double(),
                                                     `Reverse dot product` = readr::col_double(),
                                                     `Fragment presence %` = readr::col_double(),
                                                     `S/N average` = readr::col_double(),
                                                     `Spectrum reference file name` = readr::col_character(),
                                                     `MS1 isotopic spectrum` = readr::col_character(),
                                                     `MS/MS spectrum` = readr::col_character()),
                             skip = 4,
                             show_col_types = FALSE)
  }
  
  
  
  return(res)
}


#' @title Clean up tibble from MS-DIAL
#'
#' @description Clean up the columns and column names of the tibble after
#'     reading the MS-DIAL result files.
#'
#' @param met_data tibble, containing all the data (from read_msdial()).
#'
#' @return Returns a tibble
#'
#' @importFrom dplyr rename mutate select
#' @importFrom tidyselect matches
#' @importFrom tidyr unnest
#' @importFrom rlang .data
#'
#' @author Rico Derks
#'
clean_up <- function(met_data) {
  # make a single dataframe
  met_data <- met_data |>
    dplyr::select(polarity, raw_data) |>
    tidyr::unnest(cols = c(polarity, raw_data))
  
  column_names <- colnames(met_data)
  
  if("Simple dot product" %in% column_names) {
    # version > 5.1
    # rename some columns in the data frame for ease of access later on.
    met_data <- met_data |>
      dplyr::rename(AlignmentID = `Alignment ID`,
                    AverageRT = `Average Rt(min)`,
                    AverageMZ = `Average Mz`,
                    ion = `Adduct type`,
                    MetaboliteName = `Metabolite name`,
                    DotProduct = `Simple dot product`,
                    RevDotProduct = `Reverse dot product`,
                    TotalScore = `Total score`,
                    FragPresence = `Matched peaks percentage`,
                    RefFile = `Spectrum reference file name`,
                    MSMSspectrum = `MS/MS spectrum`) |>
      dplyr::mutate(DotProduct = DotProduct * 100,
                    RevDotProduct = RevDotProduct * 100,
                    scale_DotProduct = DotProduct / 10,
                    scale_RevDotProduct = RevDotProduct / 10,
                    my_id = paste(polarity, AlignmentID, sep = "_"),
                    ion = factor(ion),
                    polarity = factor(polarity)) |>
      dplyr::select(my_id, AlignmentID, AverageRT, AverageMZ, ion, MetaboliteName, 
                    DotProduct, scale_DotProduct, RevDotProduct, scale_RevDotProduct,
                    FragPresence, TotalScore, polarity, MSMSspectrum,
                    tidyselect::matches("^([qQ][cC]pool|[sS]ample|[bB]lank)_?.*[0-9]{3}$"))
  } else {
    # version 4.x
    # rename some columns in the data frame for ease of access later on.
    met_data <- met_data |>
      dplyr::rename(AlignmentID = `Alignment ID`,
                    AverageRT = `Average Rt(min)`,
                    AverageMZ = `Average Mz`,
                    ion = `Adduct type`,
                    MetaboliteName = `Metabolite name`,
                    DotProduct = `Dot product`,
                    RevDotProduct = `Reverse dot product`,
                    TotalScore = `Total score`,
                    FragPresence = `Fragment presence %`,
                    RefFile = `Spectrum reference file name`,
                    MSMSspectrum = `MS/MS spectrum`) |>
      dplyr::mutate(scale_DotProduct = DotProduct / 10,
                    scale_RevDotProduct = RevDotProduct / 10,
                    my_id = paste(polarity, AlignmentID, sep = "_"),
                    ion = factor(ion),
                    polarity = factor(polarity)) |>
      dplyr::select(my_id, AlignmentID, AverageRT, AverageMZ, ion, MetaboliteName, 
                    DotProduct, scale_DotProduct, RevDotProduct, scale_RevDotProduct,
                    FragPresence, TotalScore, polarity, MSMSspectrum,
                    matches("^([qQ][cC]pool|[sS]ample|[bB]lank)_?.*[0-9]{3}$"))
  }
  

  return(met_data)
}


#' @title Make the tibble in tidy format
#'
#' @description Filter the tibble to keep only the identified lipids.
#'
#' @param met_data tibble, containing data in wide format (from clean_data()).
#'
#' @details After making the tibble in long format also some additional columns are added.
#'
#' @return Returns a tibble in tidy (long) format
#'
#' @importFrom dplyr mutate n group_by ungroup arrange if_else filter
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches
#' @importFrom rlang .data
#' @importFrom stringr str_extract str_replace
#'
#' @author Rico Derks
#'
tidy_data <- function(met_data) {
  df_long <- met_data |>
    tidyr::pivot_longer(cols = tidyselect::matches("^([sS]ample|[qQ][cC]pool|[bB]lank).*"),
                        names_to = "SampleName",
                        values_to = "PeakArea")  |> 
    dplyr::mutate(
      SampleType = factor(tolower(stringr::str_extract(string = SampleName,
                                                       pattern = "([bB]lank|[qQ][cC]pool|[sS]ample)")))) |> 
    dplyr::relocate(SampleType, .after = SampleName) |> 
    dplyr::filter(!grepl(x = SampleName,
                         pattern = "\\.\\.\\."))
  
  # ### rename duplicate lipids
  # df_long <- df_long |>
  #   # determine what are the duplicates
  #   dplyr::group_by(ShortLipidName, sample_name) |>
  #   dplyr::arrange(AverageRT) |>
  #   dplyr::mutate(count_duplicates = dplyr::n(),
  #                 append_name = paste0("_", 1:dplyr::n())) |>
  #   dplyr::ungroup() |>
  #   # rename them
  #   dplyr::mutate(ShortLipidName = dplyr::if_else(count_duplicates > 1,
  #                                                 paste0(ShortLipidName, append_name),
  #                                                 ShortLipidName)) |>
  #   # sort back
  #   dplyr::arrange(LipidClass, ShortLipidName) |>
  #   dplyr::select(-count_duplicates, -append_name)
  
  return(df_long)
}


#' @title RSD histogram
#' 
#' @description
#' Show a histogram of all RSD values.
#' 
#' @param data tibble, containing all data.
#' @param rsd_limit numeric(1), the rsd limit, between 0 and 1.
#' @param type character(1), show the overall rsd plot or the rsd plot per 
#'     sequence.
#' 
#' @return A ggplot2 object.
#' 
#' @author Rico Derks
#' 
rsd_plot <- function(data = NULL,
                     rsd_limit = 0.3,
                     type = c("overall", "sequence")) {
  # sanity checking
  if(is.null(data)) {
    stop("No data to show!")
  }
  
  if(rsd_limit <= 0 & rsd_limit >= 1) {
    stop("'rsd_limit' should be between 0 and 1!")
  }
  
  type <- match.arg(arg = type,
                    choices = c("overall", "sequence"))
  
  # do stuff
  res <- switch(
    type,
    "overall" = rsd_plot.overall(data = data,
                                 rsd_limit = rsd_limit),
    "sequence" = rsd_plot.rsd(data = data,
                              rsd_limit = rsd_limit)
  )
  
  return(res)
}

#' @title Overall RSD histogram
#' 
#' @description
#' Show a histogram of all RSD values.
#' 
#' @param data tibble, containing all data.
#' @param rsd_limit numeric(1), the rsd limit, between 0 and 1.
#'     sequence.
#' 
#' @return A ggplot2 object.
#' 
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal
#'     theme annotate
#' 
#' @author Rico Derks
#'
rsd_plot.overall <- function(data = data,
                             rsd_limit = 0.3) {
  # sanity checking
  if(is.null(data)) {
    stop("No data to show!")
  }
  
  if(rsd_limit <= 0 & rsd_limit >= 1) {
    stop("'rsd_limit' should be between 0 and 1!")
  }
  
  # make plot
  p <- data |> 
    ggplot2::ggplot(ggplot2::aes(x = rsd,
                                 fill = polarity)) +
    ggplot2::geom_histogram(bins = 50,
                            alpha = 0.5) +
    ggplot2::geom_vline(xintercept = rsd_limit,
                        color = "red",
                        linetype = 2) +
    ggplot2::annotate(geom = "label",
                      x = rsd_limit,
                      y = +Inf,
                      size = 3,
                      label = rsd_limit,
                      vjust = 1) +
    ggplot2::labs(title = "Overall RSD",
                  x = "RSD",
                  y = "count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  return(p)
}


#' @title RSD histogram per sequence
#' 
#' @description
#' Show a histogram of all RSD values.
#' 
#' @param data tibble, containing all data.
#' @param rsd_limit numeric(1), the rsd limit, between 0 and 1.
#'     sequence.
#' 
#' @return A ggplot2 object.
#' 
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal
#'     theme element_text
#' 
#' @author Rico Derks
#'
rsd_plot.rsd <- function(data = data,
                         rsd_limit = 0.3) {
  # sanity checking
  if(is.null(data)) {
    stop("No data to show!")
  }
  
  if(rsd_limit <= 0 & rsd_limit >= 1) {
    stop("'rsd_limit' should be between 0 and 1!")
  }
  
  # create strip labels
  strip_label <- function(string) {
    string <- paste0("Seq. ", string)
    return(string)
  }
  
  # make plot
  p <- data |> 
    ggplot2::ggplot(ggplot2::aes(x = rsd,
                                 fill = polarity)) +
    ggplot2::geom_histogram(bins = 50,
                            alpha = 0.5) +
    ggplot2::geom_vline(xintercept = rsd_limit,
                        color = "red",
                        linetype = 2) +
    ggplot2::labs(title = "RSD per sequence",
                  y = "count") +
    ggplot2::facet_wrap(. ~ Sequence,
                        labeller = ggplot2::labeller(Sequence = strip_label)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text.x = ggplot2::element_text(face = "bold"))
  
  return(p)
}


#' @title Show the trend of several features
#' 
#' @description
#' Show the trend of several features in the QCpool samples. The raw data is used.
#' 
#' @param data data.frame with all data
#' 
#' @return ggplot2 object showing the trend of the features.
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_y_continuous 
#'     facet_wrap guides guide_legend labs theme_minimal theme element_text
#' @importFrom dplyr mutate
#' 
#' @author Rico Derks
#' 
plot_trend <- function(data = NULL) {
  p <- data |> 
    dplyr::mutate(StripName = paste0(MetaboliteName, " (", polarity, ")")) |> 
    ggplot2::ggplot(ggplot2::aes(x = SampleName,
                                 y = PeakArea,
                                 colour = Sequence)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(aes(group = MetaboliteName)) +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                                labels = function(x) format(x, scientific = TRUE)) +
    ggplot2::facet_wrap(~ StripName,
                        ncol = 2,
                        scales = "free") +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Batch")) +
    ggplot2::labs(x = "Acquisition order",
                  y = "Peak area",
                  title = "Trend over all batches",
                  subtitle = "raw data") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1))
  
  return(p)
}


#' @title Show the trend of all features
#' 
#' @description
#' Show the trend of all features in the QCpool samples. All QCpool samples are 
#' compared to the first QCpool sample.
#' 
#' @param data data.frame with all data
#' 
#' @return ggplot2 object showing the trend of the features.
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs theme_minimal 
#'     theme element_text
#' 
#' @author Rico Derks
#' 
plot_all_trends <- function(data = NULL) {
  p <- data |> 
    ggplot(aes(x = SampleName,
               y = log2fc)) +
    geom_line(aes(group = my_id),
              alpha = 0.3) +
    geom_hline(yintercept = c(-0.5, 0.5),
               colour = "red",
               linetype = 2) +
    labs(x = "Sample name",
         y = "Log2(fold change)",
         title = "Trend of all metabolites") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  
  return(p)
}


#' @title Perform QC-RLSC for batch correction of chromatographic signal
#' 
#' @description Perform QC-RLSC for batch correction of chromatographic signal.
#'
#' @param tab table N*K (row * column) with N samples and K variables.
#' @param colv vector N of numbers: 1 for QC samples and 2 for other samples.
#' @param or vector of measuring order (see details).
#' @param span the parameter \alpha which controls the degree of smoothing.
#' @param verbose print which variable has been corrected to monitor the process (default = FALSE).
#'
#' @return corrected table N*K
#' 
#' @details Make sure that everything is sorted in measurement order!!!
#' 
#' @export
#' @importFrom stats loess approx
#'
#' @author E. Nevedomskaya
#' @author Rico Derks
#' @references Dunn et al. Nature Protocols 6, 1060-1083 (2011)
qc_rlsc <- function(tab, colv, or, span = 0.75, verbose = FALSE) {
  # create table of the same size as initial
  tab_corr <- tab
  
  # For each variable (columns) in the initial table
  for (i in 1:ncol(tab)) {
    # fit loess curve to the QCs
    ll <- loess(tab[which(colv == 1), i] ~ or[which(colv == 1)], span = span)
    
    # approximate the curve for all the samples
    aa <- approx(x = or[which(colv == 1)],
                 y = ll$fitted, 
                 xout = or)
    
    # correct the variable according to the curve for all the samples
    tab_corr[, i] <- tab[, i] / aa$y
    
    # print which variable has been corrected in order to monitor the progress
    if(verbose == TRUE) {
      print(i)  
    }
    
  }
  
  return(tab_corr)
}


#' @title Perform Probabilistic Quotient Normalization
#'
#' @description Perform Probabilistic Quotient Normalization
#'
#' @param X matrix to normalize samples * variables (rows * columns)
#' @param n normalization reference: "mean" for using the overall average of variables as reference
#' or "median" (default) for using the overall median of variables as reference
#' @param QC vector of number(s) to specify samples which average to use as reference
#' (e.g. QC samples)
#'
#' @return Normalized table samples * variables (rows * columns)
#'
#' @details First a total area normalization should be done before PQN is applied.
#'
#' @export
#' @importFrom stats median
#'
#' @author E. Nevedomskaya
#' @author Rico Derks
#'
#' @references Dieterle, F., Ross, A., Schlotterbeck, G. & Senn, H. Probabilistic Quotient
#' Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures.
#' Application in H1 NMR Metabonomics. Anal. Chem. 78, 4281-4290 (2006).
pqn <- function(X, n = "median", QC = NULL) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)

  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }

  # do the actual normalisation
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
  }

  return(X.norm)
}
