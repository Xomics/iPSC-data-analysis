## Create named lists of ranking metric
get_ranks_vector <- function(df, factor, col=NULL) {
  df <- df[order(df[, factor], decreasing = T),]
  ranks <- df[, factor]
  if (!is.null(col)) {
    names(ranks) <- df[, col]  
  } else {
    names(ranks) <- rownames(df)
  }
  return(ranks)
}


## Pathway enrichment analysis
multiGSEA_new <- function(pathways, ranks, eps = 1e-100, scoreType = "std") {
  # function modifed from multiGSEA::multiGSEA 
  # Canzler, Sebastian and Hackermüller, Jörg: multiGSEA: A GSEA-based pathway 
  # enrichment analysis for multi-omics data. BMC Bioinformatics 21, 561 (2020).
  # https://doi.org/10.1186/s12859-020-03910-x
  # allows to change parameters of fgsea::fgseaMultilevel because
  # eps=0 (boundary to calculate p-value) can result in very long running times
  # note that scoreType = "std" test both directions of ranked list
  require(fgsea)
  
  # Go through all omics layers
  es <- lapply(names(pathways), function(omics) {
    fgseaMultilevel(pathways = pathways[[omics]], stats = ranks, 
                    eps = eps, scoreType = scoreType)
  })
  
  names(es) <- names(pathways)
  
  return(es)
}


## Plot significant pathways
plot_patwhays <- function(df.signif, title) {
  if ("Database" %in% colnames(df.signif)) {
    filename.out <- paste0(title, "_", df.signif$Database[[1]])
    title <- paste0(title, "  (", df.signif$Database[[1]], ")")
    
  }
  # Plot NES
  ggp <- ggplot(df.signif) +
    geom_segment(aes(y = pathway, x = GeneRatio, xend = 0, yend = pathway)) +
    geom_point(aes(y = pathway, x = GeneRatio, fill=padj, size = size), shape = 21) +
    scale_fill_viridis_c() +
    facet_grid(Enrichment ~ ., scales = "free", space = "free") +
    ggtitle(title)
  
  filename.out <- paste0('Top_pathways_', filename.out, '.png')
  
  ggsave(ggp, filename = filename.out)
  
  return(ggp)
}


# Select only signfiic
prepare_sign_hits_df <- function(df) {
  # extract significantly enriched pathways
  df.signif <- df[df$padj < 0.05, ]
  return(df.signif)
}


# subset results for plotting if too many significant pathways
prepare_plot_df <- function(df.signif, db = FALSE) {
  
  # Remove pathways that share the same leading edge genes
  df.signif <- df.signif[!duplicated((as.character(df.signif$leadingEdge))),]
  
  # Add column indicating negative or positive enrichment
  df.signif$Enrichment <- NA
  enrichment_col <- lapply(df.signif$NES, function(x){
    if (x > 0) {
      value <- 'Positively associated'
    }
    else {
      value <- 'Negatively associated'
    }
    value
  })
  df.signif$Enrichment <- as.character(enrichment_col)
  
  ## Add column indicating database used
  df.signif$Database <- sapply(strsplit(as.character(df.signif$pathway), "\\)"), function(x) sub(".*\\(", "", x[1]))
  
  if (db != FALSE) {
    df.signif <- df.signif[df.signif$Database == db,]
  }
  
  if (nrow(df.signif) > max_num_pathways) {
    # Make separate df for negative and positive enrichments
    df.signif_pos <- df.signif[df.signif$NES > 0, ]
    df.signif_neg <- df.signif[df.signif$NES < 0, ]
    
    # Order based on NES
    df.signif_pos <- df.signif_pos[order(NES, decreasing = T),][1:5,]
    df.signif_neg <- df.signif_neg[order(NES, decreasing = F),][1:5,]
    
    # Bind data frames
    df.signif <- rbind(df.signif_pos, df.signif_neg)
    df.signif <- na.omit(df.signif)
  }
  
  # Add Gene ratio column
  gene_ratios <- lapply(seq(1, length(df.signif$leadingEdge)), function(x) { 
    leadingedge <- df.signif[x,8]
    leadingedge_n <- length(leadingedge[[1]][[1]])
    size <- df.signif[x,7]
    gene_ratio <- as.numeric(leadingedge_n) / as.numeric(size) 
    round(gene_ratio, digits = 3)
  })
  df.signif$GeneRatio <- as.numeric(as.character(gene_ratios))
  
  
  #order pathways by normalized enrichment score
  df.signif$pathway <- factor(
    stringr::str_wrap(df.signif$pathway, width = 50),
    levels = stringr::str_wrap(df.signif$pathway, width = 50)[order(df.signif$NES)],
    ordered = T)
  return(df.signif)
}


## LeadingEdge genes, translate ENSEMBL to Gene Symbols (mRNA-seq)
translateLeadingEdgeGenes <- function(df, feature_metadata, col_nr) {
  leadingEdge_col <- df$leadingEdge
  
  leadingEdge_symb <- lapply(leadingEdge_col, function(x) {
    
    ## Remove "\n" that is present in some ENSEMBL IDs
    x <- (gsub("\n", "", x))
    
    leading_genes <- x
    symbols <- lapply(leading_genes, function(x) {
      sym <- feature_metadata[feature_metadata$feature==x, col_nr] # Use result dataframe with gene symbols
      sym <- (gsub("\n", "", sym))
      sym
    }) 
    list(as.character(symbols))
    
  })
  df$leadingEdge_Symbols <- as.character(leadingEdge_symb)
  return(df)
}