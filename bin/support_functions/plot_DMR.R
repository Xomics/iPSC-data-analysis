## Plot DMR
## @params df Is data frame with CpG locations, beta-value, and sample group (Clone)
## Return ggp object
plot_DMR <- function(df, coord_start, coord_end, gene_name) {
  df_plot <- df[df$Genome_Location > coord_start & df$Genome_Location < coord_end, ] 
  n_CpG <- length(unique(df_plot$Genome_Location))
  max_beta <- max(df_plot$Avr_Beta + 0.05)
  ggp <- ggplot(df_plot, aes(x = Genome_Location, y = Avr_Beta, color = Clone, group = Clone)) +
    geom_line(size=1) +
    labs(title = paste0(n_CpG, ' CpGs, DMR mapped to ', gene_name), 
         x= 'GenomeLocation',
         y ='Beta-value') +
    coord_cartesian(ylim =c(0, max_beta)) +
    theme_bw()
  return(ggp)
}