five_prime_plot <- function(gr, min=NULL, max=NULL){
  #df <- data.frame(width=width(two_mismatches), strand=strand(two_mismatches), five=mcols(two_mismatches)$five)
  df <- data.frame(width=width(gr), strand=strand(gr), five=mcols(gr)$five, stringsAsFactors = FALSE)
  p <- ggplot(data = df) +
    geom_line(aes(x=df$width,
                  group=df$five,
                  color=df$five),
              stat = "count") +
    xlab("Length") +
    ylab("Count") +
    scale_colour_manual(values = c("blue","red","darkgreen","purple")) +
    theme(legend.title=element_blank()) +
    scale_y_continuous(labels=comma) +
    facet_grid(. ~ df$strand, labeller = as_labeller(c("+"="Plus Strand", "-"="Minus Strand")))
  return(p)
}


length_scatter_plot <- function(df, comparison_col, min=NULL, max=NULL){
  melted_df <- data.table::melt(df, id.vars=c(comparison_col, "Gene_strand"))
  p <- ggplot(data = melted_df, aes(x=melted_df[as.character(comparison_col)], y=melted_df$value)) +
    geom_point(inherit.aes = TRUE, size=.5) +
    xlab(paste(comparison_col, "nt reads", sep = "")) +
    ylab("Count") +
    #scale_colour_manual(values = c("blue","red","darkgreen","purple")) +
    #theme(legend.title=element_blank()) +
    scale_x_log10(labels=comma) +
    scale_y_log10(labels=comma) +
    facet_grid(melted_df$variable ~ melted_df$Gene_strand,
               scales = "fixed",
               space = "fixed",
               labeller = as_labeller(c("+"="Plus Strand", "-"="Minus Strand"))) +
    #stat_smooth(geom="text",method="lm",hjust=0,parse=TRUE, inherit.aes = TRUE, fullrange = TRUE) +
    geom_smooth(method="rq", se=FALSE, inherit.aes = TRUE, size=.5)
  ggsave('testplot.png', height = 25, width = 8)
  return(p)
}
