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


length_scatter_plot <- function(gr, regions, x_length, overlap = "sense", min=NULL, max=NULL){

  df <- count_overlaps_by_width(gr = gr, regions = regions, overlap = overlap)
  melted_df <- data.table::melt(df, id.vars=c(x_length, "Gene_strand"))
  p <- ggplot(data = melted_df) +
    geom_point(aes(x=melted_df[as.character(x_length)], y=melted_df$value)) +
    xlab(paste(x_length, "nt reads", sep = "")) +
    ylab("Count") +
    #scale_colour_manual(values = c("blue","red","darkgreen","purple")) +
    #theme(legend.title=element_blank()) +
    #scale_x_continuous(labels=comma) +
    #scale_y_continuous(labels=comma) +
    facet_grid(melted_df$variable ~ melted_df$Gene_strand,
               space = "fixed",
               widths = 1:4, heights = 4:1,
               labeller = as_labeller(c("+"="Plus Strand", "-"="Minus Strand")))
    return(p)
}
