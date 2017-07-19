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

sequence_logo_comparison <- function(gr, method="bits", flanks=0){
  colors_scheme = make_col_scheme(chars=c('A', 'C', 'G', 'T'),
                                  groups=c('A', 'C', 'G', 'T'),
                                  cols=c('blue', 'red', 'green', 'purple'))
  p <- ggplot() +
    geom_logo(gr, col_scheme=colors_scheme, method = method) +
    #theme_logo() +
    theme(panel.grid.minor.x = element_blank()) +
    #coord_cartesian(ylim = c(0, 2)) +
    #coord_cartesian(ylim = c(0, 2), xlim = c(0-flanks, 22+flanks)) +
    #scale_x_discrete(breaks = seq(from = 0,to = 22, by = 5)) +
    #scale_y_discrete(breaks = seq(0, 2, by = .5)) +
    facet_grid(.~seq_group) +
    scale_y_continuous(limits = c(0, 2)) +
    scale_x_continuous(limits = c(0, 24), breaks = pretty_base1(1,2), labels = seq(1,22, by=5))
  return(p)



  #p <- ggplot() + geom_logo(sequences, col_scheme=colors_scheme, method = "bits") + theme(panel.grid.minor.x = element_blank()) + facet_grid(.~seq_group) + scale_y_continuous(limits = c(0, 2)) + scale_x_continuous(limits = c(0, 24), breaks = pretty_base1(1,22))
  #q <- ggplot() + geom_logo(sequences2, col_scheme=colors_scheme, method = "bits") + theme(panel.grid.minor.x = element_blank()) + facet_grid(.~seq_group) + scale_y_continuous(limits = c(0, 2)) + scale_x_continuous(limits = c(1, 42), breaks = pretty_base1(1,45)+1, labels = pretty_base1(-9,35))




  }


