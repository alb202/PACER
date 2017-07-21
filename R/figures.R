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
    facet_wrap(.~seq_group) #+
    scale_y_continuous(limits = c(0, 2)) +
    scale_x_continuous(limits = c(0, 24), breaks = pretty_base1(1,2), labels = seq(1,22, by=5))
  return(p)



  #p <- ggplot() + geom_logo(sequences, col_scheme=colors_scheme, method = "bits") + theme(panel.grid.minor.x = element_blank()) + facet_grid(.~seq_group) + scale_y_continuous(limits = c(0, 2)) + scale_x_continuous(limits = c(0, 24), breaks = pretty_base1(1,22))
  #q <- ggplot() + geom_logo(sequences2, col_scheme=colors_scheme, method = "bits") + theme(panel.grid.minor.x = element_blank()) + facet_grid(.~seq_group) + scale_y_continuous(limits = c(0, 2)) + scale_x_continuous(limits = c(0, 42), breaks = pretty_base1(1,45)+1, labels = pretty_base1(-9,35))

  sequences <- list("plus_seq"=as.character(mcols(gr[width(gr)==22 & strand(gr)=="+"])$seq),
                    "minus_seq"=as.character(mcols(gr[width(gr)==22 & strand(gr)=="-"])$seq))
  sequences2 <- list("plus_int"=as.character(mcols(gr[width(gr)==22 & strand(gr)=="+"])$with_flanks),
                     "minus_int"=as.character(mcols(gr[width(gr)==22 & strand(gr)=="-"])$with_flanks))
  sequences3 <- c(sequences, sequences2)

seq_logos <- function(gr, length){
  colors_scheme = make_col_scheme(chars=c('A', 'C', 'G', 'T'),
                                  groups=c('A', 'C', 'G', 'T'),
                                  cols=c('blue', 'red', 'green', 'purple'))

  interval_plus_data <- as.character(mcols(gr[width(gr)==22 & strand(gr)=="+"])$seq)
  interval_minus_data <- as.character(mcols(gr[width(gr)==22 & strand(gr)=="-"])$seq)
  flanks_plus_data <- as.character(mcols(gr[width(gr)==22 & strand(gr)=="+"])$with_flanks)
  flanks_minus_data <- as.character(mcols(gr[width(gr)==22 & strand(gr)=="-"])$with_flanks)

  interval_plus <- ggplot() +
    geom_logo(interval_plus_data,
              col_scheme=colors_scheme,
              method = "bits") +
    theme(panel.grid.minor.x = element_blank(),
          plot.title = element_text(lineheight=.5, hjust = 0.5)) +
    ggtitle(label = "+ Strand") +
    #geom_tile(mapping = aes(x = 0:43, y = 2, color = "gray")) +
    #coord_fixed(ratio=5) +
    #facet_grid(.~seq_group,space = "free_x", scales = "free_x") +
    geom_rect(mapping = aes(xmin=0.5, xmax=22+0.5, ymin=0, ymax=2), fill="yellow", alpha=0.1, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 2),
                       breaks = c(0,1,2)) +
    scale_x_continuous(limits = c(0-10, 22+1+10),
                       breaks = pretty_base1(0,22))

  interval_minus <- ggplot() +
    geom_logo(interval_minus_data,
              col_scheme=colors_scheme,
              method = "bits") +
    theme(panel.grid.minor.x = element_blank(),
          #axis.text.y=element_text("RNA Sequence"),
          #axis.title.y=element_blank(),
          #axis.text.y=element_text(),
          #axis.text.y = element_text("RNA Sequence"),
          axis.ticks.y = element_blank(),
          #axis.title.y.right = element_text("RNA Sequence"),
          axis.title.y = element_text(family = "Sans",
                                      size = 12),
          plot.title = element_text(lineheight=.5,
                                    hjust = 0.5,
                                    family = "Sans",
                                    size = 12)) +
    ylab("RNA Sequence") +
    ggtitle(label = "- Strand") +
    #coord_fixed(ratio=5) +
    #facet_grid(.~seq_group,space = "free_x", scales = "free_x") +
    geom_rect(mapping = aes(xmin=0.5, xmax=22+0.5, ymin=0, ymax=2), fill="yellow", alpha=0.1, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 2),
                       breaks = c(0,1,2),
                       position = "right",
                       labels = NULL) +
    scale_x_continuous(limits = c(0-10, 22+1+10),
                       breaks = pretty_base1(0,22))

  flanks_plus <- ggplot() +
    geom_logo(flanks_plus_data,
            col_scheme=colors_scheme,
            method = "bits") +
    theme(panel.grid.minor.x = element_blank()) +
    #facet_grid(.~seq_group) +
    #coord_fixed(ratio=8) +
    geom_rect(mapping = aes(xmin=10.5, xmax=22+10.5, ymin=0, ymax=2), fill="yellow", alpha=0.1, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 2),
                       breaks = c(0,1,2)) +
    scale_x_continuous(limits = c(0, 22+20+1),
                       #breaks = pretty_base1(0,22+20+1, c(0,10), c(1,11)),
                       breaks = pretty_base1(0,22+20+1, c(0,10), c(1,11)),
                       labels = pretty_base1(0-10,22+20+1-10, c(-9, 0), c(-10, 1)))

  flanks_minus <- ggplot() +
    geom_logo(flanks_minus_data,
              col_scheme=colors_scheme,
              method = "bits") +
    theme(panel.grid.minor.x = element_blank(),
          #axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_text(family = "Sans",
                                      size = 12)) +
    ylab("DNA Sequence") +
    #facet_grid(.~seq_group) +
    #coord_fixed(ratio=8) +
    geom_rect(mapping = aes(xmin=10.5, xmax=22+10.5, ymin=0, ymax=2), fill="yellow", alpha=0.1, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 2),
                       breaks = c(0,1,2),
                       position = "right",
                       labels = NULL) +
    scale_x_continuous(limits = c(0, 22+20+1),
                       #breaks = pretty_base1(0,22+20+1, c(0,10), c(1,11)),
                       breaks = pretty_base1(0,22+20+1, c(0,10), c(1,11)),
                       labels = pretty_base1(0-10,22+20+1-10, c(-9, 0), c(-10, 1)))

   ip <- ggplot_gtable(ggplot_build(interval_plus))
   im <- ggplot_gtable(ggplot_build(interval_minus))
   fp <- ggplot_gtable(ggplot_build(flanks_plus))
   fm <- ggplot_gtable(ggplot_build(flanks_minus))
   #c1_label <- rectGrob() + textGrob("column1")
   full <- grid.arrange(textGrob(label = paste(as.character(22), "nt Reads", sep = ""),
                                 gp=gpar(fontfamily = "Sans", size=24), vjust = 0),
                        arrangeGrob(ip, im, ncol = 2),
                        arrangeGrob(fp, fm, ncol = 2),
                        textGrob(label = "Position along RNA or Genome (1 is 5' end)",
                                 gp=gpar(fontfamily = "Sans", size=12), vjust = 0),
                        nrow = 4, ncol = 1, heights = c(.15,1,1,.07))
   #full <- grid.arrange(ip, im, fp, fm, nrow = 2, ncol = 2)
   gtable::gtable_show_layout(full)
   grid.draw(full)
  # maxWidth <- unit.pmax(pgt$widths[2:3], qgt$widths[2:3])
  # pgt$widths[2:3] <- maxWidth
  # qgt$widths[2:3] <- maxWidth
  # #grid.arrange(pgt, qgt, heights = c(7, 3), nrow = 2, ncol = 2)
  # gtable::gtable_show_layout(arrangeGrob(pgt, qgt,
  #                       heights = c(5, 5), #,
  #                       widths = c(6, 6)))
  #                       nrow = 2,
  #                       ncol = 2))
  # gtable::gtable_show_layout(arrangeGrob(pgt, qgt))
  # q + annotate(geom = "text",
  #              x = c(11),
  #              y = c(),
  #              label = c("5'"),
  #              color=c("black"),
  #              size=c(5)) +
  #   annotate(geom = "text",
  #            x = c(33),
  #            y = c(2),
  #            label = c("3'"),
  #            color=c("black"),
  #            size=c(5))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 100)))
  # print(interval_plus, vp = vplayout(1, (20-5):42))
  # print(interval_minus, vp = vplayout(1, 62:(62+22)))
  print(interval_plus, vp = vplayout(1, (10-5):52))
  print(interval_minus, vp = vplayout(1, 52:94))
  print(flanks_plus, vp = vplayout(2, (10-5):52))
  print(flanks_minus, vp = vplayout(2, 52:94))


}
}


