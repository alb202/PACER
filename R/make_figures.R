source("../R/figures.R")

##########################################################################################
### The main workflow ###

## Make plots and save
save_plot <- function(p, path, label){
  print('ggsaving the plot')
  ggsave(filename = paste(label, ".svg", sep = ""),
         plot = p,
         device = "svg",
         path = path)
  return(p)
}




## Make plots and save
make_length_plots <- function(gr, path, label){
#make_length_plots <- function(gr, regions=NULL, overlap="both", invert=FALSE, gene_list=NULL, invert_gl=FALSE, ensembl_id=path, label){
  # if(!is.null(regions) & !is.null(gene_list))
  #   regions <- filter_by_gene(gr = regions, gene_list = gene_list, invert = invert_gl, ensembl_id = ensembl_id)
  # if(!is.null(regions))
  # gr <- filter_by_regions(alignments = gr, regions = regions, type = type, invert = invert)
  print('making the five_prime_plot')
  p <- five_prime_plot(gr = gr)
  print('ggsaving the five_prime_plot')
  ggsave(filename = paste("five_prime_plot__", label, ".svg", sep = ""),
         plot = p,
         device = "svg",
         path = path)
  return(p)
}

make_length_scatter_plots <- function(gr, primary_length=22, regions, overlap="antisense", path, label){
  df <- count_overlaps_by_width(gr = gr,
                                regions = regions,
                                overlap = "antisense",
                                normalized = TRUE)
  p <- scatter_plot(df = df, comparison_col = as.character(primary_length))
  ggsave(filename = paste("scatter_plot__", label, sep = ""),
         plot = p,
         device = "png",
         path = path)
  return(p)
}



#
#
# ### 5' plots
# p <- five_prime_plot(gr = two_mm)
# ggsave(filename = "five_prime_plot__two_mm__all.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm,
#                                             regions = genome_data$gene_intervals,
#                                             type = "sense",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm__sense_to_genes.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm,
#                                             regions = genome_data$gene_intervals,
#                                             type = "antisense",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm__antisense_to_genes.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
#
# ### 5' plots - 5' shared with 22nt removed
# two_mm__22nt_5prime_filtered <- assign_5prime_to_a_length(gr = two_mm, primary_length = 22)
# p <- five_prime_plot(gr = two_mm__22nt_5prime_filtered, min = NULL, max = NULL)
# ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__all.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm__22nt_5prime_filtered,
#                                             regions = genome_data$gene_intervals,
#                                             type = "sense",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__sense_to_genes.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm__22nt_5prime_filtered,
#                                             regions = genome_data$gene_intervals,
#                                             type = "antisense",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__antisense_to_genes.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
#
# ### 5' plots - 5' shared with 22nt removed - MUT targets
# #two_mm__22nt_5prime_filtered <- assign_5prime_to_a_length(gr = two_mm, primary_length = 22)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm__22nt_5prime_filtered,
#                                             regions = filter_by_gene(gr = genome_data$gene_intervals,
#                                                                      gene_list = gene_lists$mut,
#                                                                      invert = FALSE,
#                                                                      ensembl_id = FALSE),
#                                             type = "both",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__MUT_targets.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm__22nt_5prime_filtered,
#                                             regions = filter_by_gene(gr = genome_data$gene_intervals,
#                                                                      gene_list = gene_lists$mut,
#                                                                      invert = FALSE,
#                                                                      ensembl_id = FALSE),
#                                             type = "sense",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__sense_MUT_targets.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
# p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm__22nt_5prime_filtered,
#                                             regions = filter_by_gene(gr = genome_data$gene_intervals,
#                                                                      gene_list = gene_lists$mut,
#                                                                      invert = FALSE,
#                                                                      ensembl_id = FALSE),
#                                             type = "antisense",
#                                             invert = FALSE))
# ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__antisense_MUT_targets.svg",
#        plot = p,
#        device = "svg",
#        path = dataset_info$figure_dir)
#
#
#
# # mut filter_by_gene(gr = genome_data$gene_intervals, gene_list = gene_lists$mut, invert = FALSE, ensembl_id = FALSE)
# # csr1 filter_by_gene(gr = genome_data$gene_intervals, gene_list = gene_lists$csr1, invert = FALSE, ensembl_id = TRUE)
# # wago1 filter_by_gene(gr = genome_data$gene_intervals, gene_list = gene_lists$wago1, invert = FALSE, ensembl_id = TRUE)
#
#
# ## Make graphs
# # 5' by length
# p <- five_prime_plot(gr = two_mm)
#
# # 22 vs non-22
# df <- count_overlaps_by_width(gr = two_mm, regions = genome_data$gene_intervals, overlap = "antisense", normalized = FALSE)
# df_norm <- count_overlaps_by_width(gr = two_mm, regions = genome_data$gene_intervals, overlap = "antisense", normalized = TRUE)
# #p <- scatter_plot(df = df, comparison_col = "22")
# p <- scatter_plot(df = df, x = 22, y = 23)
# p_norm <- scatter_plot(df = df_norm, x = 22, y = 23)
#
# # 22A vs 22C/G/T
# df <- count_overlaps_by_width_and_base(gr = two_mm, regions = genome_data$exon_intervals, alignment_width = 22, base_col = "five", overlap =  "antisense", normalized = TRUE )
# p <- scatter_plot(df = df, x = "G", y = "C")
#
# # Heatmaps
# df_hm_two_mm <- calculate_heatmaps(gr = two_mm, length = 22, strand = "+")
# df_hm_two_mm_shuffled <- calculate_heatmaps(gr = two_mm_shuffled, length = 22, strand = "+")
# df_hm_bgremoved_two_mm <- subtract_heatmap_background(gr = df_hm_two_mm, bg = df_hm_two_mm_shuffled)
# p <- heatmap_plot(heatmap_data = df_hm_two_mm)
# p <- heatmap_plot(heatmap_data = df_hm_two_mm_shuffled)
# p <- heatmap_plot(heatmap_data = df_hm_bgremoved_two_mm)
#
# ### Sequence Logos
# p <- seq_logo_comparisons(gr = two_mm, shuffled_gr = two_mm_shuffled, length = 15, five_prime_base = "G")
#
# # Check to make sure the 5' end are mostly G in the real datasets and random in the shuffled datasets
# qplot(x = mcols(two_mm[width(two_mm)==22])$five)
# qplot(x = mcols(two_mm_shuffled[width(two_mm_shuffled)==22])$five)

# asdf
# asdf
# asdf

