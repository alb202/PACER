##########################################################################################
### The main workflow ###

### 5' plots
p <- five_prime_plot(gr = two_mm, min = NULL, max = NULL)
ggsave(filename = "five_prime_plot__two_mm__all.svg",
       plot = p,
       device = "svg",
       path = dataset_info$figure_dir)
p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm,
                                            regions = genome_data$gene_intervals,
                                            type = "sense",
                                            invert = FALSE))
ggsave(filename = "five_prime_plot__two_mm__sense_to_genes.svg",
       plot = p,
       device = "svg",
       path = dataset_info$figure_dir)
p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm,
                                            regions = genome_data$gene_intervals,
                                            type = "antisense",
                                            invert = FALSE))
ggsave(filename = "five_prime_plot__two_mm__antisense_to_genes.svg",
       plot = p,
       device = "svg",
       path = dataset_info$figure_dir)

### 5' plots - 5' shared with 22nt removed
two_mm_22nt_5prime_filtered <- assign_5prime_to_a_length(gr = two_mm, primary_length = 22)
p <- five_prime_plot(gr = two_mm_22nt_5prime_filtered, min = NULL, max = NULL)
ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__all.svg",
       plot = p,
       device = "svg",
       path = dataset_info$figure_dir)
p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm_22nt_5prime_filtered,
                                            regions = genome_data$gene_intervals,
                                            type = "sense",
                                            invert = FALSE))
ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__sense_to_genes.svg",
       plot = p,
       device = "svg",
       path = dataset_info$figure_dir)
p <- five_prime_plot(gr = filter_by_regions(alignments = two_mm_22nt_5prime_filtered,
                                            regions = genome_data$gene_intervals,
                                            type = "antisense",
                                            invert = FALSE))
ggsave(filename = "five_prime_plot__two_mm_22nt5prime_filtered__antisense_to_genes.svg",
       plot = p,
       device = "svg",
       path = dataset_info$figure_dir)

### 5' plots - 5' shared with 22nt removed
# mut filter_by_gene(gr = genome_data$gene_intervals, gene_list = gene_lists$mut, invert = FALSE, ensembl_id = FALSE)
# csr1 filter_by_gene(gr = genome_data$gene_intervals, gene_list = gene_lists$csr1, invert = FALSE, ensembl_id = TRUE)
# wago1 filter_by_gene(gr = genome_data$gene_intervals, gene_list = gene_lists$wago1, invert = FALSE, ensembl_id = TRUE)


## Make graphs
# 5' by length
p <- five_prime_plot(gr = two_mm)

# 22 vs non-22
df <- count_overlaps_by_width(gr = two_mm, regions = genome_data$gene_intervals, overlap = "antisense", normalized = TRUE)
p <- length_scatter_plot(df = df, comparison_col = "22")

# 22A vs 22C/G/T
df <- count_overlaps_by_width_and_base(gr = two_mm, regions = genome_data$gene_intervals, alignment_width = 22, base_col = "five", overlap =  "antisense", normalized = TRUE )
df <- count_overlaps_by_width_and_base(gr = two_mm, regions = genome_data$exon_intervals, alignment_width = 22, base_col = "five", overlap =  "antisense", normalized = TRUE )
p <- length_scatter_plot(df = df, comparison_col = "G")

two_mm_hm <- calculate_heatmaps(gr = two_mm, length = 22, strand = "+")
two_mm_shuffled_hm <- calculate_heatmaps(gr = two_mm_shuffled, length = 22, strand = "+")
two_mm_bgremoved <- subtract_heatmap_background(gr = two_mm_hm, bg = two_mm_shuffled_hm)

qplot(x = mcols(two_mm[width(two_mm)==22])$five)
qplot(x = mcols(two_mm_shuffled[width(two_mm_shuffled)==22])$five)
