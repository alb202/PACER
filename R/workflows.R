# line_plots <- function()






filter_alignments <- function(alignments, regions, regions_filter="both", minimum=10, maximum=30, cutoff=.001){
  alignments <- filter_by_regions(alignments = alignments, regions = regions, type = regions_filter)
  alignments <- filter_alignments_by_size_range(alignments = alignments, minimum = minimum, maximum = maximum)
  alignments <- remove_overrepresented_sequences(alignments = alignments, cutoff = cutoff)
  return(sort.GenomicRanges(alignments))
}

set_dataset_names <- function(input_dir, output_dir, dataset_info){
  dataset_names <- list(fastq_file=NA, # The full name of the file (basename.ext)
                     basename=NA, # The basename of the file (basename)
                     name=NA, # The descriptive name of the file (basename or other)
                     output_dir=NA # The output directory (.../outputdir/basename/)
  )
  dataset_names["fastq_file"] <- dataset_info[[1]]
  dataset_names["basename"] <- strsplit(dataset_info, "\\.")[[1]][1]
  if (is.na(names(dataset_info))){
    dataset_names["name"] <- dataset_names[["basename"]]
  } else{
      dataset_names["name"] <- names(dataset_info[1])}
  dataset_names["output_dir"] <- create_output_dirs(out_dir = output_dir, name = dataset_names[["name"]])
  dataset_names["figure_dir"] <- create_output_dirs(out_dir = dataset_names[["output_dir"]], name =c("figures"))
  return(dataset_names)
}


#
# main_workflow <- function()
# {
#   incProgress(amount = .2, detail = "Loading intervals")
#   values$genome_data <- load_genome_data(path = values$genomes_dir,
#                                          genome = values$selected_genome[["Version"]])
#   print(values$genome_data)
#   print(values$selected_genome[["Gene.Sets"]])
#   incProgress(amount = .2, detail = "Loading gene sets")
#   values$gene_sets <- load_gene_sets(gene_sets = values$selected_genome[["Gene.Sets"]])
#
#   ### Load alignment
#   incProgress(amount = .2, detail = "Loading alignments")
#   values$alignments <- load_alignments(path = values$bam_path)
#
#   incProgress(amount = .2, detail = "Filtering regions")
#   values$alignments <- filter_by_regions(alignments = values$alignments,
#                                          regions = values$genome_data[["gene_intervals"]],
#                                          type = "both")
#   incProgress(amount = .1, detail = "Filtering alignment sizes")
#   values$alignments <- filter_alignments_by_size(alignments = values$alignments,
#                                                  minimum = input$get_range[[1]],
#                                                  maximum = input$get_range[[2]])
#   incProgress(amount = .1, detail = "Filtering overrpreseented reads")
#   values$alignments <- remove_overrepresented_sequences(alignments = values$alignments,
#                                                         cutoff = input$read_cutoff)
#
#   values$alignments <- get_genome_sequence(gr = values$alignments,
#                                            genome_sequence = load_fasta_genome(
#                                              path = values$selected_genome[['Genome.FASTA']]
#                                            ))
#
#   print('print(values$alignments)')
#   print(values$alignments)
#
#   incProgress(amount = .1, detail = "Filtering mismatches")
#   mismatch_indexes <- filter_BAM_tags(values$alignments)
#   values$two_mm <- values$alignments[mismatch_indexes$two_mm]
#   values$no_mm <- values$alignments[mismatch_indexes$no_mm]
#   values$no_mm_in_seed <- values$alignments[mismatch_indexes$no_mm_seed]
#   values$shuffled <- shuffle_alignments(alignments = values$two_mm,
#                                         intervals = values$genome_data[["gene_intervals"]],
#                                         antisense = TRUE)
#   #values$two_mm <- filter_alignments(alignments = values$two_mm, regions = )
#
#   print('The length of values$two_mm is: ')
#   print(length(x = print(values$two_mm)))
#   #print(values$two_mm)
#   #print(values$no_mm)
#   #print(values$no_mm_in_seed)
#   print('making the plots')
#   print('making the plots two_mm')
#   print(width(values$two_mm))
#   print(strand(values$two_mm))
#   print(mcols(values$two_mm))
#
#   create_output_dirs(out_dir = values$output_dir,
#                      name = 'five_prime')
#   p <- make_length_plots(gr = values$two_mm,
#                          path = paste(values$output_dir,
#                                       '/five_prime',
#                                       sep = ''),
#                          label = "two_mm__all_genes")
#   output$fivepp_all <- renderPlot({p})
#   print('making the plots no_mm')
#   make_length_plots(gr = values$no_mm,
#                     path = paste(values$output_dir,
#                                  '/five_prime',
#                                  sep = ''),
#                     label = "no_mm__all_genes")
#   print('making the plots no_mm_in_seed')
#   make_length_plots(gr = values$no_mm_in_seed,
#                     path = paste(values$output_dir,
#                                  '/five_prime',
#                                  sep = ''),
#                     label = "no_mm_in_seed__all_genes")
# }


# process_fastq <- function(input_dir, output_dir, dataset_info, adapter_file, genome, alignment_settings){
#
#   # Trim the adapter sequences
#   run_cutadapt()
#
#   #Align the reads using bowtie
#   for(j in 1:length(alignment_settings)){
#     dataset_names[names(alignment_settings[j])] <- run_bowtie(
#       alignment_settings[j],   # The options for this alignment file
#       names(alignment_settings[j]), # The name of the sam file configuration
#       genome # ID of genome
#     )
#   }
#   system(command = paste("gzip -f",
#                          paste(dataset_names["output_dir"],"/", dataset_names["basename"], ".trimmed.fastq", sep = ""),
#                          wait = TRUE)
#   )
# return(dataset_names)
# }


