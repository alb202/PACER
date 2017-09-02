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


figure_workflow <- function(dataset1, dataset2, dataset3, dataset_shuffled, figure_dir){

}




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


