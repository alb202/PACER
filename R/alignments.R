# alignments.R

run_bowtie <- function(bt_options, options_name, genome) {

  alignment_file <- paste(dataset_names["output_dir"],'/',options_name,"_", dataset_names["basename"],".bam", sep = "")
  trimmed_fastq_file <- paste(dataset_names["output_dir"],'/',dataset_names["basename"],".trimmed.fastq", sep = "")
  bowtie_cmd <- make_bowtie_command(
    #output_dir = output_dir, # Directory of trimmed fastq file and sam output
    #input_file = trimmed_fastq, #
    trimmed_fastq_file = trimmed_fastq_file, # Path and name of fastq file
    alignment_file = alignment_file, # Path and name of output alignment file
    bt_options = bt_options, # Aignment options
    #options_name = options_name, # Options name
    #output_name = output_name, # Basename of output
    genome = genome # ID of genome
    )
  #print(bowtie_cmd)
  write(x = "\n", file = paste(dataset_names["output_dir"], "bowtie.log", sep = "/"),
        sep = "", append = TRUE)
  write(x = paste(bt_options, options_name, dataset_names["basename"], sep = " | "),
        file = paste(dataset_names["output_dir"], "bowtie.log", sep = "/"),
        sep = "",
        append = TRUE)

  alignment_output <- system2(command = "bowtie", args = bowtie_cmd, wait = TRUE)
  # write(x = alignment_output,
  #       file = paste(dataset_names["output_dir"], "bowtie.log", sep = "/"),
  #       sep = "\t",
  #       append = TRUE)

  return(alignment_file)
}

make_bowtie_command <- function(trimmed_fastq_file, alignment_file, bt_options, genome){
  #print(input_file)
  #output_sam
  #index_location <- paste(getwd(), "indexes", genome, genome, sep="/")
  #bowtie_cmd <- paste("bowtie", bt_options, genome_indexes[genome], trimmed_fastq_file, "| samtools view -bo", alignment_file,  sep=" ")
  bowtie_cmd <- c(bt_options, genome_indexes[genome], trimmed_fastq_file, "2>>", paste(dataset_names["output_dir"], "/bowtie.log", sep = ""), "| samtools view -bo", alignment_file, "-")
  return(bowtie_cmd)
}
