# alignments.R
run_bowtie <- function(bt_options, options_name, genome, dataset_info) {

  alignment_file <- paste(dataset_info["output_dir"],'/',options_name,"_", dataset_info["basename"],".bam", sep = "")
  trimmed_fastq_file <- paste(dataset_info["output_dir"],'/',dataset_info["basename"],".trimmed.fastq", sep = "")
  bowtie_cmd <- make_bowtie_cmd(
    trimmed_fastq_file = trimmed_fastq_file, # Path and name of fastq file
    alignment_file = alignment_file, # Path and name of output alignment file
    bt_options = bt_options, # Aignment options
    genome = genome, # ID of genome
    dataset_info = dataset_info) # The dataset names

  write(x = "\n", file = paste(dataset_info["output_dir"], "bowtie.log", sep = "/"),
        sep = "", append = TRUE)
  write(x = paste(bt_options, options_name, dataset_info["basename"], sep = " | "),
        file = paste(dataset_info["output_dir"], "bowtie.log", sep = "/"),
        sep = "",
        append = TRUE)

  alignment_output <- system2(command = "bowtie", args = bowtie_cmd, wait = TRUE)
  return(alignment_file)
}

make_bowtie_cmd <- function(trimmed_fastq_file, alignment_file, bt_options, genome, dataset_info){
  bowtie_cmd <- c(bt_options,
                  genome_indexes[genome],
                  trimmed_fastq_file, "2>>",
                  paste(dataset_info["output_dir"],
                        "/bowtie.log",
                        sep = ""),
                  "| samtools sort - | samtools view -bo",
                  alignment_file,
                  "-")
  return(bowtie_cmd)
}

