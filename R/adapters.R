
make_adapter_string <- function(adapter_file){
  adapter_list <- read.delim(adapter_file, sep = " ", header = FALSE, comment.char = "#")[,1]
  if (length(adapter_list)==0)
    stop("No adapter sequences were found")
  adapter_string <- ""
  for (i in adapter_list){
    adapter_string <- paste(adapter_string, "-a", i, sep = " ")
  }
  return(adapter_string)
}

make_cutadapt_command <- function(trimmed_fastq_file, adapter_file, dataset_info){
  adapter_string <- make_adapter_string(adapter_file)
  file_path <- paste(input_dir, "/", dataset_info["fastq_file"], sep="")
  cutadapt_cmd <- paste("cutadapt","-m10", adapter_string, "-o", paste(dataset_info["output_dir"],'/',trimmed_fastq_file, sep = ""), file_path)
  return(cutadapt_cmd)
}

run_cutadapt <- function(adapter_file, dataset_info){
  trimmed_fastq_file <- paste(dataset_info["basename"],
                              ".trimmed.fastq", sep = "")
  cutadapt_cmd <- make_cutadapt_command(trimmed_fastq_file, adapter_file, dataset_info)
  cutadapt_output <- system(command = cutadapt_cmd,
                            wait = TRUE,
                            intern = TRUE)
  write(x = cutadapt_output,
        file = paste(dataset_info["output_dir"],
                     "cutadapt_log.txt", sep = "/"),
        sep = "\t" )
  return(trimmed_fastq_file)
}
