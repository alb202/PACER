run_cutadapt <- function(){
  trimmed_output <- paste(dataset_names["basename"],
                          ".trimmed.fastq", sep = "")
  cutadapt_cmd <- make_cutadapt_command(trimmed_output)
  cutadapt_output <- system(command = cutadapt_cmd,
                            wait = TRUE,
                            intern = TRUE)
  write(x = cutadapt_output,
        file = paste(dataset_names["output_dir"],
                     "cutadapt_log.txt", sep = "/"),
        sep = "\t" )
  return(trimmed_output)
}

make_cutadapt_command <- function(trimmed_output){
  adapter_string <- make_adapter_string()
  file_path <- paste(input_dir, "/", dataset_names["file"], sep="")
  cutadapt_cmd <- paste("cutadapt","-m10", adapter_string, "-o", paste(dataset_names["output_dir"],'/',trimmed_output, sep = ""), file_path)
  return(cutadapt_cmd)
}

make_adapter_string <- function(){
  adapter_list <- read.delim(adapter_file, sep = " ", header = FALSE, comment.char = "#")[,1]
  if (length(adapter_list)==0)
    stop("No adapter sequences were found")
  adapter_string <- ""
  for (i in adapter_list){
    adapter_string <- paste(adapter_string, "-a", i, sep = " ")
  }
  return(adapter_string)
}
