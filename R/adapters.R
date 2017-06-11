run_cutadapt <- function(){
  cutadapt_cmd <- make_cutadapt_command()
  cutadapt_output <- system(command = cutadapt_cmd,
                            wait = TRUE,
                            intern = TRUE)
  write(x = cutadapt_output,
        file = paste(output_dir, "cutadapt_log.txt", sep = "/"), sep = "\t" )
  return(cutadapt_cmd[2])
}

make_cutadapt_command <- function(){
  adapter_string <- make_adapter_string()
  new_filename <- paste(output_dir,'/',filename,".trimmed.fastq", sep = "")
  cutadapt_cmd <- paste("cutadapt","-m10",adapter_string, "-o", new_filename, localfile)
  return(c(cutadapt_cmd, new_filename))
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
