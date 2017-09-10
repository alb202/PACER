# # Create the adapter part of the cutadapt command
# make_adapter_string <- function(adapter_file){
#   adapter_list <- read.delim(file = adapter_file,
#                                             strip.white = TRUE,
#                                             sep = " ",
#                                             header = FALSE,
#                                             comment.char = "#",
#                                             stringsAsFactors = FALSE)[,1]
#   if (length(adapter_list)==0)
#     stop("No adapter sequences were found")
# return(paste(" -a", adapter_list, collapse = ""))
# }

# Create the full cutadapt command
make_trim_command <- function(input_dir, output_dir, input_file, dataset_ID, adapter_file, min_length){
  adapter_string <- paste(" -a", read.delim(file = adapter_file,
                                            strip.white = TRUE,
                                            sep = " ",
                                            header = FALSE,
                                            comment.char = "#",
                                            stringsAsFactors = FALSE)[,1],
                          collapse = "")
  file_path <- getAbsolutePath(paste(input_dir, input_file, sep="/"))
  output_file <- paste(output_dir,"/",dataset_ID,
                       ".trimmed.fastq",
                       sep = "")
  trimmer_cmd <- paste("cutadapt",
                        "--max-n=0",
                        paste("-m", min_length, sep = ""),
                        adapter_string,
                        "-o",
                        output_file,
                        file_path, sep = " ")
  return(trimmer_cmd)
}

run_trimmer <- function(output_dir, dataset_ID, trim_cmd){
  # trimmed_fastq_file <- paste(dataset_info["basename"],
  #                             ".trimmed.fastq", sep = "")
  #cutadapt_cmd <- make_cutadapt_command(trimmed_fastq_file, adapter_file, dataset_info)
  trim_output <- system(command = trim_cmd,
                            wait = TRUE,
                        intern = TRUE)
  #return(trim_output)#, stringsAsFactors = FALSE))
  write(x = trim_output,
        file = paste(output_dir,
                     "/",
                     dataset_ID,
                     ".trim.log",
                     sep = ""),
        sep = "\t" )
  #return(trimmed_fastq_file)
}
