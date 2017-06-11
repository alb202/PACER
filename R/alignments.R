# alignments.R

run_bowtie <- function(trimmed_fastq, bt_options, output_name) {
  bowtie_cmd <- make_bowtie_command(input_file = trimmed_fastq,
                                    bt_options = bt_options,
                                    output_name = output_name)
  print(bowtie_cmd)
  bowtie_output <- system(command = bowtie_cmd,
                          wait = TRUE,
                          intern = TRUE)

  write(x = paste(bt_options, output_name, sep = " | "),
        file = paste(output_dir, "bowtie_log.txt", sep = "/"),
        sep = "\t",
        append = TRUE)
  write(x = bowtie_output,
        file = paste(output_dir, "bowtie_log.txt", sep = "/"),
        sep = "\t",
        append = TRUE)
}

make_bowtie_command <- function(input_file, bt_options, output_name ){
  print(input_file)
  new_filename <- paste(output_dir,'/',output_name,"_",filename,".sam", sep = "")
  index_location <- paste(getwd(), "indexes", genome, genome, sep="/")
  bowtie_cmd <- paste("bowtie", bt_options, index_location, input_file, new_filename, "2>&1", sep=" ")
  return(c(bowtie_cmd, new_filename))
}
