source("http://bioconductor.org/biocLite.R")
#biocLite("SRAdb")
#library(SRAdb)
library(R.utils)
library(ShortRead)

adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1="SRR5023999_100K_sample.fastq.gz")
output_dir <- paste(getwd(), "/output", sep="")
sam_files <- c(two_mismatch="-v2 -k4 --best -S",
               no_mismatch="-v0 -k4 --best -S",
               no_seed_mismatch="-n0 -e1000 -l22 -k4 --best -S")
genome = "ce10"

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


make_cutadapt_command <- function(){
  adapter_string <- make_adapter_string()
  new_filename <- paste(output_dir,'/',filename,".trimmed.fastq", sep = "")
  cutadapt_cmd <- paste("cutadapt","-m10",adapter_string, "-o", new_filename, localfile)
  return(c(cutadapt_cmd, new_filename))
}

create_output_dirs <- function(){
  if(!dir.exists(output_dir))
    dir.create(output_dir)
  if(!dir.exists(paste(output_dir, name, sep = "/")))
    dir.create(paste(output_dir, name, sep = "/"))
  return(paste(output_dir, name, sep = "/"))
}

make_bowtie_command <- function(input_file, bt_options, output_name ){
  print(input_file)
  new_filename <- paste(output_dir,'/',output_name,"_",filename,".sam", sep = "")
  index_location <- paste(getwd(), "indexes", genome, genome, sep="/")
  bowtie_cmd <- paste("bowtie", bt_options, index_location, input_file, new_filename, "2>&1", sep=" ")
  return(c(bowtie_cmd, new_filename))
}

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
run_cutadapt <- function(){
  cutadapt_cmd <- make_cutadapt_command()
  cutadapt_output <- system(command = cutadapt_cmd,
                            wait = TRUE,
                            intern = TRUE)
  write(x = cutadapt_output,
        file = paste(output_dir, "cutadapt_log.txt", sep = "/"), sep = "\t" )
  return(cutadapt_cmd[2])
}

for(i in 1:length(datasets)){
  file <- datasets[i]
  name <- names(datasets[i])
  filename <- strsplit(datasets[i], "\\.")[[1]][1]
  if (name=="")
    name <- filename
  localfile <- file.path(input_dir, file)
  print(localfile)
  print(filename)
  print(name)
  output_dir <- create_output_dirs()
  # Trim the adapter sequences
  trimmed_fastq <- run_cutadapt()
  #Align the reads using bowtie
  for(i in 1:length(sam_files)){
    run_bowtie(trimmed_fastq, sam_files[i], names(sam_files[i]))
  }
  system(command = paste("gzip -f", trimmed_fastq, sep = " "), wait = TRUE)

}


