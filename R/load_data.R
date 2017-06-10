source("http://bioconductor.org/biocLite.R")
#biocLite("SRAdb")
#library(SRAdb)
library(R.utils)
library(ShortRead)

adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1="SRR5023999_100K_sample.fastq.gz")
output_dir <- paste(getwd(), "/output", sep="")
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
  new_filename <- paste(input_dir,'/',filename,".trimmed.fastq.gz", sep = "")
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

make_bowtie_command <- function(input_file){
  print(input_file)
  new_filename <- paste(output_dir,'/',filename,".sam", sep = "")
  index_location <- paste(getwd(), "indexes", genome, genome, sep="/")
  bowtie_cmd <- paste("gzip -d -c", input_file, "|", "bowtie", "-v2 -k4 --best -S", index_location, "-", new_filename, "2>&1", sep=" ")
  #gzip -d -c SRR5023999_100K_sample.trimmed.fastq.gz | bowtie -v2 -k4 --best -S ../indexes/ce10/ce10  "-" SRR5023999_100K_sample.trimmed.sam
  return(c(bowtie_cmd, new_filename))
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
  cutadapt_cmd <- make_cutadapt_command()
  cutadapt_output <- system(command = cutadapt_cmd,
                            wait = TRUE,
                            intern = TRUE)
  write(x = cutadapt_output,
        file = paste(output_dir, "cutadapt_log.txt", sep = "/"), sep = "\t" )

  #Align the reads using bowtie
  bowtie_cmd <- make_bowtie_command(cutadapt_cmd[2])
  print(bowtie_cmd)
  bowtie_output <- system(command = bowtie_cmd,
                            wait = TRUE,
                            intern = TRUE)
  write(x = bowtie_output,
        file = paste(output_dir, "bowtie_log.txt", sep = "/"), sep = "\t" )
  }

