source("http://bioconductor.org/biocLite.R")
source("R/adapters.R")
source("R/alignments.R")
source("R/settings.R")
#biocLite("SRAdb")
#library(SRAdb)
library(R.utils)
library(ShortRead)

print("settings")
adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1="SRR5023999_100K_sample.fastq.gz")
output_dir <- paste(getwd(), "/output", sep="")
sam_files <- c(two_mismatch="-v2 -k4 --best -S",
               no_mismatch="-v0 -k4 --best -S",
               no_seed_mismatch="-n0 -e1000 -l22 -k4 --best -S")
genome = "ce10"

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

# Other helper functions below
create_output_dirs <- function(){
  if(!dir.exists(output_dir))
    dir.create(output_dir)
  if(!dir.exists(paste(output_dir, name, sep = "/")))
    dir.create(paste(output_dir, name, sep = "/"))
  return(paste(output_dir, name, sep = "/"))
}
