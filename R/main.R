#source("http://bioconductor.org/biocLite.R")
source("R/adapters.R")
source("R/alignments.R")
source("R/settings.R")
#biocLite("SRAdb")
#library(SRAdb)
library(R.utils)
#library(ShortRead)

print("settings")
dataset_names <- ""
genome_indexes <- c(ce10 = paste(getwd(), "indexes", "ce10", "ce10", sep="/"))
adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1="SRR5023999_100K_sample.fastq.gz")
output_dir <- paste(getwd(), "/output", sep="")
sam_files <- c(two_mm="-v2 -k4 --best -S",
               zero_mm="-v0 -k4 --best -S",
               zero_seed_mm="-n0 -e1000 -l22 -k4 --best -S")
genome = "ce10"

for (i in 1:length(datasets)){
  dataset_names <- c(file=NA, # The full name of the file (basename.ext)
                     basename=NA, # The basename of the file (basename)
                     name=NA, # The descriptive name of the file (basename or other)
                     output_dir=NA # The output directory (.../outputdir/basename/)
                     )
  dataset_names["file"] <- datasets[i]

  dataset_names["basename"] <- strsplit(datasets[i], "\\.")[[1]][1]
  print("1")
  if (is.na(names(datasets[i])))
    dataset_names["name"] <- dataset_names["basename"]
  else
    dataset_names["name"] <- names(datasets[i])
  #dataset_names["name"] <- names(datasets[i])
  #if (dataset_names["name"] == "")
  #  dataset_names["name"] <- dataset_names["basename"]
  #localfile <- file.path(input_dir, file)
  #print(localfile)
  #print(filename)
  #print(name)
  print(dataset_names)
  print("1.5")
  dataset_names["output_dir"] <- create_output_dirs(dataset_names["name"])
  print("2")
  # Trim the adapter sequences
  run_cutadapt()
  print("3")
  print(dataset_names)
  #dataset_names[trimmed_fastq] <- run_cutadapt()
  #Align the reads using bowtie
  for(j in 1:length(sam_files)){
    dataset_names[names(sam_files[j])] <- run_bowtie(
      #dataset_names["output_dir"],    # The directory of the trimmed fastq
      #paste(dataset_names[basename],".trimmed.fastq", sep=""), # The trimmed fastq file name
      sam_files[j],   # The options for this alignment file
      names(sam_files[j]), # The name of the sam file configuration
      dataset_names["basename"], # The basename of the file
      genome # ID of genome
    )
  print("4")
  }
  system(command = paste("gzip -f",
                         paste(dataset_names["output_dir"],"/", dataset_names["basename"], ".trimmed.fastq", sep = ""),
                         wait = TRUE)
  )

}
