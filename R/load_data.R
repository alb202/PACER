source("http://bioconductor.org/biocLite.R")
biocLite("SRAdb")
#library(SRAdb)
library(R.utils)
library(ShortRead)

input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1="SRR5023999.fastq.gz")

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
  fq <- readFastq(localfile)
  print(class(fq))
  print(object.size(fq))
  print(fq)


  }
