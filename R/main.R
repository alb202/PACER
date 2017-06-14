#source("http://bioconductor.org/biocLite.R")
source("R/adapters.R")
source("R/alignments.R")
source("R/other.R")
#biocLite("SRAdb")
#library(SRAdb)
library(R.utils)
library(ShortRead)
library(GenomicAlignments)
#bowtie-build Caenorhabditis_elegans.WBcel235.dna.chromosome.fa ../../indexes/WBcel235/WBcel235

print("settings")
dataset_names <- ""
genome_indexes <- c(ce10 = paste(getwd(), "indexes", "ce10", "ce10", sep="/"),
                    WBcel235 = paste(getwd(), "indexes", "WBcel235", "WBcel235", sep="/"))
adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1_full="SRR5023999.fastq.gz")
output_dir <- paste(getwd(), "/output", sep="")
bam_files <- c(two_mm="-v2 -k4 --best --strata -S",
               zero_mm="-v0 -k4 --best --strata -S",
               zero_seed_mm="-n0 -e1000 -l22 -k4 --best --strata -S")
genome <- "WBcel235"
genome_files <- tibble(type = c("genome", "sizes", "genes", "filter", "filter", "filter", "filter"),
                       WBcel235 = c("Caenorhabditis_elegans.WBcel235.dna.chromosome.fa",
                                    "ce11.chrom.sizes","Caenorhabditis_elegans.WBcel235.89.gff3",
                                    "Caenorhabditis_elegans.WBcel235.89.snoRNA.gff3",
                                    "Caenorhabditis_elegans.WBcel235.89.rRNA.gff3",
                                    "Caenorhabditis_elegans.WBcel235.89.tRNA.gff3",
                                    "Caenorhabditis_elegans.WBcel235.89.miRNA.gff3"))
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
  for(j in 1:length(bam_files)){
    dataset_names[names(bam_files[j])] <- run_bowtie(
      bam_files[j],   # The options for this alignment file
      names(bam_files[j]), # The name of the sam file configuration
      genome # ID of genome
    )
  print("4")
  }
  system(command = paste("gzip -f",
                         paste(dataset_names["output_dir"],"/", dataset_names["basename"], ".trimmed.fastq", sep = ""),
                         wait = TRUE)
  )


}

