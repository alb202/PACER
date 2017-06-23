library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)
library(tidyverse)
source("http://bioconductor.org/biocLite.R")
source("R/adapters.R")
source("R/alignments.R")
source("R/filters.R")
source("R/load.R")
source("R/other.R")


#biocLite("SRAdb")
#library(SRAdb)
#library(R.utils)
#library(dplyr)
#library(ShortRead)

#library(rtracklayer)
#library(tibble)
#library(valr)

#library(seqbias)

#library(readr)
#bowtie-build Caenorhabditis_elegans.WBcel235.dna.chromosome.fa ../../indexes/WBcel235/WBcel235

print("settings")
dataset_names <- ""
genome_indexes <- c(ce10 = paste(getwd(), "indexes", "ce10", "ce10", sep="/"),
                    WBcel235 = paste(getwd(), "indexes", "WBcel235", "WBcel235", sep="/"))
adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "/raw_data", sep="")
datasets <- c(WT_early_rep1_full="SRR5023999.fastq.gz")
output_dir <- paste(getwd(), "/output", sep="")
alignment_settings <- c(two_mm="-v2 -k4 --best --strata -S",
               zero_mm="-v0 -k4 --best --strata -S",
               zero_seed_mm="-n0 -e1000 -l22 -k4 --best --strata -S")
genome <- "WBcel235"
genome_files <- tibble(type = c("genome", "sizes", "genes"),
                       WBcel235 = c("Caenorhabditis_elegans.WBcel235.dna.chromosome.fa",
                                    "ce11.chrom.sizes",
                                    "Caenorhabditis_elegans.WBcel235.89.gff3"))
length_range <- c(minimum=10, maximum=30)
sequence_cutoff <- 0.001

for (i in 1:length(datasets)){
  dataset_names <- c(file=NA, # The full name of the file (basename.ext)
                     basename=NA, # The basename of the file (basename)
                     name=NA, # The descriptive name of the file (basename or other)
                     output_dir=NA # The output directory (.../outputdir/basename/)
                     )
  dataset_names["file"] <- datasets[i]
  dataset_names["basename"] <- strsplit(datasets[i], "\\.")[[1]][1]
  #print("1")
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
  #print(dataset_names)
  #print("1.5")
  dataset_names["output_dir"] <- create_output_dirs(dataset_names["name"])
  #print("2")
  # Trim the adapter sequences
  run_cutadapt()
  #print("3")
  #print(dataset_names)
  #dataset_names[trimmed_fastq] <- run_cutadapt()
  #Align the reads using bowtie
  for(j in 1:length(alignment_settings)){
    dataset_names[names(alignment_settings[j])] <- run_bowtie(
      alignment_settings[j],   # The options for this alignment file
      names(alignment_settings[j]), # The name of the sam file configuration
      genome # ID of genome
    )
  #print("4")
  }
  system(command = paste("gzip -f",
                         paste(dataset_names["output_dir"],"/", dataset_names["basename"], ".trimmed.fastq", sep = ""),
                         wait = TRUE)
  )
}

#chrom_sizes <- tbl_genome(x = read.table(file = paste(getwd(), "genomes", genome, genome_files[[2,genome]], sep = "/"), sep = "\t", col.names = c("chrom", "size")))
#gene_intervals <- load_gene_intervals(path = paste(output_dir, genome_files[[3,genome]], sep = "/"))
mart <- get_mart(genome)
chromosome_sizes <- load_genome_data(genome = genome)
gene_intervals <- get_genes_from_biomart(mart = mart)
gene_intervals <- filter_RNA_from_intervals(gene_intervals)
exon_intervals <- get_exons_from_biomart(mart = mart)
exon_intervals <- filter_by_metadata(target = exon_intervals, source = gene_intervals, column = "ensembl_gene_id")
genome_sequence <- Biostrings::readDNAStringSet(filepath = paste("genomes", genome, genome_files$WBcel235[1], sep = "/"), format = "fasta", use.names = TRUE)

two_mismatches <- load_alignments(path = "output/WT_early_rep1_full/two_mm_SRR5023999.bam")
print(length(two_mismatches))
two_mismatches_filtered <- filter_reads_by_regions(alignments = two_mismatches, regions = gene_intervals, type = "both")
print(length(two_mismatches_filtered))
two_mismatches_filtered <- filter_alignments_by_size_range(alignments = two_mismatches_filtered, minimum = length_range["minimum"], maximum = length_range["maximum"])
print(length(two_mismatches_filtered))
two_mismatches_filtered <- remove_overrepresented_sequences(alignments = two_mismatches_filtered, cutoff = sequence_cutoff)
print(length(two_mismatches_filtered))

#two_mismatches_shuffled <- shuffle_intervals(alignments = two_mismatches_filtered, intervals = exon_intervals)
two_mismatches_shuffled <- shuffle_intervals(alignments = two_mismatches_filtered, intervals = exon_intervals)
print(length(two_mismatches_shuffled))
