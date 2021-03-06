library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)
library(tidyverse)
library(Rcpp)
library(data.table)
library(scales)
library(svglite)
library(quantreg)
library(ggseqlogo)
library(geometry)
library(grid)
library(gridExtra)
source("http://bioconductor.org/biocLite.R")
source("R/adapters.R")
source("R/alignments.R")
source("R/annotations.R")
source("R/filters.R")
source("R/load.R")
source("R/other.R")
source("R/sequences.R")
source("R/workflows.R")
source("R/calculations.R")
source("R/figures.R")
source("R/make_figures.R")
source("R/new_genome_data.R")

Rcpp::sourceCpp(file = "cpp/cpp_functions.cpp")

#bowtie-build Caenorhabditis_elegans.WBcel235.dna.chromosome.fa ../../indexes/WBcel235/WBcel235

#### Load settings
genome_indexes <- c(ce10 = paste(getwd(), "data", "indexes", "ce10", "ce10", sep="/"),
                    WBcel235 = paste(getwd(), "data", "indexes", "WBcel235", "WBcel235", sep="/"))
adapter_file <- paste(getwd(), "data/adapters/adapters.txt", sep="")
input_dir <- paste(getwd(), "data/input", sep="")
#datasets <- c(WT_early_rep1="SRR5023999.fastq.gz")
datasets <- c(WT_early_rep1_TEST="SRR5023999_100K_sample.fastq.gz")
output_dir <- paste(getwd(), "/data/output", sep="")
alignment_settings <- c(two_seed_mm="-n2 -e1000 -l22 -k4 --best --strata -S")
genome <- "WBcel235"
genome_files <- list(WBcel235 = list(genome_file="Caenorhabditis_elegans.WBcel235.dna.chromosome.fa",
                                     gene_list_files=list(mut="celegans_mut_ensembl.txt",
                                                     csr1="celegans_csr1_ensembl.txt",
                                                     wago1="celegans_wago1_ensembl.txt")))
# genome_files <- tibble(type = c("genome", "sizes", "genes"),
#                        WBcel235 = c("Caenorhabditis_elegans.WBcel235.dna.chromosome.fa",
#                                     "ce11.chrom.sizes",
#                                     "Caenorhabditis_elegans.WBcel235.89.gff3"))
length_range <- c(minimum=10, maximum=30)
sequence_cutoff <- 0.001

#### Start the workflow
## Set the names of the files that will be created
dataset_info <- set_dataset_names(input_dir = input_dir,
                                  output_dir = output_dir,
                                  dataset_info = datasets)

## Trim the adapter sequences
trimmed_fastq_file <- run_trimmer(adapter_file = adapter_file,
                                   dataset_info = dataset_info)

## Align the reads using bowtie
alignment_file <- run_bowtie(
  alignment_settings[[1]],   # The options for this alignment file
  names(alignment_settings), # The name of the sam file configuration
  genome, # ID of genome
  dataset_info) # Dataset names
## Gzip the fastq file
gzip_a_file(dir = dataset_info["output_dir"], file = trimmed_fastq_file)

#### Import other genomic data
## Get genomic features

if(!check_for_intervals(path = "data/genomes", genome = "WBcel235")){
  get_ensembl_intervals(path = "data/genomes", genome = "WBcel235")
}
genome_data <- load_genome_data(path = "data/genomes", genome = "WBcel235")
#genome_data <- load_genome_data(path = paste(getwd(),"data/genomes", sep="/"), genome = genome)
genome_sequence <- load_fasta_genome(path = paste(getwd(),"genomes",genome,as.character(genome_files[[genome]]$genome_file), sep="/"))
gene_lists <- load_gene_lists(path = paste(getwd(),"genomes",genome, sep = "/"), gene_lists_files = genome_files[[genome]]$gene_list_files)
## Load and filter the main alignment file
#alignments <- load_alignments(path = alignment_file)
alignments <- load_alignments(path = "data/output/WT_early_rep1/two_seed_mm_SRR5023999.bam")
alignments <- filter_alignments(alignments = alignments,
                  regions = genome_data[["gene_intervals"]],
                  regions_filter = "both",
                  minimum = length_range["minimum"],
                  maximum = length_range["maximum"],
                  cutoff = sequence_cutoff)

mismatch_indexes <- filter_BAM_tags(alignments)
values$two_mm <- alignments[mm_indexes$two_mm]
values$no_mm <- alignments[mm_indexes$no_mm]
values$no_mm_in_seed <- alignments[mm_indexes$no_mm_seed]


## Create a shuffled version of the file
two_mm_shuffled <- shuffle_intervals(alignments = two_mm,
                                             intervals = genome_data[["exon_intervals"]],
                                             antisense = TRUE)

## Get the sequences for the actual and the shuffled alignments
two_mm <- get_genome_sequence(
  gr = two_mm, genome_sequence = genome_sequence)
two_mm_shuffled <- get_genome_sequence(
  gr = two_mm_shuffled, genome_sequence = genome_sequence)

no_mm <- get_genome_sequence(
  gr = no_mm, genome_sequence = genome_sequence)
