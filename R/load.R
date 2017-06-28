
load_gene_intervals <- function(path){
 results <- import(con = "genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.gff3")
 return(sort.GenomicRanges(results))
}

filter_RNA_from_intervals <- function(intervals){
  results <- subset(x = intervals,
                    gene_biotype!="snoRNA" &
                      gene_biotype!="miRNA" &
                      gene_biotype!="rRNA" &
                      gene_biotype!="tRNA" &
                      gene_biotype!="snRNA")
  return(sort.GenomicRanges(results))
}

load_alignments <- function(path, params=ScanBamParam(reverseComplement = FALSE, what = c("seq", "qname", "flag"))){
  results <- readGAlignments(file = path, param = params)
  return(sort.GenomicRanges(results))
}


load_genome_data <- function(genome){
  marts <- listMarts()
  useMart(marts[[1,1]])
  mart_datasets <- listDatasets(useMart(biomart=marts[[1,1]], host="www.ensembl.org"))
  dataset_id <- mart_datasets[mart_datasets$version==genome,][[1]]
  chromosomes <- getChromInfoFromBiomart(biomart=marts[[1,1]],
                                       dataset=dataset_id)
  return(chromosomes)
}

# get exons by gene
# txdb <- makeTxDbFromBiomart(biomart="ensembl",dataset="celegans_gene_ensembl")

# library(ensembldb)
# library(AnnotationHub)
# source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

# marts <- listMarts()
# useMart(marts[[1,1]])
# GenomicFeatures::getChromInfoFromBiomart(dataset = dataset_id)
# chrom_sizes <- GenomicFeatures::getChromInfoFromBiomart(dataset = dataset_id)
# um <- useMart(marts[[1,1]], dataset = c(dataset_id))
# biomaRt::useEnsembl("ensembl", a)
# biomaRt::listFilters(um)
# bm_exons <- biomaRt::getBM(
#   attributes = c("external_gene_name", "external_gene_source","ensembl_gene_id",
#                  "ensembl_exon_id", "strand", "start_position", "end_position",
#                  "exon_chrom_start", "exon_chrom_end","rank",
#                  "5_utr_start" ,"5_utr_end" ,"3_utr_start" ,"3_utr_end"), mart = um)
#
# bm_genes <- biomaRt::getBM(
#   attributes = c("description", "external_gene_name", "external_gene_source",
#                  "ensembl_gene_id", "chromosome_name", "strand", "start_position",
#                  "end_position", "source", "status", "gene_biotype"), mart = um)
# bm_genes$strand <- replace(x = bm_genes$strand, list = bm_genes$strand=="-1", "-")
# bm_genes$strand <- replace(x = bm_genes$strand, list = bm_genes$strand=="1", "+")
# GenomicRanges::makeGRangesFromDataFrame(
#   df=genes, seqnames.field = "chromosome_name", start.field = "start_position",
#   end.field = "end_position", strand.field = "strand", keep.extra.columns = TRUE)
# useEnsembl(biomart = c(marts[[1,1]]))

get_mart <- function(genome){
  mart_name <- c(listMarts()[[1,1]])
  mart <- useMart(biomart = mart_name)
  datasets_list <- listDatasets(mart = mart)
  dataset <- datasets_list[datasets_list["version"]==genome][1]
  mart <- useDataset(dataset = dataset, mart = mart)
  return(mart)
}

get_genes_from_biomart <- function(mart){
  genes <- biomaRt::getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "strand",
                   "description", "external_gene_name", "external_gene_source",
                   "ensembl_gene_id", "source", "status", "gene_biotype"), mart = mart)
  genes$strand <- swap_values(x = genes$strand, old = c(1, -1), new = c("+", "-"))
  new_grange <- makeGRangesFromDataFrame(
    df=genes, seqnames.field = "chromosome_name", start.field = "start_position",
    end.field = "end_position", strand.field = "strand", keep.extra.columns = TRUE)
  return(sort.GenomicRanges(new_grange))
}

get_exons_from_biomart <- function(mart){
  exons <- exons <- biomaRt::getBM(
    attributes = c("strand", "exon_chrom_start", "exon_chrom_end", "chromosome_name",
                   "external_gene_name", "external_gene_source", "ensembl_gene_id",
                   "ensembl_exon_id", "rank"), mart = mart)
  exons$strand <- swap_values(x = exons$strand, old = c(1, -1), new = c("+", "-"))
  new_grange <- makeGRangesFromDataFrame(
    df=exons, seqnames.field = "chromosome_name", start.field = "exon_chrom_start",
    end.field = "exon_chrom_end", strand.field = "strand", keep.extra.columns = TRUE)
  return(sort.GenomicRanges(new_grange))
}

filter_by_metadata <- function(target, source, column){

  matches <- mcols(target)[,column] %in% mcols(source)[,column]
  results <- target[matches]
  return(sort.GenomicRanges(results))
}

load_fasta_genome <- function(path){
  #TRy importing using DNAString() then converting to data.table
  #data.table::fread(input = "genomes/WBcel235/Caenorhabditis_elegans.WBcel235.dna.chromosome.TEST.fa")
  genome_fasta <- Biostrings::readDNAStringSet(filepath = path, format = "fasta", use.names = TRUE)
  genome_sequence <- list()
  #genome_fasta[strsplit(x = names(a)[1], split = " " )[[1]][1]] <- as.character(a[names(a)[[1]]])
  for (i in 1:length(genome_fasta)){
    genome_sequence[strsplit(x = names(genome_fasta)[i], split = " " )[[1]][1]] <-
      strsplit(x = as.character(genome_fasta[i]), split = '')
  }
  return(genome_sequence)
}

load_fasta_genome2 <- function(path){
  #TRy importing using DNAString() then converting to data.table
  #data.table::fread(input = "genomes/WBcel235/Caenorhabditis_elegans.WBcel235.dna.chromosome.TEST.fa")
  genome_fasta <- Biostrings::readDNAStringSet(filepath = path, format = "fasta", use.names = TRUE)
  genome_sequence <- list()
  #genome_fasta[strsplit(x = names(a)[1], split = " " )[[1]][1]] <- as.character(a[names(a)[[1]]])
  for (i in 1:length(genome_fasta)){
    genome_sequence[strsplit(x = names(genome_fasta)[i], split = " " )[[1]][1]] <-
      unlist(as.character(genome_fasta[i]))
  }
  return(genome_sequence)
}
