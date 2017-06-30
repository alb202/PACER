get_exons_from_biomart <- function(mart){
  exons <- biomaRt::getBM(
    attributes = c("strand", "exon_chrom_start", "exon_chrom_end", "chromosome_name",
                   "external_gene_name", "external_gene_source", "ensembl_gene_id",
                   "ensembl_exon_id", "rank"), mart = mart)
  exons$strand <- swap_values(x = exons$strand, old = c(1, -1), new = c("+", "-"))
  new_grange <- makeGRangesFromDataFrame(
    df=exons, seqnames.field = "chromosome_name", start.field = "exon_chrom_start",
    end.field = "exon_chrom_end", strand.field = "strand", keep.extra.columns = TRUE)
  return(sort.GenomicRanges(new_grange))
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

get_mart <- function(genome){
  mart_name <- c(listMarts()[[1,1]])
  mart <- useMart(biomart = mart_name)
  datasets_list <- listDatasets(mart = mart)
  dataset <- datasets_list[datasets_list["version"]==genome][1]
  mart <- useDataset(dataset = dataset, mart = mart)
  return(mart)
}

load_alignments <- function(path, params=ScanBamParam(reverseComplement = FALSE, what = c("seq", "qname", "flag"), tag = c("XA", "MD", "NM"))){
  return(sort.GenomicRanges(readGAlignments(file = path, param = params)))
}

load_chromosome_sizes <- function(genome){
  mart_name <- c(listMarts()[[1,1]])
  mart_datasets <- listDatasets(useMart(biomart=mart_name, host="www.ensembl.org"))
  dataset_id <- mart_datasets[mart_datasets$version==genome,][[1]]
  #chromosomes <- getChromInfoFromBiomart(biomart=mart_name,
  #                                     dataset=dataset_id)
  return(getChromInfoFromBiomart(biomart=mart_name, dataset=dataset_id))
}

load_fasta_genome <- function(path){
  # Load the genome as a Biostrings object
  genome_fasta <- Biostrings::readDNAStringSet(filepath = path, format = "fasta", use.names = TRUE)
  # Create an empty list to store the sequences of each chromosome
  genome_sequence <- list()
  # Loop over each DNAstring, split the chromosome name, and save it in the list
  for (i in 1:length(genome_fasta)){
    genome_sequence[strsplit(x = names(genome_fasta)[i], split = " " )[[1]][1]] <-
      unlist(as.character(genome_fasta[i]))
  }
  return(genome_sequence)
}

load_genome_data <- function(genome){
  genomes_path <- paste(getwd(),"genomes",genome, sep="/")
  files <- dir(path = genomes_path)

  gene_file_name <- paste(genome,"_gene_intervals.tsv", sep = "")
  if(gene_file_name %in% files){
    print("Found")
    gene_intervals <- makeGRangesFromDataFrame(read_tsv(paste(genomes_path,
                                                              gene_file_name,
                                                              sep = "/"),
                                                        col_names = TRUE),
                                               keep.extra.columns = TRUE)
  } else{
    if(!exists("mart")){
      mart <- get_mart(genome = genome)
      "new mart"
    }
    gene_intervals <- get_genes_from_biomart(mart = mart)
    gene_intervals <- filter_RNA_from_intervals(gene_intervals)
    write_data_to_TSV(data = data.frame(gene_intervals),
                      path = genomes_path,
                      filename = gene_file_name)
  }

  exon_file_name <- paste(genome,"_exon_intervals.tsv", sep = "")
  if(exon_file_name %in% files){
    exon_intervals <- makeGRangesFromDataFrame(df = read_tsv(file = paste(genomes_path,
                                                                          exon_file_name,
                                                                          sep = "/"),
                                                             col_names = TRUE),
                                               keep.extra.columns = TRUE)
  } else{
    if(!exists("mart")){
      mart <- get_mart(genome = genome)
    }
    exon_intervals <- get_exons_from_biomart(mart = mart)
    exon_intervals <- filter_by_metadata(target = exon_intervals, source = gene_intervals, column = "ensembl_gene_id")
    write_data_to_TSV(data = data.frame(exon_intervals),
                      path = genomes_path,
                      filename = exon_file_name)
  }

  chromosome_sizes_file <- paste(genome,"_chromosome_sizes.tsv", sep = "")
  if(chromosome_sizes_file %in% files){
    chromosome_sizes <- read_tsv(paste(genomes_path,
                                       chromosome_sizes_file,
                                       sep = "/"),
                                 col_names = TRUE)
  } else{
    if(!exists("mart")){
      mart <- get_mart(genome = genome)
    }
    chromosome_sizes <- load_chromosome_sizes(genome)
    write_data_to_TSV(data = chromosome_sizes,
                      path = genomes_path,
                      filename = chromosome_sizes_file)
  }

  return(list(chromosome_sizes=chromosome_sizes,
              gene_intervals=gene_intervals,
              exon_intervals=exon_intervals))
}

write_data_to_TSV <- function(data, path, filename){
  write_delim(
    data,
    path = paste(path,filename, sep="/"),
    delim = "\t",
    col_names = TRUE)
}
