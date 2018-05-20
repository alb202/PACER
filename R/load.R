# get_exons_from_biomart <- function(mart){
#   exons <- biomaRt::getBM(
#     attributes = c("strand", "exon_chrom_start", "exon_chrom_end", "chromosome_name",
#                    "external_gene_name", "external_gene_source", "ensembl_gene_id",
#                    "ensembl_exon_id", "rank"), mart = mart)
#   exons$strand <- swap_values(x = exons$strand, old = c(1, -1), new = c("+", "-"))
#   new_grange <- makeGRangesFromDataFrame(
#     df=exons, seqnames.field = "chromosome_name", start.field = "exon_chrom_start",
#     end.field = "exon_chrom_end", strand.field = "strand", keep.extra.columns = TRUE)
#   return(sort.GenomicRanges(new_grange))
# }

# get_genes_from_biomart <- function(mart){
#   genes <- biomaRt::getBM(
#     attributes = c("chromosome_name", "start_position", "end_position", "strand",
#                    "description", "external_gene_name", "external_gene_source",
#                    "ensembl_gene_id", "source", "status", "gene_biotype"), mart = mart)
#   genes$strand <- swap_values(x = genes$strand, old = c(1, -1), new = c("+", "-"))
#   new_grange <- makeGRangesFromDataFrame(
#     df=genes, seqnames.field = "chromosome_name", start.field = "start_position",
#     end.field = "end_position", strand.field = "strand", keep.extra.columns = TRUE)
#   return(sort.GenomicRanges(new_grange))
# }


check_for_intervals <- function(path, genome){
  genome_path <- getAbsolutePath(paste(path, genome, sep = '/'))
  if(dir.exists(genome_path)){
    files <- dir(path = genome_path)
    if((paste(genome, "_chr.tsv", sep = "") %in% files) &
       (paste(genome, "_genes.tsv", sep = "") %in% files) &
       (paste(genome, "_exons.tsv", sep = "") %in% files)){
      return(TRUE)
    }
  }
  return(FALSE)
}

get_ensembl_intervals <- function(path, genome, updateProgress = NULL){
  genome_path <- getAbsolutePath(paste(path, genome, sep = '/'))
  # If the directory doesn't exist, make it
  if(!dir.exists(genome_path))
    dir.create(genome_path)

  # Create the Mart and get dataset ID
  if (is.function(updateProgress)) {updateProgress(detail = "Creating BioMart", value = .0)}
  mart_name <- listMarts()[[1,1]]
  mart <- useMart(biomart = mart_name)
  datasets_list <- listDatasets(mart = mart)
  dataset <- datasets_list[datasets_list["version"]==genome][1]
  mart <- useDataset(dataset = dataset, mart = mart)


  # Get chromosome sizes
  if (is.function(updateProgress)) {updateProgress(detail = "Getting chromosome lengths", value = .2)}
  chromosome_sizes <- getChromInfoFromBiomart(biomart = mart_name, dataset = dataset)
  # Get gene intervals and filter non-siRNA
  if (is.function(updateProgress)) {updateProgress(detail = "Getting gene info", value = .4)}
  gene_intervals <- get_genes_from_biomart(mart = mart)
  gene_intervals <- filter_RNA_from_intervals(gene_intervals)
  # Get exon intervals and filter intervals
  if (is.function(updateProgress)) {updateProgress(detail = "Getting exon info", value = .6)}
  exon_intervals <- get_exons_from_biomart(mart = mart)
  exon_intervals <- filter_by_metadata(target = exon_intervals, source = gene_intervals, column = "ensembl_gene_id")
  # Calculate size of intervals
  if (is.function(updateProgress)) {updateProgress(detail = "Processing intervals", value = .8)}
  gene_intervals <- calculate_exonic_bp_of_gene(exons = exon_intervals, genes = gene_intervals)
  mcols(exon_intervals) <- data.frame(mcols(exon_intervals),
                                      exon_bp=as.integer(end(exon_intervals)-start(exon_intervals)+1),
                                      stringsAsFactors = FALSE)
  # Write genome lengths to file
  if (is.function(updateProgress)) {updateProgress(detail = "Saving files", value = .9)}
  write_data_to_TSV(data = chromosome_sizes,
                    path = genome_path,
                    filename = paste(genome, "_chr.tsv", sep = ""))
  # Write gene intervals to file
  write_data_to_TSV(data = data.frame(gene_intervals),
                    path = genome_path,
                    filename = paste(genome, "_genes.tsv", sep = ""))
  # Write exon intervals to file
  write_data_to_TSV(data = data.frame(exon_intervals),
                    path = genome_path,
                    filename = paste(genome, "_exons.tsv", sep = ""))
  #if (is.function(updateProgress)) {updateProgress(detail = "Done", value = 1)}
}

get_datasets_from_biomart <- function(updateProgress = NULL){
  if (is.function(updateProgress)) {updateProgress(detail = "Getting biomarts", value = .25)}
  mart_info <- listMarts()[1,]
  if (is.function(updateProgress)) {updateProgress(detail = "Selecting biomart", value = .5)}
  biomart <- useMart(biomart = mart_info[[1,1]])
  if (is.function(updateProgress)) {updateProgress(detail = "Getting datasets", value = .75)}
  ensembl_genome_index <- listDatasets(mart = biomart)
  return(ensembl_genome_index)
}

get_exons_from_biomart <- function(mart){
  exons <- biomaRt::getBM( mart = mart,
                           attributes = c("strand",
                                          "exon_chrom_start",
                                          "exon_chrom_end",
                                          "chromosome_name",
                                          "external_gene_name",
                                          "external_gene_source",
                                          "ensembl_gene_id",
                                          "ensembl_exon_id",
                                          "rank"))
  exons$strand <- swap_values(x = exons$strand,
                              old = c(1, -1),
                              new = c("+", "-"))
  return(sort.GenomicRanges(makeGRangesFromDataFrame(df=exons,
                                                     seqnames.field = "chromosome_name",
                                                     start.field = "exon_chrom_start",
                                                     end.field = "exon_chrom_end",
                                                     strand.field = "strand",
                                                     keep.extra.columns = TRUE)))
}

get_genes_from_biomart <- function(mart){
  genes <- biomaRt::getBM(mart = mart,
                          attributes = c("chromosome_name",
                                         "start_position",
                                         "end_position",
                                         "strand",
                                         "description",
                                         "external_gene_name",
                                         "external_gene_source",
                                         "ensembl_gene_id",
                                         "source", # "status",
                                         "gene_biotype"))
  genes$strand <- swap_values(x = genes$strand,
                              old = c(1, -1),
                              new = c("+", "-"))
  return(sort.GenomicRanges(makeGRangesFromDataFrame(df=genes,
                                                     seqnames.field = "chromosome_name",
                                                     start.field = "start_position",
                                                     end.field = "end_position",
                                                     strand.field = "strand",
                                                     keep.extra.columns = TRUE)))
}


load_genome_data <- function(path, genome){
  genome_path <- getAbsolutePath(paste(path, genome, sep = '/'))
  if(dir.exists(genome_path)){
    files <- dir(path = genome_path)
    if(!(paste(genome, "_chr.tsv", sep = "") %in% files &
         paste(genome, "_genes.tsv", sep = "") %in% files &
         paste(genome, "_exons.tsv", sep = "") %in% files)){
      return(NULL)
    }
  } else {
    return(NULL)
  }

  chromosome_sizes <- read_tsv(file = paste(genome_path,
                                            "/",
                                            genome,
                                            "_chr.tsv",
                                            sep = ""),
                               col_names = TRUE)
  exons <- makeGRangesFromDataFrame(df = read_tsv(file = paste(genome_path,
                                                               "/",
                                                               genome,
                                                               "_exons.tsv",
                                                               sep = ""),
                                                  col_names = TRUE),
                                    keep.extra.columns = TRUE)
  genes <- makeGRangesFromDataFrame(df = read_tsv(file = paste(genome_path,
                                                               "/",
                                                               genome,
                                                               "_genes.tsv",
                                                               sep = ""),
                                                  col_names = TRUE),
                                    keep.extra.columns = TRUE)
  return(list("chromosome_sizes"=chromosome_sizes,
              "gene_intervals"=genes,
              "exon_intervals"=exons))
}


# get_mart <- function(genome){
#   mart_name <- c(listMarts()[[1,1]])
#   mart <- useMart(biomart = mart_name)
#   datasets_list <- listDatasets(mart = mart)
#   dataset <- datasets_list[datasets_list["version"]==genome][1]
#   mart <- useDataset(dataset = dataset, mart = mart)
#   return(mart)
# }

load_alignments <- function(path,
                            params=ScanBamParam(reverseComplement = TRUE,
                                                what = c("seq", "qname", "flag"),
                                                tag = c("XA", "MD", "NM"))){
  return(sort.GenomicRanges(readGAlignments(file = path, param = params)))
}

# load_chromosome_sizes <- function(genome){
#   mart_name <- c(listMarts()[[1,1]])
#   mart_datasets <- listDatasets(useMart(biomart=mart_name, host="www.ensembl.org"))
#   dataset_id <- mart_datasets[mart_datasets$version==genome,][[1]]
#   #chromosomes <- getChromInfoFromBiomart(biomart=mart_name,
#   #                                     dataset=dataset_id)
#   return(getChromInfoFromBiomart(biomart=mart_name, dataset=dataset_id))
# }

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

load_gene_sets <- function(gene_sets){
  gene_sets <- strsplit(x = gene_sets, split = ";", fixed = FALSE)[[1]]
  #genome_path <- getAbsolutePath(paste(path, genome, sep = '/'))
  # Create an empty list
  results <- list()
  results["Full"] <- ""
  # Loop through the filenames, load the list of genes, and add as a named list
  for(i in 1:length(gene_sets)){
    filename <- get_filename(path = gene_sets[[i]], extension = FALSE)
    # split_path <- strsplit(x = , split = "/", fixed = TRUE)
    results[filename] <- list(scan(file = gene_sets[[i]], what = character(), quote = ""))
    #results[names(gene_sets)[i]] <- list(scan(paste(path, gene_sets[i], sep = "/")), character(), quote = "")
  }
  return(results)
}
#
# load_genome_data <- function(path, genome){
#   genomes_path <- paste(path, genome, sep="/")
#   files <- dir(path = genomes_path)
#
#   gene_file_name <- paste(genome,"_gene_intervals.tsv", sep = "")
#   if(gene_file_name %in% files){
#     print("Found")
#     gene_intervals <- makeGRangesFromDataFrame(read_tsv(paste(genomes_path,
#                                                               gene_file_name,
#                                                               sep = "/"),
#                                                         col_names = TRUE),
#                                                keep.extra.columns = TRUE)
#   } else{
#     if(!exists("mart")){
#       mart <- get_mart(genome = genome)
#       "new mart"
#     }
#     gene_intervals <- get_genes_from_biomart(mart = mart)
#     gene_intervals <- filter_RNA_from_intervals(gene_intervals)
#     write_data_to_TSV(data = data.frame(gene_intervals),
#                       path = genomes_path,
#                       filename = gene_file_name)
#   }
#   exon_file_name <- paste(genome,"_exons.tsv", sep = "")
#   exon_file_name <- paste(genome,"_exon_intervals.tsv", sep = "")
#   if(exon_file_name %in% files){
#     aexon_intervals <- makeGRangesFromDataFrame(df = read_tsv(file = paste(genomes_path,
#                                                                           exon_file_name,
#                                                                           sep = "/"),
#                                                              col_names = TRUE),
#                                                keep.extra.columns = TRUE)
#   } else{
#     if(!exists("mart")){
#       mart <- get_mart(genome = genome)
#     }
#     exon_intervals <- get_exons_from_biomart(mart = mart)
#     exon_intervals <- filter_by_metadata(target = exon_intervals, source = gene_intervals, column = "ensembl_gene_id")
#     write_data_to_TSV(data = data.frame(exon_intervals),
#                       path = genomes_path,
#                       filename = exon_file_name)
#   }
#
#   chromosome_sizes_file <- paste(genome,"_chromosome_sizes.tsv", sep = "")
#   if(chromosome_sizes_file %in% files){
#     chromosome_sizes <- read_tsv(paste(genomes_path,
#                                        chromosome_sizes_file,
#                                        sep = "/"),
#                                  col_names = TRUE)
#   } else{
#     if(!exists("mart")){
#       mart <- get_mart(genome = genome)
#     }
#     chromosome_sizes <- load_chromosome_sizes(genome)
#     write_data_to_TSV(data = chromosome_sizes,
#                       path = genomes_path,
#                       filename = chromosome_sizes_file)
#   }
#   # Get the total exonic bp length for each gene
#   gene_intervals <- calculate_exonic_bp_of_gene(exons = exon_intervals,
#                                              genes = gene_intervals)
#   # Get the length of each exon
#   mcols(exon_intervals) <- data.frame(mcols(exon_intervals),
#                                       exon_bp=end(exon_intervals)-start(exon_intervals)+1,
#                                       stringsAsFactors = FALSE)
#   return(list(chromosome_sizes=chromosome_sizes,
#               gene_intervals=gene_intervals,
#               exon_intervals=exon_intervals))
# }

write_data_to_TSV <- function(data, path, filename){
  write_delim(
    data,
    path = paste(path,filename, sep="/"),
    delim = "\t",
    col_names = TRUE)
}
