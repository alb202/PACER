
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
  # results <- subset(x = intervals, type!="snRNA" & type!="miRNA_gene" & biotype!="rRNA" & biotype!="tRNA")
  return(sort.GenomicRanges(results))
}

### Chromosome lengths
#chrom_sizes <- seqlengths(WBcel235_txdb)

#subset_gene_intervals <- function(intervals, )
#WBcel235_genes <- import(con = "genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.gff3")
#WBcel235_otherRNA <- subset(x = WBcel235_genes, biotype=="snoRNA" | biotype=="miRNA" | biotype=="rRNA" | biotype=="tRNA")
#WBcel235_genes_filtered <- subset(x = WBcel235_genes, biotype!="snoRNA" & biotype!="miRNA" & biotype!="rRNA" & biotype!="tRNA")
#WBcel235_genes_filtered <- subset(x = WBcel235_genes, biotype!="snoRNA" & biotype!="miRNA" & biotype!="rRNA" & biotype!="tRNA")


load_alignments <- function(path, params=ScanBamParam(reverseComplement = FALSE, what = c("seq", "qname", "flag"))){
  results <- readGAlignments(file = path, param = params)
  return(sort.GenomicRanges(results))
}





#two_mismatches <- subsetByOverlaps(query = two_mismatches, subject = WBcel235_otherRNA, invert = TRUE)
#two_mismatches <- readGAlignments(file = "output/WT_early_rep1_full/two_mm_SRR5023999.bam", param = ScanBamParam(reverseComplement = FALSE, what = c("seq", "qname", "flag")))






#library(biomaRt)
#mart <- useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")

#dplyr::filter(features_files, Name=="miRNA")[["WBcel235"]]
# dplyr::filter(features_files, Name=="miRNA")[[genome]]

#readGFF(filepath = paste(getwd(),"genomes", genome, dplyr::filter(features_files,
#                                                                  Name=="miRNA")[[genome]], sep="/"))
#bam <- readGAlignments(file = "output/WT_early_rep1/zero_mm_SRR5023999_100K_sample.bam")
#filter(.data = genome_files, type=="filter")$WBcel235


#new_bam <- subsetByOverlaps(subject = miRNA, query = bam, invert = TRUE)
#load_intervals <- function(filepath, filetype){
#  a <- readGAlignments(file = "output/WT_early_rep1/zero_mm_SRR5023999_100K_sample.bam")
#a <- scanBam(file = "output/WT_early_rep1/zero_mm_SRR5023999_100K_sample.bam",
#             param = ScanBamParam(what = scanBamWhat(), flag = scanBamFlag(isUnmappedQuery = FALSE)))
#a <- import(con = "output/WT_early_rep1/two_mm_SRR5023999_100K_sample.bam")
#}


# miRNA <- import("genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.miRNA.gff3")
# # snoRNA <- import("genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.snoRNA.gff3")
# rRNA <- import("genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.rRNA.gff3")
# tRNA <- import("genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.tRNA.gff3")










#two_mismatches_5prime_assigned_to_longer <- assign_5prime_to_longer(alignments = two_mismatches_filtered, use_longer = TRUE)
#write.table(x = two_mismatches_5prime_assigned_to_longer, file = "two_mismatches_5prime_assigned_to_longer.tsv", sep="\t", quote = FALSE)







# print(length(positive))
# print(length(negative))
# print(length(two_mismatches_sorted[full_results]))
# sum(full_results)
# write.table(x = negative, file = "negative.tsv", sep="\t", quote = FALSE)
# write.table(x = positive, file = "positive.tsv", sep="\t", quote = FALSE)


#two_mismatches_overlap_22 <- subset(two_mismatches_filtered, width == 22)
#new_counts <- rle(sort(as.character(mcols(two_mismatches_filtered)$seq)))
#qplot(x = new_counts$values, y = new_counts$lengths)

#all <- data.frame(seqnames(two_mismatches_filtered), strand(two_mismatches_filtered), start(two_mismatches_filtered), width(two_mismatches_filtered))
#plus <- subset(x = all,  width.two_mismatches_filtered.=="22" & strand.two_mismatches_filtered.=="+")
#minus <- subset(x = all,  width.two_mismatches_filtered.=="22" & strand.two_mismatches_filtered.=="-")
#two_mismatches[sample(c(TRUE, FALSE), size = length(two_mismatches_filtered), replace = TRUE, prob = c(.3,.7))]

# assign_5prime_to_a_length <- function(alignments, aln_length){
#   df_length <- length(alignments)
#   print(df_length)
#   alignments <- sort.GenomicRanges(alignments)
#   read_info <- data.frame(chr=seqnames(alignments), str=strand(alignments), start=start(alignments), end=end(alignments), len=width(alignments), stringsAsFactors = FALSE)
#   print(length(read_info$chr))
#   print(names(read_info))
#   reference_reads <- subset(x = read_info, len==aln_length)
#   print(names(reference_reads))
#   print(reference_reads)
#   print(length(reference_reads$chr))
#   results <- rep(TRUE, df_length)
#   for(i in 1:length(read_info[,1])){
#     #print(i)
#     for(j in 1:length(reference_reads[,1])){
#       #print(j)
#       if (read_info$chr[i] == reference_reads$chr[j]){
#         #print("chr match")
#         if (((read_info$str[i] == "+") & (read_info$str[i] == reference_reads$str[j]) & (read_info$start[i] == reference_reads$start[j]) & (reference_reads$len[j] != aln_length)) |
#             ((read_info$str[i] == "-") & (read_info$str[i] == reference_reads$str[j]) & (read_info$end[i] == reference_reads$end[j]) & (reference_reads$len[j] != aln_length))){
#               results[i] <- FALSE
#               print(FALSE)
#       }}}}
#   print(df_length)
#   return(results)
# }
#results <- assign_5prime_to_a_length(alignments = test, aln_length = 22)
#sum(results)

#continue
#results <- rep(TRUE, length(a[,1]))
# for (i in 1:length(a[,1])){
#   for (j in 1:length(b[,1])){
#     results[i] <- (a[i,]!=b[j,])[2] & results[i]
#   }
# }
#results
#minus <- subset(x = read_info, len=="22" & str=="-")

### Find by overlaps does NOT switch the start and end with reverse strand reads
# a <- two_mismatches[seqnames(two_mismatches)=="I" & strand(two_mismatches)=="-"]
# a
# a[319146:319147]
# b <- a[319146:319147]
# findOverlaps(query = b[1], subject = b[2], type = "start", ignore.strand=FALSE)
# findOverlaps(query = b[1], subject = b[2], type = "end", ignore.strand=FALSE)

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

#getChromInfoFromBiomart(biomart=mart, dataset=genome)
