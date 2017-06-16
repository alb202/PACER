library(GenomicAlignments)
#library(VariantAnnotation)
library(rtracklayer)
library(dplyr)
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


WBcel235_genes <- import(con = "genomes/WBcel235/Caenorhabditis_elegans.WBcel235.89.gff3")
WBcel235_otherRNA <- subset(x = WBcel235_genes, biotype=="snoRNA" | biotype=="miRNA" | biotype=="rRNA" | biotype=="tRNA")
WBcel235_genes_filtered <- subset(x = WBcel235_genes, biotype!="snoRNA" & biotype!="miRNA" & biotype!="rRNA" & biotype!="tRNA")
two_mismatches <- readGAlignments(file = "output/WT_early_rep1_full/two_mm_SRR5023999.bam", param = ScanBamParam(reverseComplement = FALSE, what = c("seq", "qname", "flag")))
two_mismatches <- subsetByOverlaps(query = two_mismatches, subject = WBcel235_otherRNA, invert = TRUE)
two_mismatches <- two_mismatches[qwidth(two_mismatches)>=length_range["minimum"] & qwidth(two_mismatches)<=length_range["maximum"]]
counts <- rle(sort(as.character(mcols(two_mismatches)$seq)))
counts_df <- data.frame(values=counts$values, lengths=counts$lengths)
overreppresented <- subset(counts_df, lengths>(0.001*length(two_mismatches)))
'%nin%' <- Negate('%in%')
two_mismatches_filtered <- subset(two_mismatches, seq %nin% overreppresented$values)








two_mismatches_5prime_assigned_to_longer <- assign_5prime_to_longer(alignments = two_mismatches_filtered, use_longer = TRUE)
write.table(x = two_mismatches_5prime_assigned_to_longer, file = "two_mismatches_5prime_assigned_to_longer.tsv", sep="\t", quote = FALSE)







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

