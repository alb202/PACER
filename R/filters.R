assign_5prime_to_a_length <- function(alignments, aln_length){
  positive <- alignments[qwidth(alignments)==aln_length & strand(alignments)=="+"]
  negative <- alignments[qwidth(alignments)==aln_length & strand(alignments)=="-"]
  positive_results <- findOverlaps(query = alignments, subject = positive, type = "start", select = "first", ignore.strand=FALSE)
  negative_results <- findOverlaps(query = alignments, subject = negative, type = "end", select = "first", ignore.strand=FALSE)
  positive_results <- replace(x = positive_results, !is.na(positive_results), FALSE)
  positive_results <- replace(x = positive_results, is.na(positive_results), TRUE)
  negative_results <- replace(x = negative_results, !is.na(negative_results), FALSE)
  negative_results <- replace(x = negative_results, is.na(negative_results), TRUE)
  full_results <- positive_results & negative_results
  return(sort.GenomicRanges(c(alignments[full_results], positive, negative)))
}

assign_5prime_to_longer <- function(alignments, use_longer=c(TRUE, FALSE)){

  lengths <- sort(unique.Vector(qwidth(alignments)), decreasing = use_longer)
  results <- c(alignments[qwidth(alignments)==lengths[1]])
  for (i in 1:(length(lengths)-1)){
    #higher <- alignments[qwidth(alignments) %in% lengths[1:i]]
    lower <- alignments[qwidth(alignments) %in% lengths[(i+1)]]
    positive <- results[strand(results)=="+"]
    negative <- results[strand(results)=="-"]
    positive_results <- findOverlaps(query = lower, subject = positive , type = "start", select = "first", ignore.strand=FALSE)
    negative_results <- findOverlaps(query = lower, subject = negative, type = "end", select = "first", ignore.strand=FALSE)
    positive_results <- replace(x = positive_results, !is.na(positive_results), FALSE)
    positive_results <- replace(x = positive_results, is.na(positive_results), TRUE)
    negative_results <- replace(x = negative_results, !is.na(negative_results), FALSE)
    negative_results <- replace(x = negative_results, is.na(negative_results), TRUE)
    full_results <- positive_results & negative_results
    results <- c(results, lower[full_results])
  }
  return(sort.GenomicRanges(results))
}

filter_alignments_by_size_range <- function(alignments, minimum=10, maximum=30){
  results <- alignments[qwidth(alignments)>=minimum & qwidth(alignments)<=maximum]
  return(sort.GenomicRanges(results))
}


filter_BAM_tags <- function(gr){
  # Get the index for alignments with no mismatches
  no_mismatches_index <- mcols(gr)$NM==0

  # Get the index for alignments with up to 2 mismatches
  two_mismatches_index <- mcols(gr)$NM<=2

  # Get the index for alignments with no mismatches in the first 22 bases
  MD_split <- strsplit(mcols(gr)$MD, split = "[A-Z]")
  setA <- ifelse(as.numeric(mcols(gr)$NM)<=0, TRUE, FALSE)
  setB <- ifelse(strand(gr)=="+"&unlist(lapply(X = MD_split, FUN = function(x) as.numeric(x[[1]][1])>=22)), TRUE, FALSE)
  setC <- ifelse(strand(gr)=="-"&unlist(lapply(X = MD_split, FUN = function(x) as.numeric(x[[length(x)]][1])>=22)), TRUE, FALSE)
  no_mismatches_in_seed_index <- setA | setB | setC
  return(list(no_mm=no_mismatches_index, two_mm=two_mismatches_index, no_mm_seed=no_mismatches_in_seed_index))
}
#
# filter_MD_tag <- function(strand, NM, MD){
#   if(NM==0)
#     return(TRUE)
#   MD_split <- strsplit(x = MD, split = "[A-Z]")
#   # return(filter_MD_tags3(strand = strand, MD = MD_split))
#   if(strand=="-")
#     MD_split <- rev(MD_split[[1]])
#   if((as.numeric(MD_split[[1]][1])>=as.numeric(22)))
#     return(TRUE)
#   else
#     return(FALSE)
  # MD <- strsplit(x = mcols(gr)$MD, split = "[A-Z]")
  # ((strand(gr)=="+" & unlist(lapply(X = MD, FUN = function(x) as.numeric(x[[1]])))>=22) |
  #    (strand(gr)=="-" & unlist(lapply(X = rev(MD), FUN = function(x) as.numeric(x[[1]])))>=22) | )
  #else()
#   #  return(FALSE)
# }


filter_by_metadata <- function(target, source, column){

  matches <- mcols(target)[,column] %in% mcols(source)[,column]
  results <- target[matches]
  return(sort.GenomicRanges(results))
}

filter_reads_by_regions <- function(alignments, regions, type=c("both", "sense", "antisense"), invert=FALSE){
  if (type=="both") {
    results <- subsetByOverlaps(query = alignments, subject = regions, invert = invert, ignore.strand=TRUE)
  }
  if (type=="sense") {
    results <- subsetByOverlaps(query = alignments, subject = regions, invert = invert, ignore.strand=FALSE)
  }
  if (type=="antisense") {
    strand_info <- as.character(strand(regions))
    strand(regions) <- invert_vector(strand_info)
    results <- subsetByOverlaps(query = alignments, subject = regions, invert = invert, ignore.strand=FALSE)
  }
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

remove_overrepresented_sequences <- function(alignments, cutoff=0.001){
  counts <- rle(sort(as.character(mcols(alignments)$seq)))
  counts_df <- data.frame(values=counts$values, lengths=counts$lengths)
  overreppresented <- subset(counts_df, lengths>(cutoff*length(alignments)))
  '%nin%' <- Negate('%in%')
  results <- subset(alignments, seq %nin% overreppresented$values)
  return(sort.GenomicRanges(results))
}
