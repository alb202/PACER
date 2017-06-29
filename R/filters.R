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
  no_mismatches_index <- mcols(gr)$NM==0
  two_mismatches_index <- mcols(gr)$NM<=2
  MD <- strsplit(x = x, split = "[A-Z]")
  no_mismatches_in_seed_index <- (strand(gr)=="+" & unlist(lapply(X = MD, FUN = function(x) x[1]))>=22) |
    (strand(gr)=="-" & unlist(lapply(X = MD, FUN = function(x) x[length(x)]))>=22)
  return(list(no_mm=no_mismatches_index, two_mm=two_mismatches_index, no_mm_seed=no_mismatches_in_seed_index))
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
