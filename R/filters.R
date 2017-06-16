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
  results <- c(alignments[full_results], positive, negative)
  # print(length(positive))
  # print(length(negative))
  # print(length(two_mismatches_sorted[full_results]))
  # sum(full_results)
  # write.table(x = negative, file = "negative.tsv", sep="\t", quote = FALSE)
  # write.table(x = positive, file = "positive.tsv", sep="\t", quote = FALSE)
  return(sort.GenomicRanges(results))
}
two_mismatches_5prime_filtered <- assign_5prime_to_a_length(two_mismatches_filtered, 22)



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
    # print(c("results", length(results)))
    # print(c("lower", length(lower)))
    # print(c("full_results", length(full_results)))
    # print(c("positive_results", length(positive_results)))
    # print(c("negative_results", length(negative_results)))
    results <- c(results, lower[full_results])
  }
  return(sort.GenomicRanges(results))
}
