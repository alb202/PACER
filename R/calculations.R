
count_overlaps_by_width <- function(gr, regions, overlap = "sense", normalized=FALSE){
  widths <- sort(base::unique(width(gr)))
  results <- as.data.frame(matrix(ncol=length(widths), nrow=length(regions)))
  names(results) <- widths
  original_region_strand <-as.character(strand(regions))
  if (overlap=="sense")
    ignore_strand <- FALSE
  if (overlap=="both")
    ignore_strand <- TRUE
  if (overlap=="antisense"){
    ignore_strand <- FALSE
    strand(regions) <- invert_vector(as.character(strand(regions)))
  }
  if (normalized==TRUE)
    interval_widths <- (end(regions) - start(regions))
  for (i in 1:length(widths)){
    results[,i] <- GenomicRanges::countOverlaps(subject = gr[width(gr)==widths[i]],
                                            query = regions,
                                            type = "any",
                                            ignore.strand=ignore_strand)
    if (normalized==TRUE)
      results[,i] <- results[,i] / interval_widths
  }
  return(cbind.data.frame("Gene_strand"=original_region_strand, results))
}

count_overlaps_by_width_and_base <- function(gr, regions, alignment_width, base_col, overlap = "sense", normalized=FALSE){
  gr <- gr[width(gr)==alignment_width]
  bases <- sort(base::unique(mcols(gr)[[base_col]]))
  results <- as.data.frame(matrix(ncol=length(bases), nrow=length(regions)))
  names(results) <- bases
  original_region_strand <-as.character(strand(regions))
  if (overlap=="sense")
    ignore_strand <- FALSE
  if (overlap=="both")
    ignore_strand <- TRUE
  if (overlap=="antisense"){
    ignore_strand <- FALSE
    strand(regions) <- invert_vector(as.character(strand(regions)))
  }
  if (normalized==TRUE)
    interval_widths <- (end(regions) - start(regions))
  for (i in 1:length(bases)){
    results[,i] <- GenomicRanges::countOverlaps(subject = gr[mcols(gr)[[base_col]]==bases[i]],
                                                query = regions,
                                                type = "any",
                                                ignore.strand=ignore_strand)
    if (normalized==TRUE)
      results[,i] <- results[,i] / interval_widths
  }
  return(cbind.data.frame("Gene_strand"=original_region_strand, results))
}


calculate_offsets <- function(gr, primary_length, overlap_type="sense", maximum_offset=10){
  # If the type is 'sense', the 5' ends of reads from the same strands are compared
  if(overlap_type=="sense"){
    strands <- list(c("+", "+"), c("-", "-"))
    ignore_strand = FALSE
  }
  # If the type is 'antisense', the 5' end of the primary set is compared to the 3' end of the secondary set
  if(overlap_type=="antisense"){
    strands <- list(c("+", "-"), c("-", "+"))
    ignore_strand = TRUE
  }
  # Make the 4 sets of GenomicRanges
  set1_1 <- granges(gr[width(gr)==primary_length & strand(gr)==strands[[1]][1]])
  set1_2 <- granges(gr[width(gr)!=primary_length & strand(gr)==strands[[1]][2]])
  set2_1 <- granges(gr[width(gr)==primary_length & strand(gr)==strands[[2]][1]])
  set2_2 <- granges(gr[width(gr)!=primary_length & strand(gr)==strands[[2]][2]])

  # Find the overlaps that are within the maximum offset
  set1 <- GenomicRanges::findOverlaps(query = set1_1, subject = set1_2, maxgap = maximum_offset, type = "start", select = "all", ignore.strand = ignore_strand)
  set2 <- GenomicRanges::findOverlaps(query = set2_1, subject = set2_2, maxgap = maximum_offset, type = "end", select = "all", ignore.strand = ignore_strand)

  # Use the 'hits' lists to make list of the GenomicRanges
  pos_query_grs <- set1_1[queryHits(set1)]
  pos_subject_grs <- set1_2[subjectHits(set1)]
  neg_query_grs <- set2_1[queryHits(set2)]
  neg_subject_grs <- set2_2[subjectHits(set2)]

  # Create a data frame with the values from each strand, and bind them into a single data frame
  total_results <- rbind.data.frame(data.frame("offsets" = start(pos_subject_grs) - start(pos_query_grs),
                                               "widths"= width(pos_subject_grs),
                                               #"chromosomes"= seqnames(pos_subject_grs),
                                               "strands"= strand(pos_query_grs)),
                                    data.frame("offsets" = end(neg_query_grs) - end(neg_subject_grs),
                                               "widths"= width(neg_subject_grs),
                                               #"chromosomes"= seqnames(neg_subject_grs),
                                               "strands"= strand(neg_query_grs)))

  # Count the unique rows and return a data table
  # res_1 <- data.table(dplyr::count(x = total_results, offsets, widths, chromosomes, strands, sort = TRUE))
  # res_2 <- data.table(dplyr::count(x = total_results, widths, chromosomes, strands, sort = TRUE))
  # results <- dplyr::left_join(x = res_1, y = res_2, by = c("widths", "chromosomes", "strands"), copy = TRUE, suffix = c(".x", ".y"))

  res_1 <- data.table(dplyr::count(x = total_results, offsets, widths, strands, sort = TRUE))
  res_2 <- data.table(dplyr::count(x = total_results, widths, strands, sort = TRUE))
  results <- dplyr::left_join(x = res_1, y = res_2, by = c("widths", "strands"), copy = TRUE, suffix = c(".x", ".y"))


  results["ratio"] <- results$n.x / results$n.y
  return(results)
}


find_overlaps <- function(A, B){
  a <- subsetByOverlaps(query = A, subject = B, maxgap = 10, minoverlap = 1, type = "start", invert = FALSE)
  return(start(A) - start(a))
}

find_minimum <- function(A, B){
  return(B[which.min(abs(B-A))]-A)
}

gr <- sort.GenomicRanges(two_mm[sample(x = 1:3499785, size = 1000000, replace = FALSE)])

#p <- ggplot(data = a, mapping = aes(x = offsets, by = as.factor(a$widths))) + geom_freqpoly(binwidth = 1)
#p <- ggplot2::ggplot() + geom_line(mapping = aes(x = results$offsets, y = results$n, group = results$widths, color = results$widths))

ggplot(data=results3, aes(x=offsets, y=ratio)) + geom_line(aes(group = widths), color=results3$widths, inherit.aes = TRUE, show.legend = TRUE) + facet_grid(chromosomes~strands) + ggplot2::theme(legend.position = "right")

adjust_ggplot_units <- function(maxWidths, x1, x2){
  widths <- as.numeric(maxWidth2[[2]][[1]][[2]])
  y1 <- maxWidth[[2]][[1]][[1]]*(x1/x2)
  y2 <- maxWidth[[2]][[2]][[1]]
  new_units <- unit.pmax(unit(x = list(y1), units = list("cm")),
                         unit(x = list(y2), units = list("cm")))
  return(new_units)
}
