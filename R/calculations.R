
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

