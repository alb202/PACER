
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

calculate_offsets <- function(gr, primary_length, overlap_type="sense"){
  if(overlap_type=="sense")
    strands <- list(c("+", "+"), c("-", "-"))
  if(overlap_type=="antisense")
    strands <- list(c("+", "-"), c("-", "+"))
  gr <- gr[seqnames(gr)=="X"]
  set1_1 <- granges(gr[width(gr)==primary_length & strand(gr)==strands[[1]][1]])
  set1_2 <- granges(gr[width(gr)!=primary_length & strand(gr)==strands[[1]][2]])
  set2_1 <- granges(gr[width(gr)==primary_length & strand(gr)==strands[[2]][1]])
  set2_2 <- granges(gr[width(gr)!=primary_length & strand(gr)==strands[[2]][2]])

  z <- mapply(FUN = CalculateOffset,
              as.character(seqnames(set1_1)),
              as.numeric(start(set1_1)),
              MoreArgs = list(as.character(seqnames(set1_2)),
                              as.numeric(start(set1_2)),
                              as.numeric(width(set1_2)),
                              10),
              USE.NAMES = FALSE,
              SIMPLIFY = TRUE)
  results <- data.frame("chromosomes"=unlist(z[1, ]), "offsets"=unlist(z[2, ]), "widths"=unlist(z[3, ]))
  return(results)
}
              # Rcpp::List CalculateOffset(std::string primaryChromosome,
              #                            int primaryPosition,
              #                            Rcpp::CharacterVector secondaryChromosome,
              #                            Rcpp::IntegerVector secondaryPosition,
              #                            Rcpp::IntegerVector secondaryWidth,
              #                            const int& maxOffset){


  ### This is all wrong. Go interval by interval with find overlap and then for each interval subtract the start positions to get an integer for each one.

### for each comparison, send the two sets to C++. then compare and return a list containing a length vector and an offset vector, and maybe a Chromasome vector

  # A <- c(1,4,6,8,25,58,87)
  # B <- c(8,9,23,34,84,6,128)
  #calculate_distance(A = A, B = B)
  #  system.time(X <-  mapply(FUN = find_minimum, start(set1_1), MoreArgs = list(start(set1_2)), USE.NAMES = FALSE))

  #  unlist(mapply(FUN = find_overlaps, as.list(set1_1), MoreArgs = list(set1_2), USE.NAMES = FALSE, SIMPLIFY = TRUE))
  #set1 <- GenomicRanges::findOverlaps(query = set1_1, subject = set1_2, maxgap = 10, type = "any")
  #set2 <- GenomicRanges::findOverlaps(query = set2_1, subject = set2_2, maxgap = 10, type = "any")
  # set1_offsets <- CalculateOffset(x = as.matrix(set1),A = start(set1_1), B = start(set1_2))
  # set2_offsets <- CalculateOffset(x = as.matrix(set2),A = end(set2_1), B = end(set2_2))
  # set1_secondary_widths <- width(set1_2[subjectHits(set1)])
  # set1_secondary_strands <- strand(set1_2[subjectHits(set1)])
  # set1_secondary_chromosome <- seqnames(set1_2[subjectHits(set1)])
  # set2_secondary_widths <- width(set2_2[subjectHits(set2)])
  # set2_secondary_strands <- strand(set2_2[subjectHits(set2)])
  # set2_secondary_chromosome <- seqnames(set2_2[subjectHits(set2)])
  # set1_df <- data.frame(stringsAsFactors = FALSE,
  #                       set1_offsets,
  #                       width(set1_2[subjectHits(set1)]),
  #                       strand(set1_2[subjectHits(set1)]),
  #                       seqnames(set1_2[subjectHits(set1)]))
  # set2_df <- data.frame(stringsAsFactors = FALSE,
  #                       set2_offsets,
  #                       width(set2_2[subjectHits(set2)]),
  #                       strand(set2_2[subjectHits(set2)]),
  #                       seqnames(set2_2[subjectHits(set2)]))
  # names(set1_df) <- c("offsets", "widths", "strands", "seqname")
  # names(set2_df) <- c("offsets", "widths", "strands", "seqname")
  # return(rbind(set1_df, set2_df))
  # }

find_overlaps <- function(A, B){
  a <- subsetByOverlaps(query = A, subject = B, maxgap = 10, minoverlap = 1, type = "start", invert = FALSE)
  return(start(A) - start(a))
}

find_minimum <- function(A, B){
  return(B[which.min(abs(B-A))]-A)
}


DFTest(primaryChromosome = c("I", "II", "III"), primaryStrand = c("+", "-", "+"), primaryPosition = c(100, 200, 300), secondaryChromosome = c("IV", "V", "X"), secondaryStrand = c("-", "+", "-"), primaryPosition = c(400, 500, 600))
a <- data.frame(mapply(FUN = CalculateOffset, list("I", "II", "III"), list(1000, 1002, 1004), MoreArgs = list(secondaryChromosome = c("I", "IV", "II", "I", "X", "III", 'IV', "I"), secondaryPosition = c(1001, 1002, 1003, 1004, 1005, 997, 998, 999), secondaryWidth = c(18, 19, 20, 21, 22, 25, 28, 18), maxOffset = 10)))
