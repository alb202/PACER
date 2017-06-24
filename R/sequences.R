get_sequence <- function(genome, gr){
  # Run the lookup on each interval, then add as metadata column
  mcols(gr)["DNAseq"] <- mcmapply(FUN = get_genome_sequence,
                                  as.character(seqnames(gr)),
                                  as.numeric(start(gr)-10),
                                  as.numeric(end(gr)+10),
                                  as.character(strand(gr)),
                                  MoreArgs = list(sequences = genome),
                                  USE.NAMES = FALSE)
  # Sort the data and return
  return(sort.GenomicRanges(gr))
}

rev_comp_DNA <- function(x){
  # Create the lookup list
  RC <- list(A="T", C="G", G="C", T="A")
  #Split the string into a vector
  split_string <- strsplit(x = x, split = "")[[1]]
  #Run lapply on the vector to swap values, then reverse it, then turn back into string
  result <- paste(rev(unname(unlist(mclapply(X = x, FUN = function(x) RC[x])))), collapse = '')
  return(result)
}

get_genome_sequence <- function(chromosome, start, end, strand, sequences){
  result <- as.character(sequences[[chromosome]][start:end])
  if (strand=="-")
    result <- rev_comp_DNA(result)
  return(result)
}


get_DNA_sequence <- function(genome, gr, start, end, type=c("5", "3", "full"), name){

  if (start< (-1) | end< (-1))
    stop("The start and end values must be integers -1 or greater. They indicate how many additional bases should be included on each end.")
  # Split the file by strand
  plus <- gr[strand(gr)=="+"]
  minus <- gr[strand(gr)=="-"]
  # Calculate regions to extract from genome
  if (type=="5"){
    plus_start <- start(plus)-start     #
    plus_end <- start(plus)+end         #
    minus_start <- end(minus)-end      #
    minus_end <- end(minus)+start      #
  } else if (type=="3"){
    plus_start <- end(plus)-start       #
    plus_end <- end(plus)+end           #
    minus_start <- start(minus)-end    #
    minus_end <- start(minus)+start        #
  } else if (type=="full"){
    plus_start <- start(plus)-start     #
    plus_end <- end(plus)+end           #
    minus_start <- start(minus)-end    #
    minus_end <- end(minus)+start          #
  }

  # Use mcMapply to get get each region
  if (length(plus) > 0){
    mcols(plus)[name] <- unlist(mclapply(FUN = as.character, X = mcmapply(FUN = get_genome_sequence,
                                                        as.character(seqnames(plus)),
                                                        as.numeric(plus_start),
                                                        as.numeric(plus_end),
                                                        as.character(strand(plus)),
                                                        MoreArgs = list(sequences = genome),
                                                        USE.NAMES = FALSE)))
  }
  if (length(minus) > 0){
    mcols(minus)[name] <- unlist(mclapply(FUN = as.character, X = mcmapply(FUN = get_genome_sequence,
                                                        as.character(seqnames(minus)),
                                                        as.numeric(minus_start),
                                                        as.numeric(minus_end),
                                                        as.character(strand(minus)),
                                                        MoreArgs = list(sequences = genome),
                                                        USE.NAMES = FALSE)))
  }
  # Recombine the plus and minus strand results
  results <- c(plus,minus)
  return(sort.GenomicRanges(results))
}

### Old version using DNAString - much slower
# get_sequence_from_genome <- function(chromosome, start, end, strand, sequences){
#   result <- DNAString(sequences[[chromosome]][start:end])
#   if (strand=="-")
#     result <- reverseComplement(result)
#   return(result)
# }

