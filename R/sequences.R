get_DNA_sequence <- function(genome, gr, start, end, type=c("5", "3", "full"), name){

  if (start< (-1) | end< (-1))
    stop("The start and end values must be integers -1 or greater. They indicate how many additional bases should be included on each end.")

  plus <- gr[strand(gr)=="+"]
  minus <- gr[strand(gr)=="-"]
  # process each
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

  # mapply each
  if (length(plus) > 0){
    #mcols(plus)[name] <- mapply(FUN = get_sequence_from_genome, as.character(seqnames(plus)), as.numeric(plus_start), as.numeric(plus_end), as.character(strand(plus)), MoreArgs = list(sequences = genome))
    mcols(plus)[name] <- unlist(unname(lapply(FUN = as.character, X = mcmapply(FUN = get_sequence_from_genome, as.character(seqnames(plus)),
                                           as.numeric(plus_start), as.numeric(plus_end), as.character(strand(plus)), MoreArgs = list(sequences = genome)))))
  }
  if (length(minus) > 0){
    #mcols(minus)[name] <- mapply(FUN = get_sequence_from_genome, as.character(seqnames(minus)), as.numeric(minus_start), as.numeric(minus_end), as.character(strand(minus)), MoreArgs = list(sequences = genome))
    mcols(minus)[name] <- unlist(unname(lapply(FUN = as.character, X = mcmapply(FUN = get_sequence_from_genome, as.character(seqnames(minus)),
                                           as.numeric(minus_start), as.numeric(minus_end), as.character(strand(minus)), MoreArgs = list(sequences = genome)))))
  }
  # combine the results
  results <- c(plus,minus)
  return(sort.GenomicRanges(results))
  }

get_sequence_from_genome <- function(chromosome, start, end, strand, sequences){
  result <- DNAString(sequences[[chromosome]][start:end])
  if (strand=="-")
    result <- reverseComplement(result)
  return(result)
}
