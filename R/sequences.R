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
  # print(plus)
  # print(data.frame(seqnames(plus), plus_start, plus_end, strand(plus), stringsAsFactors = FALSE))
  # print(length(as.character(seqnames(plus))))
  # print(length(as.numeric(plus_start)))
  # print(length(as.numeric(plus_end)))
  # print(length(as.character(strand(plus))))



  # mapply each
  if (length(plus) > 0){
    #mcols(plus)[name] <- mapply(FUN = get_sequence_from_genome, as.character(seqnames(plus)), as.numeric(plus_start), as.numeric(plus_end), as.character(strand(plus)), MoreArgs = list(sequences = genome))
    mcols(plus)[name] <- unlist(unname(lapply(FUN = as.character, X = mapply(FUN = get_sequence_from_genome, as.character(seqnames(plus)),
                                           as.numeric(plus_start), as.numeric(plus_end), as.character(strand(plus)), MoreArgs = list(sequences = genome)))))
  }
  if (length(minus) > 0){
    #mcols(minus)[name] <- mapply(FUN = get_sequence_from_genome, as.character(seqnames(minus)), as.numeric(minus_start), as.numeric(minus_end), as.character(strand(minus)), MoreArgs = list(sequences = genome))
    mcols(minus)[name] <- unlist(unname(lapply(FUN = as.character, X = mapply(FUN = get_sequence_from_genome, as.character(seqnames(minus)),
                                           as.numeric(minus_start), as.numeric(minus_end), as.character(strand(minus)), MoreArgs = list(sequences = genome)))))
  }
  # print(class(a))
  # print(a)
  # print(class(b))
  # print(b)
  # combine the results
  results <- c(plus,minus)
  #return(sort(sortSeqlevels(results)))
  return(sort.GenomicRanges(results))
  }

get_sequence_from_genome <- function(chromosome, start, end, strand, sequences){
  # print(paste(chromosome, start, end, strand))
  # print(sequences)
  result <- DNAString(sequences[[chromosome]][start:end])
  if (strand=="-")
    result <- reverseComplement(result)
  # print(class(result))
  # print(class(result[1]))
  # print(result)
  # print(result[1])
  # print(class(result[[1]]))
  return(result)
}

# One function for the mapply function
# split the plus and minus strand reads and operate on them separately

