get_genome_sequence <- function(gr, genome_sequence){
  mcols(gr) <- cbind.data.frame(
    mcols(gr),
    data.frame(
      get_sequence(
        chrom = as.character(seqnames(gr)),
        start = as.numeric(start(gr)),
        end = as.numeric(end(gr)),
        strand = as.character(strand(gr)),
        genome = genome_sequence),
      stringsAsFactors = FALSE))
  return(gr)
}

# get_sequence <- function(genome, gr){
#   plus <- gr[strand(gr)=="+"]
#   minus <- gr[strand(gr)=="-"]
#   # Run the lookup on each interval, then add as metadata column
#   mcols(plus)["DNAseq"] <- mcmapply(FUN = get_genome_sequence_plus,
#                                     as.character(seqnames(plus)),
#                                     as.numeric(start(plus)-10),
#                                     as.numeric(end(plus)+10),
#                                     #as.character(strand(plus)),
#                                     MoreArgs = list(sequences = genome),
#                                     USE.NAMES = FALSE)
#
#   mcols(minus)["DNAseq"] <- unlist(mclapply(mcmapply(FUN = get_genome_sequence_plus,
#                                                      as.character(seqnames(minus)),
#                                                      as.numeric(start(minus)-10),
#                                                      as.numeric(end(minus)+10),
#                                                      #as.character(strand(minus)),
#                                                      MoreArgs = list(sequences = genome),
#                                                      USE.NAMES = FALSE), FUN = function(x) reverse_complement(x)))
#                                                     #USE.NAMES = FALSE), FUN = function(x) as.character(reverseComplement(DNAString(x)))))
#   results <- c(plus,minus)
#
#   # Sort the data and return
#   #results <-
#   #return(sort.GenomicRanges(c(plus,minus)))
#   #results <- lapply(X = mcols(results[mcols(z)$DNAseq!=FALSE])$DNAseq, FUN = sequence_slice)
#   return(sort.GenomicRanges(results[mcols(results)$DNAseq!=FALSE]))
# }

split_sequences <- function(gr, column="DNAseq"){
  tmp <- mcols(gr)
  mcols(gr) <- cbind(tmp,
                     rbindlist(l = ,
                               lapply(X = tmp[[column]], FUN = sequence_slice),
                               use.names = TRUE),
                     deparse.level = 0)
  return(gr)
}
#
# get_genome_sequence <- function(chromosome, start, end, strand, sequences){
#   #result <- DNAString("AAA")
#   result <- sequences[[chromosome]][start:end]
#   if (strand=="-")
#     result <- reverseComplement(result)
#   return(as.character(result))
# }
#
# get_genome_sequence_plus <- function(chromosome, start, end, sequences){
#   return(as.character(sequences[[chromosome]][start:end]))
# }
#
# get_genome_sequence_minus <- function(chromosome, start, end, sequences){
#   return(as.character(reverseComplement(sequences[[chromosome]][start:end])))
# }
#
# get_genome_sequence_plus2 <- function(chromosome, start, end, sequences){
#   if (start<1) return(FALSE)
#   if (end>length(sequences[[chromosome]])) return(FALSE)
#   return(paste(sequences[[chromosome]][start:end], collapse = ''))
# }

# get_genome_sequence_plus2_cmp <- cmpfun(get_genome_sequence_plus2)
#
# get_genome_sequence_minus2 <- function(chromosome, start, end, sequences){
#   return(as.character(reverseComplement(DNAString(paste(sequences[[chromosome]][start:end], collapse = '')))))
# }


#
# rev_comp_DNA <- function(old_seq){
#   # Create the lookup list
#   RC <- list(A="T", C="G", G="C", T="A")
#   #Split the string into a vector
#   #split_seq <- strsplit(x = old_seq, split = "")[[1]]
#   #Run lapply on the vector to swap values, then reverse it, then turn back into string
#   #return(paste(rev(unname(unlist(mclapply(X = old_seq, FUN = function(x) RC[x])))), collapse = ''))
#   return(rev(unname(unlist(mclapply(X = old_seq, FUN = function(x) RC[x])))))
#   #return(result)
# }

# get_DNA_sequence <- function(genome, gr, start, end, type=c("5", "3", "full"), name){
#
#   if (start< (-1) | end< (-1))
#     stop("The start and end values must be integers -1 or greater. They indicate how many additional bases should be included on each end.")
#   # Split the file by strand
#   plus <- gr[strand(gr)=="+"]
#   minus <- gr[strand(gr)=="-"]
#   # Calculate regions to extract from genome
#   if (type=="5"){
#     plus_start <- start(plus)-start     #
#     plus_end <- start(plus)+end         #
#     minus_start <- end(minus)-end      #
#     minus_end <- end(minus)+start      #
#   } else if (type=="3"){
#     plus_start <- end(plus)-start       #
#     plus_end <- end(plus)+end           #
#     minus_start <- start(minus)-end    #
#     minus_end <- start(minus)+start        #
#   } else if (type=="full"){
#     plus_start <- start(plus)-start     #
#     plus_end <- end(plus)+end           #
#     minus_start <- start(minus)-end    #
#     minus_end <- end(minus)+start          #
#   }
#
#   # Use mcMapply to get get each region
#   if (length(plus) > 0){
#     mcols(plus)[name] <- unlist(mclapply(FUN = as.character, X = mcmapply(FUN = get_genome_sequence,
#                                                         as.character(seqnames(plus)),
#                                                         as.numeric(plus_start),
#                                                         as.numeric(plus_end),
#                                                         as.character(strand(plus)),
#                                                         MoreArgs = list(sequences = genome),
#                                                         USE.NAMES = FALSE)))
#   }
#   if (length(minus) > 0){
#     mcols(minus)[name] <- unlist(mclapply(FUN = as.character, X = mcmapply(FUN = get_genome_sequence,
#                                                         as.character(seqnames(minus)),
#                                                         as.numeric(minus_start),
#                                                         as.numeric(minus_end),
#                                                         as.character(strand(minus)),
#                                                         MoreArgs = list(sequences = genome),
#                                                         USE.NAMES = FALSE)))
#   }
#   # Recombine the plus and minus strand results
#   results <- c(plus,minus)
#   return(sort.GenomicRanges(results))
# }

### Old version using DNAString - much slower
# get_sequence_from_genome <- function(chromosome, start, end, strand, sequences){
#   result <- DNAString(sequences[[chromosome]][start:end])
#   if (strand=="-")
#     result <- reverseComplement(result)
#   return(result)
# }

# My versions are much slower
# rev_comp_DNA <- function(old_seq){
#   # Create the lookup list
#   RC <- list(A="T", C="G", G="C", T="A")
#   #Split the string into a vector
#   split_seq <- strsplit(x = old_seq, split = "")[[1]]
#   #Run lapply on the vector to swap values, then reverse it, then turn back into string
#   result <- paste(rev(unname(unlist(mclapply(X = split_seq, FUN = function(x) RC[x])))), collapse = '')
#   return(result)
# }
#
# get_genome_sequence <- function(chromosome, start, end, strand, sequences){
#   result <- as.character(sequences[[chromosome]][start:end])
#   if (strand=="-")
#     result <- rev_comp_DNA(result)
#   return(result)
# }
# a <- names(get_sequence(chrom = as.character(seqnames(tm)), start = as.numeric(start(tm)), end = as.numeric(end(tm)), strand = as.character(strand(tm)), genome = genome_sequence2)
# b <- transpose(get_sequence(chrom = as.character(seqnames(tm)), start = as.numeric(start(tm)), end = as.numeric(end(tm)), strand = as.character(strand(tm)), genome = genome_sequence2)
# names(b) <- a
# data.frame(b

