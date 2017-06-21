# run_settings <- function(){
#   print("settings")
#   adapter_file <- paste(getwd(), "/adapters/adapters.txt", sep="")
#   input_dir <- paste(getwd(), "/raw_data", sep="")
#   datasets <- c(WT_early_rep1="SRR5023999_100K_sample.fastq.gz")
#   output_dir <- paste(getwd(), "/output", sep="")
#   sam_files <- c(two_mismatch="-v2 -k4 --best -S",
#                  no_mismatch="-v0 -k4 --best -S",
#                  no_seed_mismatch="-n0 -e1000 -l22 -k4 --best -S")
#   genome = "ce10"
# }

# Other helper functions below
create_output_dirs <- function(name){
  if(!dir.exists(output_dir))
    dir.create(output_dir)
  new_dir <- paste(output_dir, name, sep = "/")
  if(!dir.exists(new_dir))
    dir.create(new_dir)
  return(new_dir)
}

invert_vector <- function(x){
  values <- unique.Vector(x)
  value1 <- values[1]
  value2 <- values[2]
  if (is.na(value1)){
    value1_index <- is.na(x)
  } else{
    value1_index <- x==value1
  }
  if (is.na(value2)){
    value2_index <- is.na(x)
  } else{
    value2_index <- x==value2
  }
  x[value1_index] <- value2
  x[value2_index] <- value1
  return(x)
}

swap_values <- function(x, old, new){
  tmp <- list()
  if (length(old)!= length(new))
    stop("The vectors of old and new values must be the same length")
  for(i in 1:length(old)){
    tmp[[i]] <- x==old[i]
  }
  for(i in 1:length(old)){
    x[tmp[[i]]] <- new[i]
  }
  return(x)
}



shuffle_intervals <- function(alignments, chromosome_sizes, sample=NULL, include=NULL, seed=0, stranded=TRUE, chrom=TRUE){
  #eval(as.name("II"))
  strands <- c("+", "-")
  chromosomes <- levels(chromosome_sizes$chrom)
  if (!is.null(include)){
    if(class(include)[1] != "GRanges")
      stop("'include' must be a GenomicRanges object")
    for (i in chromosomes){
      assign(x = i, value = include[seqnames(include)==i])
    }
  } else {
    for (i in chromosomes){
        assign(x = i,
               value = GRanges(seqnames = rep(i,2),
                               ranges = rep(IRanges(start = 1, end = chromosome_sizes$length[chromosome_sizes$chrom==i]),2),
                               strand = strands))
    }
  }


  if (!is.null(sample)){
    if (!is.numeric(sample))
      stop("'sample' must be numeric")
    if (sample<=0)
      stop("'sample' must be a positive number")
    if (sample<1){
      sample <- length(x = alignments)*sample
    }
    alignments <- sample(x = alignments, size = sample, replace = FALSE)
  }
  if (stranded==FALSE)
    strand_probabilities <- c(sum(as.numeric(strand(two_mismatches_filtered)=="+")), sum(as.numeric(strand(two_mismatches_filtered)=="-")))
  n <- length(alignments)
  results <- GRanges()
  for (i in 1:n){
    print(i)
    row <- alignments[i]
    if (chrom==FALSE){
      new_chrom <- as.character(sample(x = as.character(chromosome_sizes$chrom), size = 1, prob = chromosome_sizes$length))
    } else {
      new_chrom <- as.character(seqnames(row))
    }
    if (stranded==FALSE){
      new_strand <- sample(x = strands, size = 1, prob = strand_probabilities)
    } else {
      new_strand <- as.character(strand(row))
    }
    width <- qwidth(alignments[i])
    #eval(as.name(new_chrom))
    intervals <- eval(as.name(new_chrom))[as.character(strand(eval(as.name(new_chrom))))==new_strand]
    new_interval <- sample(intervals, size = 1)
    tries <- 0
    repeat
    {
      new_interval <- sample(intervals, size = 1)
      if ((end(new_interval)-start(new_interval))>=width)
        break
      tries <- tries+1
      if (tries == 50)
        break
    }
    if (tries == 50)
      next
    new_start <- sample(x = start(new_interval):(end(new_interval)-width), size = 1)
    results <- c(results, makeGRangesFromDataFrame(data.frame(seqnames=new_chrom, start=new_start, end=new_start+width, strand=new_strand)))
    #a <- c(a,makeGRangesFromDataFrame(data.frame(seqnames="IV", start="438482", end="594998", strand="-")))
  }
  return(sort(sortSeqlevels(x = results)))
}

