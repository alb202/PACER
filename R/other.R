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

shuffle_intervals <- function(alignments, intervals){
  # Turn the alignments into a data frame
  alignments <- data.frame(granges(alignments))
  # Save the column names for later
  col_names <- colnames(alignments)
  # Reduce the intervals and turn into data frame
  inclusive <- data.frame(granges(GenomicRanges::reduce(x = intervals, ignore.strand=FALSE)))
  # Remove the old start and end columns from the alignments
  alignments["start"] <- NULL
  alignments["end"] <- NULL
  # Get the counts of each combination of chromosome, strand and width
  unique_rows <- plyr::count(alignments)
  # Run mapply on every unique chromosome/strand/width and
  # save results to the shuffled alignments object
  shuffled_alignments <- as_data_frame(t(data.frame(mapply(FUN = loop,
                                             as.character(unique_rows$seqnames),
                                             as.numeric(unique_rows$width),
                                             as.character(unique_rows$strand),
                                             as.numeric(unique_rows$freq)),
                                      row.names = NULL, fix.empty.names = FALSE, stringsAsFactors = FALSE)))
  # Replace the column names and remove the row names from the results
  colnames(shuffled_alignments) <- col_names
  rownames(shuffled_alignments) <- NULL
  # Remove any alignments that didn't have an exon to be placed in
  shuffled_alignments <- subset(shuffled_alignments, !(seqnames=="I" & start==1 & end==1 & width==0 & strand=="+"))
  # Convert the results dataframe to a GRange, sort and return it
  return(sort(sortSeqlevels(
    x = makeGRangesFromDataFrame(df = shuffled_alignments,
                                 seqnames.field = "seqnames",
                                 start.field = "start",
                                 end.field = "end",
                                 strand.field = "strand"))))
}

loop <- function(seqnames, width, strand, freq){
  # Get only the inclusion intervals that are relevant for this unique row
  intervals <- inclusive[inclusive$seqnames==seqnames &
                           inclusive$width>=width &
                           inclusive$strand==strand, ]
  # If there are no intervals to place the shuffled alignment, return a placeholder instead
  if (nrow(intervals) == 0)
    return(matrix(data = c("I", "1", "1", "0", "+"), nrow = 5, ncol = 1, dimnames = list(c("seqnames", "start", "end", "width", "strand"), "NA")))
  # Pick random interval indeces, weighted for the size of the interval, with replacement
    intervals <- intervals[sample(x = 1:nrow(intervals), size = freq, replace = TRUE, prob = intervals$width), ]
  # Use mapply on the random intervals to pick start and end positions within them, create a new dataframe row, and append it to the results
  shuffled_alignments <- mapply(randomize, seqnames, intervals$start, intervals$end, width, strand, SIMPLIFY = TRUE)
  # Return the new alignments
  return(shuffled_alignments)
}

# The randomize function that picks a start and end position within the intervals and returns a new position
randomize <- function(chrom,start,end,width,strand){
  # Pick the new start position
  new_start <- sample(x = start:(end-width), size = 1)
  # Return the new row
  return((c(seqnames=as.character(chrom), start=new_start, end=new_start+width, width=width, strand=as.character(strand))))
}

writeGRangesAsBED <- function(gr, filename, directory){
  df <- data.frame(gr)
  df$score <- rep(0,nrow(df))
  write_delim(x = df[c("seqnames", "start", "end", "width", "score", "strand")],
              path = paste(directory, filename, sep = "/"),
              delim = "\t",
              append = FALSE,
              col_names = FALSE)
}


# bedshuffle_compiled <- cmpfun(bedshuffle)
#
# # Redesign idea: get the exons and genes first, export them to the genome directory, then run bedtools on the raw BAM, then filter later
# # Or just do a random shuffle, then only get the ones that overlap exons. Thats going to be only 2% though, and will be too dispersed.
# # use bedtools bamtobed
# bedshuffle <- function(alignments, chromosome_sizes, samp=NULL, include=NULL, seed=0, chrom=TRUE, max_tries=10){
#   # The alignements must be a GenomicRanges or GenomicAlignments object
#   if(class(alignments)[1] != "GAlignments" & class(alignments)[1] != "GRanges")
#     stop("'alignments' must be a GenomicRanges or GenomicAlignments object")
#   mcols(alignments) <- NULL
#
#   #if there is a sampling number included, take a sample of the alignments
#   if (!is.null(samp)){
#     if (!is.numeric(samp))
#       stop("'sample' must be numeric")
#     if (samp<=0)
#       stop("'sample' must be a positive number")
#     if (samp<1){
#       samp <- length(x = alignments)*samp
#     }
#     alignments <- sample(x = alignments, size = samp, replace = FALSE)
#   }
#
#   # If there are inclusion intervals, set the chromosome intervals using those
#   if (!is.null(include)){
#       if(class(include)[1] != "GRanges")
#         stop("'include' must be a GenomicRanges object")
#       mcols(include) <- NULL
#       include_cmd <- paste("-incl include.bed")
#   } else
#   {
#     include_cmd <- ""
#   }
#   # If the alignments should stay on the same chromosome
#   if (isTRUE(chrom)){
#     chrom_cmd <- "-chrom"
#   } else{
#     chrom_cmd <- ""
#   }
#   # Make the headers
#   headers <- c("seqnames","start","end","strand")
#   # Write the genome sizes
#   write.table(data.frame(chromosome_sizes), file=paste(output_dir,"sizes.bed", sep="/"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
#   # Write the exon intervals
#   write.table(data.frame(unique(include))[headers], file=paste(output_dir,"include.bed", sep="/"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
#   # Write the alignments
#   write.table(data.frame(alignments)[headers], file=paste(output_dir,"alignments.bed", sep="/"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
#   # Create the arguments string
#   shuffle_cmd <- paste("shuffle -i alignments.bed -g sizes.bed", chrom_cmd ,include_cmd, "-maxTries", max_tries, "-seed", seed, " | bedtools sort -i > shuffled.bed")
#   # Set the output directory as the working directory
#   oldwd <- getwd()
#   setwd(output_dir)
#   output <- system2(command = "bedtools", args = shuffle_cmd, wait = TRUE)
#   results <- read.delim(file = paste(output_dir,"shuffled.bed", sep = "/"), col.names = headers, sep = "\t")
#   setwd(oldwd)
#   return(makeGRangesFromDataFrame(results))
# }
#   require(compiler)
  # #enableJIT(3)
  # shuffle_intervals_compiled <- cmpfun(shuffle_intervals)

# shuffle_intervals <- function(alignments, chromosome_sizes, samp=NULL, include=NULL, seed=0, stranded=TRUE, chrom=TRUE, max_tries=20){
#   # Setup the strand and chromosome possibilities
#   strands <- c("+", "-")
#   chromosomes <- levels(chromosome_sizes$chrom)
#
#   # If there are inclusion intervals, set the chromosome intervals using those
#   if (!is.null(include)){
#     if(class(include)[1] != "GRanges")
#       stop("'include' must be a GenomicRanges object")
#     for (i in chromosomes){
#       assign(x = i, value = include[seqnames(include)==i])
#     }
#   }
#   # If there are no inclusion intervals specified, use the entire chromosome and create an interval for both strands
#   else {
#     for (i in chromosomes){
#       assign(x = i,
#              value = GRanges(seqnames = rep(i,2),
#                              ranges = rep(IRanges(start = 1, end = chromosome_sizes$length[chromosome_sizes$chrom==i]),2),
#                              strand = strands))
#     }
#   }
#   # if there is a sampling number included, take a sample of the alignments
#   if (!is.null(samp)){
#     if (!is.numeric(samp))
#       stop("'sample' must be numeric")
#     if (samp<=0)
#       stop("'sample' must be a positive number")
#     if (samp<1){
#       samp <- length(x = alignments)*samp
#     }
#     alignments <- sample(x = alignments, size = samp, replace = FALSE)
#   }
#   # If the strands can be ignored, make the strand possibilities 50/50
#   if (stranded==FALSE)
#     strand_probabilities <- c(.5,.5)
#   # Get the number of alignments to shuffle
#   n <- length(alignments)
#   # Initialized a GenomicRanges object to store the new alignments
#   results <- GRanges()
#   # Loop through the alignments and shuffle them
#
#
#   if (chrom==FALSE){
#     new_chroms <- as.character(sample(x = as.character(chromosome_sizes$chrom), size = n, replace = TRUE, prob = chromosome_sizes$length))
#     }
#   # if (stranded==FALSE){
#   #   new_strands <- sample(x = strands, size = n, replace = TRUE, prob = strand_probabilities)
#   # }
#
#   if (chrom==FALSE & stranded==FALSE){
#     for (i in 1:n){
#       print(i)
#       # Get the alignment
#       #row <- alignments[i]
#       width <- qwidth(alignments[i])
#       #intervals <- eval(as.name(new_chroms[i]))[as.character(strand(eval(as.name(new_chrom))))==new_strand]
#       tries <- 0
#       repeat
#       {
#         new_interval <- sample(eval(as.name(new_chroms[i])), size = 1, replace = TRUE)
#         if ((end(new_interval)-start(new_interval))>=width)
#           break
#         tries <- tries+1
#         if (tries == max_tries)
#           next
#       }
#       # if (tries == max_tries)
#       #   next
#       new_start <- sample(x = start(new_interval):(end(new_interval)-width), size = 1)
#       results <- c(results, makeGRangesFromDataFrame(data.frame(seqnames=new_chroms[i], start=new_start, end=new_start+width, strand=strand(new_interval))))
#       #a <- c(a,makeGRangesFromDataFrame(data.frame(seqnames="IV", start="438482", end="594998", strand="-")))
#     }
#   }
######
#     else if (chrom==FALSE & stranded==TRUE){
#   for (i in 1:n){
#     print(i)
#     # Get the alignment
#     row <- alignments[i]
#     if (chrom==FALSE){
#       new_chrom <- as.character(sample(x = as.character(chromosome_sizes$chrom), size = 1, prob = chromosome_sizes$length))
#     } else {
#       new_chrom <- as.character(seqnames(row))
#     }
#     if (stranded==FALSE){
#       new_strand <- sample(x = strands, size = 1, prob = strand_probabilities)
#     } else {
#       new_strand <- as.character(strand(row))
#     }
#     width <- qwidth(alignments[i])
#     #eval(as.name(new_chrom))
#     intervals <- eval(as.name(new_chrom))[as.character(strand(eval(as.name(new_chrom))))==new_strand]
#     new_interval <- sample(intervals, size = 1)
#     tries <- 0
#     repeat
#     {
#       new_interval <- sample(intervals, size = 1)
#       if ((end(new_interval)-start(new_interval))>=width)
#         break
#       tries <- tries+1
#       if (tries == max_tries)
#         break
#     }
#     if (tries == max_tries)
#       next
#     new_start <- sample(x = start(new_interval):(end(new_interval)-width), size = 1)
#     results <- c(results, makeGRangesFromDataFrame(data.frame(seqnames=new_chrom, start=new_start, end=new_start+width, strand=new_strand)))
#     #a <- c(a,makeGRangesFromDataFrame(data.frame(seqnames="IV", start="438482", end="594998", strand="-")))
#   }
#
#   else if (chrom==TRUE & stranded==FALSE){
#     for (i in 1:n){
#       print(i)
#       # Get the alignment
#       row <- alignments[i]
#       if (chrom==FALSE){
#         new_chrom <- as.character(sample(x = as.character(chromosome_sizes$chrom), size = 1, prob = chromosome_sizes$length))
#       } else {
#         new_chrom <- as.character(seqnames(row))
#       }
#       if (stranded==FALSE){
#         new_strand <- sample(x = strands, size = 1, prob = strand_probabilities)
#       } else {
#         new_strand <- as.character(strand(row))
#       }
#       width <- qwidth(alignments[i])
#       #eval(as.name(new_chrom))
#       intervals <- eval(as.name(new_chrom))[as.character(strand(eval(as.name(new_chrom))))==new_strand]
#       new_interval <- sample(intervals, size = 1)
#       tries <- 0
#       repeat
#       {
#         new_interval <- sample(intervals, size = 1)
#         if ((end(new_interval)-start(new_interval))>=width)
#           break
#         tries <- tries+1
#         if (tries == max_tries)
#           break
#       }
#       if (tries == max_tries)
#         next
#       new_start <- sample(x = start(new_interval):(end(new_interval)-width), size = 1)
#       results <- c(results, makeGRangesFromDataFrame(data.frame(seqnames=new_chrom, start=new_start, end=new_start+width, strand=new_strand)))
#       #a <- c(a,makeGRangesFromDataFrame(data.frame(seqnames="IV", start="438482", end="594998", strand="-")))
#     }
#    return(sort(sortSeqlevels(x = results)))
#}

