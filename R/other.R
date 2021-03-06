### -----------------------------
### Other helper functions below

## Create the output directory for a dataset
create_output_dir <- function(out_dir, name){
  # If the main output dir does not exist, create it
  if(!dir.exists(out_dir))
    dir.create(out_dir)
  # Pick the name of the new directory
  new_dir <- paste(out_dir, name, sep = "/")
  # If the new directory doesn't exist, create it
  if(!dir.exists(new_dir))
    dir.create(new_dir)
  # Return the full path to the new directory
  return(new_dir)
}

make_filename <- function(filename, dataset_name){
  return(paste(sub(x = filename,
                   pattern = ".fastq.gz|.fastq.gzip|.fq.gzip|.fq.gz|.fastq|.fq",
                   ignore.case = TRUE,
                   replacement = ""),
               dataset_name, sep = "-"))
}

get_interval_length <- function(x){
  x <- GRanges(x)
  mcols(x)["width"] <- end(x)-start(x)
  return(x)
}


gzip_a_file <- function(dir, file){
  system(command = paste("gzip -f", paste(dir, file, sep = "/"), wait = TRUE))
  return(TRUE)
}

## Invert the locations of values in a vector of two values
invert_vector <- function(x){
  # Get the unique set of values for the input vector
  values <- unique.Vector(x)
  # Store each possible value
  value1 <- values[1]
  value2 <- values[2]
  # Find the locations of value1. NA values must be treated differently
  if (is.na(value1)){
    value1_index <- is.na(x)
  } else{
    value1_index <- x==value1
  }
  # Find the locations of value2. NA values must be treated differently
  if (is.na(value2)){
    value2_index <- is.na(x)
  } else{
    value2_index <- x==value2
  }
  # Swap the values
  x[value1_index] <- value2
  x[value2_index] <- value1
  # Return the new vector
  return(x)
}




## Swap out any number of values of a vector to with different values
swap_values <- function(x, old, new){
  # Initialize a new list to temporarily store values
  tmp <- list()
  # Make sure there are the same number of values to be replaced as there are replacements
  if (length(old)!= length(new))
    stop("The vectors of old and new values must be the same length")
  # Find the locations of the old values
  for(i in 1:length(old)){
    tmp[[i]] <- x==old[i]
  }
  # Replace them with the new values
  for(i in 1:length(old)){
    x[tmp[[i]]] <- new[i]
  }
  # Return the new list
  return(x)
}

## A helper function for shuffle_intervals
# Loops through each combination of alignment properties
loop <- function(seqnames, width, strand, freq, intervals){
  # Get only the inclusion intervals that are relevant for this unique row
  intervals_subset <- intervals[intervals$seqnames==seqnames &
                                  intervals$width>=width &
                                  intervals$strand==strand, ]
  # If there are no intervals to place the shuffled alignment, return a placeholder instead
  if (nrow(intervals_subset) == 0)
    return(matrix(data = c("I", "1", "1", "0", "+"), nrow = 5, ncol = 1, dimnames = list(c("seqnames", "start", "end", "width", "strand"), "NA")))
  # Pick random interval indeces, weighted for the size of the interval, with replacement
  intervals_subset <- intervals_subset[sample(x = 1:nrow(intervals_subset), size = freq, replace = TRUE, prob = intervals_subset$width), ]
  # Use mapply on the random intervals to pick start and end positions within them, create a new dataframe row, and append it to the results
  shuffled_alignments <- mapply(randomize, seqnames, intervals_subset$start, intervals_subset$end, width, strand, SIMPLIFY = TRUE)
  # Return the new alignments
  return(shuffled_alignments)
}

## A helper function for shuffle_intervals
# The randomize function that picks a start and end position within the intervals and returns new intervals
randomize <- function(chrom,start,end,width,strand){
  # Pick the new start position
  new_start <- sample(x = start:(end-width), size = 1)
  # Return the new row
  return((c(seqnames=as.character(chrom), start=new_start, end=new_start+width, width=width, strand=as.character(strand))))
}

## Export a GRange object as a BED file
write_granges_as_BED <- function(gr, filename, directory){
  # Convert the GRange to a dataframe
  df <- data.frame(sort(sortSeqlevels(gr)))
  # Add a score column of zeros
  df$score <- rep(0,nrow(df))
  # Rearrange the columns and export to folder as BED file
  write_delim(x = df[c("seqnames", "start", "end", "width", "score", "strand")],
              path = paste(directory, filename, sep = "/"),
              delim = "\t",
              append = FALSE,
              col_names = FALSE)
}

pretty_base1 <- function(start, end, old=c(0), new=c(1)){
  # Create the initial pretty number sequence
  #result <- pretty(start:end)
  result <- seq(start, end, 10)
  # For each set of values, replace the old with the new
  for(i in 1:length(old))
    result[ result==old[i] ] <- new[i]
  # Return the new sequence
  return(result)
}

get_filename <- function(path, extension=FALSE){
  split_path <- strsplit(x = path, split = "/", fixed = TRUE)[[1]]
  filename <- split_path[[length(split_path)]]
  if(extension==TRUE) return(filename)
  return(strsplit(x = filename, split = ".", fixed = TRUE)[[1]][[1]])
}

get_timestamp <- function()
{
  return(strsplit(x = as.character(Sys.time()), split = " "))
}
