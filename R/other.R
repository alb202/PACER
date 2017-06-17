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

shuffle_intervals <- function(alignments, genome_sizes, n=NULL, regions=NULL, invert=FALSE, seed=0, stranded=TRUE){

  strands <- c("+", "-")

  if (is.null(n))
    n <- length(alignments)
  results <- ""
  for (i in 1:n){

  }

  return(sort.GenomicRanges(results))
}
