get_DNA_sequence <- function(genome, gr, start, end, type=c("5", "3")){

  if (start>end)
    stop("The start position must be prior to the end position")

  plus <- gr[strand(gr)=="+"]
  minus <- gr[strand(gr)=="-"]

  # process each

  # mapply each

  # combine the results

  # return a new GR with the sequences as metadata? or just the data as a vector? I think you need to connect it so you can do later analysis
  }


# One function for the mapply function
# split the plus and minus strand reads and operate on them separately

get_sequence_from_genome <- function(chromosome, start, end){

}
