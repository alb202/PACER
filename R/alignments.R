# alignments.R

# trim_cmd <- "cutadapt --max-n=0 -m"
# align_cmd <- "bowtie -n2 -e1000 -l22 -k4 --best --strata -S --threads"
# input_dir <- "/media/ab/HD4/PACER/data/input"
# output_dir <- "/media/ab/HD4/PACER/data/output/SRR1820932-WT_2017-09-10_23-15-00"
# input_file <- "SRR1820932_SAMPLE.fastq.gz"
# adapters_file <- "/media/ab/HD4/PACER/data/adapters/adapters.txt"
# index_dir <- "/media/ab/HD4/PACER/data/indexes/WBcel235"
# genome <- "WBcel235"
# cores <- 2
# dataset_ID <- "SRR1820932_SAMPLE"
# min_length <- 12

align_reads <- function(trim_cmd, align_cmd, input_dir, input_file, output_dir, adapters_file, index_dir, genome, cores, dataset_ID, min_length){
  adapter_string <- paste(" -a", read.delim(file = adapters_file,
                                            strip.white = TRUE,
                                            sep = " ",
                                            header = FALSE,
                                            comment.char = "#",
                                            stringsAsFactors = FALSE)[,1],
                          collapse = "")

  full_cmd <- paste(trim_cmd,
                    min_length,
                    adapter_string,
                    getAbsolutePath(paste(input_dir, input_file, sep="/")),
                    "2>",
                    paste(output_dir, "trim.log", sep = "/"),
                    "|",
                    align_cmd,
                    cores,
                    paste(index_dir, genome, sep = "/"),
                    "- 2>",
                    paste(output_dir, "align.log", sep = "/"),
                    "|",
                    "samtools view -b -@",
                    cores,
                    "|",
                    "samtools sort -O BAM -@", cores,
                    "-o",
                    paste(output_dir, "/", dataset_ID, ".bam", sep = ""),
                    sep = " ")
  print(full_cmd)
  system(command = full_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, wait = TRUE)
  return(paste(output_dir, "/", dataset_ID, ".bam", sep = ""))
  #cutadapt --max-n=0 -m12 -a TGGAATTCTCGGGTGCCAAGG /media/ab/HD4/PACER/data/input/SRR1820932_SAMPLE.fastq.gz 2> /media/ab/HD4/PACER/data/input/cudadapt.log | bowtie -n2 -e1000 -l22 -k4 --best --strata -S --threads 2 /media/ab/HD3/align/bt1/WBcel235 - 2> /media/ab/HD4/PACER/data/input/bowtie.log | samtools view -b -@2| samtools sort -o out.bam -O BAM -@2
}












run_aligner <- function(output_dir, bt_options, cores, index_dir, genome, dataset_ID) {

  #alignment_file <- paste(dataset_info["output_dir"],'/',options_name,"_", dataset_info["basename"],".bam", sep = "")
  #trimmed_fastq_file <- paste(dataset_info["output_dir"],'/',dataset_info["basename"],".trimmed.fastq", sep = "")

  bowtie_args <- c(bt_options, "--threads", cores,
                   paste(index_dir, genome, sep = "/"),
                   paste(output_dir, "/", dataset_ID, ".trimmed.fastq", sep = ""),
                   "| samtools sort - | samtools view -bo",
                   paste(output_dir, "/", dataset_ID, ".bam", sep = ""))#,#)# >",
  #paste(output_dir, "/", dataset_ID, ".bam", sep = ""))

  # paste(dataset_info["output_dir"],
  #       "/bowtie.log",
  #       sep = ""),
  # "| samtools sort - | samtools view -bo",
  # alignment_file,
  # "-")

  # bowtie_cmd <- make_alignment_cmd(
  #   trimmed_fastq_file = trimmed_fastq_file, # Path and name of fastq file
  #   alignment_file = alignment_file, # Path and name of output alignment file
  #   bt_options = bt_options, # Aignment options
  #   genome = genome, # ID of genome
  #   dataset_info = dataset_info) # The dataset names

  # write(x = "\n", file = paste(dataset_info["output_dir"], "bowtie.log", sep = "/"),
  #       sep = "", append = TRUE)
  # write(x = paste(bt_options, options_name, dataset_info["basename"], sep = " | "),
  #       file = paste(dataset_info["output_dir"], "bowtie.log", sep = "/"),
  #       sep = "",
  #       append = TRUE)

  alignment_stderr <- system2(command = "bowtie",
                              args = bowtie_args,
                              #stdout = TRUE,
                              stdout = paste(output_dir, "/", dataset_ID, ".stdout", sep = ""),
                              stderr = paste(output_dir, "/", dataset_ID, ".stderr", sep = ""),
                              #stderr = FALSE,
                              wait = TRUE)
  return(alignment_stderr)
}
#
# make_alignment_cmd <- function(trimmed_fastq_file, alignment_file, bt_options, genome, dataset_info){
#   bowtie_cmd <- c(bt_options,
#                   genome_indexes[genome],
#                   trimmed_fastq_file, "2>>",
#                   paste(dataset_info["output_dir"],
#                         "/bowtie.log",
#                         sep = ""),
#                   "| samtools sort - | samtools view -bo",
#                   alignment_file,
#                   "-")
#   return(bowtie_cmd)
# }

# Shuffle the alignments in a GRanges or GAlignments file within given intervals
shuffle_alignments <- function(alignments, intervals, antisense=FALSE){
  # Turn the alignments into a data frame
  alignments <- data.frame(granges(alignments))
  # Save the column names for later
  col_names <- colnames(alignments)
  # Reduce the intervals and turn into data frame
  intervals <- data.frame(granges(GenomicRanges::reduce(x = intervals, ignore.strand=FALSE)))
  # Create the new alignments antisense to the given intervals
  if (antisense==TRUE)
    intervals$strand <- invert_vector(intervals$strand)

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
                                                           as.numeric(unique_rows$freq),
                                                           MoreArgs = list(intervals = intervals)),
                                                    row.names = NULL,
                                                    fix.empty.names = FALSE,
                                                    stringsAsFactors = FALSE)))
  # Replace the column names and remove the row names from the results
  colnames(shuffled_alignments) <- col_names
  rownames(shuffled_alignments) <- NULL
  # Remove any alignments that didn't have an exon to be placed in
  shuffled_alignments <- subset(shuffled_alignments, !(seqnames=="I" & start==1 & end==1 & width==0 & strand=="+"))
  # Convert the results dataframe to a GRange, sort and return it
  return(sort.GenomicRanges(
    x = makeGRangesFromDataFrame(df = shuffled_alignments,
                                 seqnames.field = "seqnames",
                                 start.field = "start",
                                 end.field = "end",
                                 strand.field = "strand")))
}

