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

