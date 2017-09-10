# Get the list of files in a directory
load_directory <- function(dir, filetypes, include_dirs=FALSE){
  return(dir(path = dir,
             recursive = FALSE,
             pattern = filetypes,
             include.dirs = include_dirs,
             full.names = FALSE))
}

# Save the adapters or the genome list
save_info <- function(x, path, delimeter, col_names){
  # If a file already exists, save a backup file
  if(file.exists(path)){
    file.copy(from = path,
              to = paste(path, ".backup", sep = ""),
              overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
  write_delim(x = x,
              path = path,
              delim = delimeter,
              col_names = col_names,
              append = FALSE)
  return(file.exists(path))
}

# Combine the numbering and the text to create a list of options for a dropdown
make_choices <- function(choices, numbered=TRUE){
  if(isTRUE(numbered)){
    return(mapply(paste(1:length(choices), choices, sep = " - "), 1:length(choices), FUN = function(x, y){x=y}))
  }else{
    return(mapply(choices, 1:length(choices), FUN = function(x, y){x=y}))
  }
}

# Update the index for the adapter modal
update_adapter_index <- function(session, adapters){
  if(nrow(adapters)>0){
    row.names(adapters) <- 1:nrow(adapters)
    row_choices <- make_choices(adapters$Description)
    toggleState(id = "remove_adapter", condition = TRUE)
  } else{
    row_choices <- c("")
    toggleState(id = "remove_adapter", condition = FALSE)
  }
  updateSelectInput(session = session,
                    inputId = "adapter_index",
                    label = NULL,
                    choices = row_choices,
                    selected = NULL)
}

# Update the index for the genome list on the genome modal
update_genome_index <- function(session, genomes){
  # Toggle buttons based on the number of items in the list
  if(nrow(genomes)>0){
    print("more than 0 genomes")
    row.names(genomes) <- 1:nrow(genomes)
    row_choices <- make_choices(genomes[,2],
                                numbered = TRUE)
    toggleState(id = "load_genome", condition = TRUE)
    toggleState(id = "remove_genome", condition = TRUE)
    toggleState(id = "view_genome", condition = TRUE)
  } else {
    print("0 genomes")
    row_choices <- c(" "=0)
    toggleState(id = "load_genome", condition = FALSE)
    toggleState(id = "remove_genome", condition = FALSE)
    toggleState(id = "view_genome", condition = FALSE)
    toggleElement(id = "genome_details", condition = FALSE)
  }
  # Update the selector
  updateSelectInput(session = session,
                    inputId = "genome_index",
                    label = NULL,
                    choices = row_choices,
                    selected = NULL)
}
# Update the ensembl genome list selector
update_ensembl_genome_index <- function(session, genomes){
  # Toggle buttons based on the number of options
  if(nrow(genomes)>0){
    row.names(genomes) <- 1:nrow(genomes)
    row_choices <- make_choices(genomes[,2],
                                numbered = FALSE)
    toggleState(id = "add_genome", condition = TRUE)
  } else {
    row_choices <- c("")
    toggleState(id = "add_genome", condition = FALSE)
  }
  # Update the dropdown list of ensembl genomes
  updateSelectInput(session = session,
                    inputId = "ensembl_genome_index",
                    label = NULL,
                    choices = row_choices,
                    selected = NULL)
}

# Check if a genome entry is complete
check_for_complete_genome <- function(genome){
  if(nrow(genome)==0) return(genome)
  result <- rep("Incomplete", nrow(genome))
  for(i in 1:nrow(genome)){
    if(genome[i, 5] != "None" &
       genome[i, 6] != "None" &
       genome[i, 7] != "None") result[i] <- "Ready"
  }
  return(result)
}

# Add a genome from the ensembl list to the genome list
add_genome <- function(genomes, new_genome){
  print(names(genomes))
  print(names(new_genome))
  names(new_genome) <- c("Dataset", "Description", "Version")
  if(nrow(genomes)==0){
    print("create new genome list")
    genomes <- data.frame(stringsAsFactors = FALSE, row.names = NULL, new_genome, "Status"="Incomplete", "Interval.Status"="None", "Genome.FASTA"="None", "Bowtie.Index"="None", "Gene Sets"="None")
  } else {
    print("append genome to list")
    genomes <- rbind(genomes, data.frame(stringsAsFactors = FALSE, row.names = NULL, new_genome, "Status"="Incomplete", "Interval.Status"="None", "Genome.FASTA"="None", "Bowtie.Index"="None", "Gene.Sets"="None"))
  }
  return(genomes)
}
