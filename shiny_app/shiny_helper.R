load_directory <- function(dir, filetypes, include_dirs=FALSE){
  return(dir(path = dir,
             recursive = FALSE,
             pattern = filetypes,
             include.dirs = include_dirs,
             full.names = FALSE))
}

save_adapters <- function(x, path){
  # If an adapter file already exists, save a backup file
  if(file.exists(path)){
    file.copy(from = path,
              to = paste(path, ".backup", sep = ""),
              overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
  write_delim(x = x,
              path = path,
              delim = " #",
              col_names = FALSE,
              append = FALSE)
  return(file.exists(path))
}

make_choices <- function(choices, numbered=TRUE){
  if(isTRUE(numbered)){
    return(mapply(paste(1:length(choices), choices, sep = " - "), 1:length(choices), FUN = function(x, y){x=y}))
  }else{
    return(mapply(choices, 1:length(choices), FUN = function(x, y){x=y}))
  }
}

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

update_genome_index <- function(session, genomes){
  if(nrow(genomes)>0){
    row.names(genomes) <- 1:nrow(genomes)
    row_choices <- make_choices(genomes$description,
                                numbered = TRUE)
    toggleState(id = "load_genome", condition = TRUE)
    toggleState(id = "view_genome", condition = TRUE)
  } else {
    row_choices <- c(NULL)
    toggleState(id = "load_genome", condition = FALSE)
    toggleState(id = "view_genome", condition = FALSE)
  }
  print("row_choices")
  print(row_choices)
  updateSelectInput(session = session,
                    inputId = "genome_index",
                    label = NULL,
                    choices = row_choices,
                    selected = NULL)
}

update_ensembl_genome_index <- function(session, genomes){
  if(nrow(genomes)>0){
    row.names(genomes) <- 1:nrow(genomes)
    row_choices <- make_choices(genomes$description,
                                numbered = FALSE)
    toggleState(id = "add_genome", condition = TRUE)
  } else {
    row_choices <- c("")
    toggleState(id = "add_genome", condition = FALSE)
  }
  updateSelectInput(session = session,
                    inputId = "ensembl_genome_index",
                    label = NULL,
                    choices = row_choices,
                    selected = NULL)
}

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
