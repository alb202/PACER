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
