#server.R
library(shiny)
library(shinyFiles)
library(shinyjs)
library(DT)
library(R.utils)

source("shiny_helper.R")
source("modals.R")
source("../R/load.R")

server <- function(input, output, session){

  values <- reactiveValues()
  observeEvent(eventExpr = values, once = TRUE, {
    # Set the initial value of the input directory to the first entry in the volumes list
    values$input_volumes <- c("Input"="../data/input","Home"="~/", "Root"="/")
    values$output_volumes <- c("Output"="../data/output","Home"="~/", "Root"="/")
    values$genome_volumes <- c("Data"="../data/genomes","Home"="~/", "Root"="/")

    #  values$genomes <- read.delim(file = "../data/genomes/genomes_new.txt", stringsAsFactors = FALSE, header = FALSE, sep = ",", row.names= NULL, col.names = c("Description", "Version", "Dataset", "Status", "Interval.Status", "Genome.FASTA", "Bowtie.Index", "Gene.Sets"))
    values$genomes_dir <- "../data/genomes/"
    values$adapters_dir <- "../data/adapters/"
    values$genomes_file <- "genomes.txt"
    values$adapters_file <- "adapters.txt"
    values$adapters <- read.delim(file = paste(values$adapters_dir, values$adapters_file, sep = ""), stringsAsFactors = FALSE, header = FALSE, sep = "#", col.names = c("Adapter", "Description"))
    values$genomes <- read.delim(file = paste(values$genomes_dir, values$genomes_file, sep = ""), stringsAsFactors = FALSE, header = TRUE, sep = "\t", row.names= NULL, col.names = c("Dataset", "Description", "Version", "Status", "Interval.Status", "Genome.FASTA", "Bowtie.Index", "Gene.Sets"))

    # Output directories
    values$input_dir <- values$input_volumes[[1]]
    values$output_dir <- values$output_volumes[[1]]


  })

  observe({
    # Only show the "Align dataset" button if a dataset has been selected
    toggleState(id = "use_dataset", condition = (!is.null(input$input_file) & !is.null(values$selected_genome)))
  })

  observe({
    # Only show the "View dataset" button if a dataset has been selected
    toggleState(id = "view_dataset", condition = !is.null(input$processed_file))
  })

  output$input_dir <- renderText({
    # Print the dataset directory
    print("render directory name")
    return(getAbsolutePath(values$input_dir))
  })

  output$output_dir <- renderText({
    # Print the Output directory
    print("render output directory name")
    return(values$output_dir)
  })

  observe({
    print("getting input directory")
    # When the directory button is pressed, these run
    shinyDirChoose(input = input,
                   id = "find_input_dir",
                   roots=isolate(values$input_volumes),
                   session=session
    )
    # Parse the result and save the output directory
    if(!is.null(input$find_input_dir))
      values$input_dir <- parseDirPath(roots = values$input_volumes,
                                       selection = input$find_input_dir)
    # Print the output directories
    print(values$input_dir)
  })

  observe({
    print("getting output directory")
    # When the directory button is pressed, these run
    shinyDirChoose(input = input,
                   id = "find_output_dir",
                   roots=isolate(values$output_volumes),
                   session=session
    )
    # Parse the result and save the output directory
    if(!is.null(input$find_output_dir))
      values$output_dir <- parseDirPath(roots = values$output_volumes,
                                        selection = input$find_output_dir)
    # Print the output directories
    print(values$output_dir)
  })

  observeEvent(values$input_dir, {
    # Update the file select box with the list of fastq files from the directory
    print("Update the file input box")
    updateSelectInput(session = session,
                      inputId = "input_file",
                      label = NULL,
                      choices = load_directory(dir = values$input_dir,
                                               filetypes = c("fastq", "fq")),
                      selected = NULL)
  })

  observeEvent(values$output_dir, {
    # Update the file select box with the list of processed datasets
    print("Update the output box")
    updateSelectInput(session = session,
                      inputId = "processed_file",
                      label = NULL,
                      choices = load_directory(dir = values$output_dir,
                                               filetypes = "",
                                               include_dirs = TRUE),
                      selected = NULL)
  })

  observeEvent(input$use_dataset, {
    print("Open the start modal")
    input_file <- input$input_file
    showModal(start_modal)
    output$start_modal_title <- renderText({paste("File: ", input$input_file, sep = "")})
  })

  ### Genome selector modal observers
  observeEvent(label = "Select genome modal",
               eventExpr = input$select_genome,
               priority = 100, {
                 print("Select genome")
                 showModal(genomes_modal)
                 output$genome_table <- renderTable(width = "100%" ,
                                                    rownames = TRUE,
                                                    spacing = "s",
                                                    expr =  {
                                                      new_table <- values$genomes
                                                      if(nrow(new_table)>0){
                                                        print("Changing values of new_table to OK")
                                                        new_table[,5:8][new_table[,5:8]!="None"] <- "OK"
                                                      }
                                                      return(new_table)
                                                    })
                 print("Update the genome index")
                 update_genome_index(session = session,
                                     genomes = values$genomes)
               })

  observeEvent(eventExpr = input$view_genome,
               label = "View genome details",
               priority = 0, {
                 print("View the genome details")
                 toggleElement(id = "genome_details")
               })

  observeEvent(eventExpr = c(values$genomes[input$genome_index, 5:8],
                             input$view_genome),
               label = "Update genome details view",
               priority = -5, {
                 print("Update the genome view")
                 if(isolate(nrow(values$genomes) > 0)){
                   # Get the genome row
                   genome_row <- values$genomes[input$genome_index, ]

                   # Indicate if the genome intervals are ready
                   output$genome_interval_status <- renderUI({
                     if(as.character(genome_row[5])=="None"){
                       return(HTML("<font color='red'>Incomplete</font>"))
                     } else {
                       return(HTML("<font color='green'>Ready</font>"))
                     }
                   })

                   # Check if gene sets are loaded
                   if(nrow(genome_row)>0){
                     gene_list_choices <- make_choices(choices = strsplit(x = genome_row[["Gene.Sets"]],
                                                                          split = ";",
                                                                          fixed = TRUE)[[1]],
                                                       numbered = FALSE)
                   }
                   else{
                     gene_list_choices <- c(" "= 0)}
                   updateSelectInput(session = session,
                                     inputId = "gene_list_status",
                                     choices = gene_list_choices)

                   # Check if the genome fasta location is set
                   output$genome_fasta_location <- renderUI({
                     if(as.character(genome_row[["Genome.FASTA"]])=="None"){
                       return(HTML("<font color='red'>Incomplete</font>"))
                     } else {
                       return(HTML("<font color='green'>",genome_row[["Genome.FASTA"]],"</font>"))
                     }
                   })

                   # Check if the genome index location is set
                   output$genome_index_location <- renderUI({
                     if(as.character(genome_row[["Bowtie.Index"]])=="None"){
                       return(HTML("<font color='red'>Incomplete</font>"))
                     } else {
                       return(HTML("<font color='green'>",genome_row[["Bowtie.Index"]],"</font>"))
                     }
                   })

                   # Find the genome fasta
                   shinyFileChoose(input = input,
                                   session = session,
                                   id = "genome_fasta_finder",
                                   roots = isolate(values$genome_volumes),
                                   filetypes = c("fasta", "fa"))

                   # Find the genome index
                   shinyDirChoose(input = input,
                                  session = session,
                                  id = "genome_index_finder",
                                  roots = isolate(values$genome_volumes))

                   # Find gene lists
                   shinyFileChoose(input = input,
                                   session = session,
                                   id = "gene_list_finder",
                                   roots = isolate(values$genome_volumes),
                                   filetypes = c("txt", "tsv"))
                 }
               })

  observeEvent(eventExpr =  input$genome_fasta_finder,
               label = "Fasta finder", {
                 print("Get the genome fasta file location")
                 if(!is.null(input$genome_fasta_finder)){
                   fasta_location <- parseFilePaths(roots = values$genome_volumes, selection = input$genome_fasta_finder)
                   values$genomes[isolate(input$genome_index), "Genome.FASTA"] <- getAbsolutePath(as.character(fasta_location[["datapath"]]))
                 }
               })

  observeEvent(eventExpr = input$genome_index_finder,
               label = "Index finder", {
                 print("Get the genome index file location")
                 if(!is.null(input$genome_index_finder)){
                   index_location <- parseDirPath(roots = values$genome_volumes, selection = input$genome_index_finder)
                   #updateTextInput(session = session, inputId = "genome_index_finder", value = NULL)
                   values$genomes[isolate(input$genome_index), "Bowtie.Index"] <- getAbsolutePath(as.character(index_location))
                 }
               })

  observeEvent(eventExpr = values$genomes,
               label = "Checking complete genome",
               priority = -10, {
                 print("Check if genome is complete")
                 if(nrow(values$genomes)>0){
                   values$genomes[,"Status"] <- check_for_complete_genome(genome = values$genomes)
                 }
               })

  observeEvent(eventExpr = input$gene_list_finder,
               label = "Gene list finder", {
                 print("Get list of genes")
                 if(!is.null(input$gene_list_finder)){
                   gene_list_location <- parseFilePaths(roots = values$genome_volumes, selection = input$gene_list_finder)
                   isolate(
                     if(values$genomes[input$genome_index, "Gene.Sets"]=="None"){
                       values$genomes[[input$genome_index, "Gene.Sets"]] <- getAbsolutePath(gene_list_location[["datapath"]])
                     } else {
                       if(getAbsolutePath(gene_list_location[["datapath"]]) %in% strsplit(x = values$genomes[input$genome_index, "Gene.Sets"],
                                                                                          split = ";",
                                                                                          fixed = TRUE)[[1]]) return(NULL)
                       values$genomes[[input$genome_index, "Gene.Sets"]] <- paste(values$genomes[input$genome_index, "Gene.Sets"],
                                                                                  getAbsolutePath(gene_list_location[["datapath"]]),
                                                                                  sep = ";")
                     }
                   )
                 }
               })

  observeEvent(eventExpr = input$remove_gene_list,
               label = "Remove gene list", {
                 print("Remove a gene list")
                 isolate(genome_row <- values$genomes[input$genome_index, ])
                 isolate(new_list <- paste(strsplit(x = genome_row[["Gene.Sets"]],
                                                    split = ";",
                                                    fixed = TRUE)[[1]][-(as.integer(input$gene_list_status))], collapse = ";"))
                 if(new_list=="") new_list <- "None"
                 isolate(values$genomes[input$genome_index, "Gene.Sets"] <- new_list)
               })


  observeEvent(eventExpr = input$get_ensembl_genomes,
               label = "Get ensembl genomes", {
                 print("Get the list of ensembl datasets")

                 # Create a Progress object
                 values$progress <- shiny::Progress$new()
                 values$progress$set(message = "Getting datasets: ", value = 0)
                 on.exit(values$progress$close())

                 # Get the list of ensembl datasets
                 values$ensembl_genome_index <- get_datasets_from_biomart(updateProgress = updateProgress)
                 update_ensembl_genome_index(session = session, genomes = values$ensembl_genome_index)
               })

  observeEvent(eventExpr = input$add_genome,
               label = "Add genome", {
                 if(values$ensembl_genome_index[as.integer(input$ensembl_genome_index), "dataset"] %in% values$genomes[, "Dataset"])
                   return(NULL)
                 values$genomes <- add_genome(genomes = values$genomes, new_genome = values$ensembl_genome_index[as.integer(input$ensembl_genome_index), ])
                 update_genome_index(session = session, genomes = values$genomes)
               })

  observeEvent(eventExpr = input$remove_genome,
               label = "Remove genome", {
                 print("Remove a genome")
                 if(!is.na(input$genome_index)){
                   values$genomes <- values$genomes[-as.integer(input$genome_index), ]
                   if(nrow(values$genomes)>0)
                     row.names(values$genomes) <- 1:nrow(values$genomes)
                 }
                 update_genome_index(session = session, genomes = values$genomes)
               })

  observeEvent(eventExpr = input$get_intervals,
               label = "Get intervals", {
                 print("Get gene intervals from ensembl")
                 genome <- values$genomes[isolate(input$genome_index), 3]
                 if(!check_for_intervals(path = values$genomes_dir, genome = genome)){
                   # Create a Progress object
                   values$progress <- shiny::Progress$new()
                   values$progress$set(message = "Getting Genome: ", value = 0)
                   on.exit(values$progress$close())

                   # Getting ensembl intervals
                   get_ensembl_intervals(path = values$genomes_dir, genome = genome, updateProgress = updateProgress)
                 }
                 values$genomes[isolate(input$genome_index), 5] <- "OK"
               })

  observeEvent(eventExpr = input$save_genomes,
               label = "Save genomes", {
                 print("Save the genomes list")
                 success <- save_info(x = values$genomes,
                                      path = paste(values$genomes_dir,
                                                   values$genomes_file,
                                                   sep = ""),
                                      delimeter = "\t",
                                      col_names = TRUE)
                 if(isTRUE(success)){
                   output$genome_changes_saved <- renderText(expr = {return("Saved")})
                 }
               })

  # Clear the indicator after saving
  observeEvent(values$genomes, {output$genome_changes_saved <- renderText({" "})})

  observeEvent(eventExpr = input$load_genome,
               label = "Load genome", {
                 print("Load the selected genome")
                 selected_genome <- values$genomes[input$genome_index, ]
                 if(selected_genome["Status"] == "Ready"){
                   # If the genome is ready
                   values$selected_genome <- selected_genome
                   output$selected_genome <- renderText({return(selected_genome[["Description"]])})
                 } else {
                   # If the genome isn't ready
                   values$selected_genome <- NULL
                   output$selected_genome <- renderText({"Select a complete genome!"})
                 }
               })


  ### Adapter selector modal observers
  observeEvent(eventExpr = input$view_adapters,
               label = "View adapters", {
                 print("Viewing adapters modal")
                 showModal(adapter_modal)
                 output$adapter_table <- renderTable(expr =  values$adapters,
                                                     rownames = TRUE)
                 update_adapter_index(session = session,
                                      adapters = values$adapters)
               })

  observeEvent(eventExpr = input$add_adapter,
               label = "Add adapter", {
                 print("Add an adapter")
                 values$adapters <- rbind(values$adapters, data.frame("Adapter"=input$add_adapter_sequence_input, "Description"=input$add_adapter_description_input))
                 update_adapter_index(session = session, adapters = values$adapters)
               })

  observeEvent(eventExpr = input$remove_adapter,
               label = "Remove adapter", {
                 print("Remove an adapter")
                 values$adapters <- values$adapters[-c(as.numeric(input$adapter_index)), ]
                 if(nrow(values$adapters)>0)
                   row.names(values$adapters) <- 1:nrow(values$adapters)
                 update_adapter_index(session = session, adapters = values$adapters)
               })

  observeEvent(eventExpr = input$save_adapters,
               label = "Save adapters", {
                 success <- save_info(x = values$adapters,
                                      path = paste(values$adapters_dir,
                                                   values$adapters_file,
                                                   sep = ""),
                                      delimeter = "#",
                                      col_names = FALSE)
                 if(isTRUE(success)){
                   output$adapter_changes_saved <- renderText(expr = {return("Saved")})
                 }
               })

  observeEvent(values$adapters, {output$adapter_changes_saved <- renderText({" "})})

  ### Common functions

  # Close the modal
  observeEvent(eventExpr = input$close_modal, {removeModal()})

  # Update progress
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- values$progress$getValue()
      value <- value + (values$progress$getMax() - value) / 5
    }
    values$progress$set(value = value, detail = detail)
  }

  ### Run the workflow on a dataset

  observeEvent(eventExpr = input$begin_processing,
               label = "Main workflow",{
    removeModal()
                 #### Toggle everything off!!
    # Create a Progress object
    values$progress <- shiny::Progress$new()
    values$progress$set(message = "Progress: ", value = 0)
    on.exit(values$progress$close())
    print(input$get_range)
    print(input$read_cutoff)

    # Get the current time
    timestamp <- strsplit(x = as.character(Sys.time()), split = " ")

    values$dataset_ID <- make_filename(filename = input$input_file,
                                       dataset_name = input$dataset_name)
    print(values$dataset_ID)
    print(values$output_dir)
    values$output_dir <- getAbsolutePath(create_output_dirs(out_dir = values$output_dir,
                                                            name = paste(values$dataset_ID,
                                                                         timestamp[[1]][1],
                                                                         gsub(pattern = ":",
                                                                              replacement = "-",
                                                                              x =timestamp[[1]][2]),
                                                                         sep = "_")))
    print(values$output_dir)
    trim_cmd <- make_trim_command(input_dir = getAbsolutePath(values$input_dir),
                                          output_dir = values$output_dir,
                                          dataset_ID = values$dataset_ID,
                                          input_file = input$input_file,
                                          adapter_file = paste(values$adapters_dir,
                                                               values$adapters_file,
                                                               sep = "/"),
                                          min_length = input$get_range[[1]])
    print(trim_cmd)
    print(input$cores)
    #trim_log <-
      run_trimmer(output_dir = values$output_dir,
                            dataset_ID = values$dataset_ID,
                            trim_cmd = trim_cmd)

    # write_data_to_TSV(data = trim_log,
    #                   path = values$output_dir,
    #                   filename = paste(values$dataset_ID,
    #                                    ".trim.log",
    #                                    sep = ""))

  })
}





### Listeners
#align_dataset_listener <- reactive({input$align_dataset})
#view_results_listener <- reactive({input$view_results})
#begin_processing_listener <- reactive({input$begin_processing})


#select_genome_listener <- reactive({input$select_genome})
#load_genome_listener <- reactive({input$load_genome})
#view_genome_listener <- reactive({input$view_genome})
#get_ensembl_genomes_listener <- reactive({input$get_ensembl_genomes})
#add_genome_listener <- reactive({input$add_genome})
#remove_genome_listener <- reactive({input$remove_genome})
#get_intervals_listener <- reactive({input$get_intervals})
#check_complete_genome_listener <- reactive({values$genomes})
#remove_gene_list_listener <- reactive({input$remove_gene_list})
#update_genome_view_listener <- reactive({values$update_genome_view})
#genome_fasta_finder_listener <- reactive({ input$genome_fasta_finder})
#genome_index_finder_listener <- reactive({ input$genome_index_finder})
#gene_list_finder_listener <- reactive({input$gene_list_finder})
#save_genomes_listener <- reactive({input$save_genomes})
#load_genome_listener <- reactive({input$load_genome})
#nrow_genomes_listener <- reactive({values$genomes})
#view_adapters_listener <- reactive({input$view_adapters})
#add_adapter_listener <- reactive({input$add_adapter})
#remove_adapter_listener <- reactive({input$remove_adapter})
#save_adapters_listener <- reactive({input$save_adapters})




#output$adapter_changes_saved <- renderText({" "})


# output$input_dir <- renderText({
#   # If a save directory is selected, display the name of the save directory
#   return(values$input_dir)
# })



# # Get the current time
# timestamp <- strsplit(x = as.character(Sys.time()), split = " ")
#
# # Make the full output directory with timestamp
# values$output_dir <- paste(values$save_dir, "/", timestamp[[1]][1], "_", timestamp[[1]][2], sep = "")
#


