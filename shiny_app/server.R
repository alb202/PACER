#server.R
library(shiny)
library(shinyFiles)
library(shinyjs)
library(DT)

source("shiny_helper.R")
source("modals.R")

server <- function(input, output, session){

  values <- reactiveValues()
  values$input_volumes <- c("Input"="../data/input","Home"="~/", "Root"="/")
  values$output_volumes <- c("Output"="../data/output","Home"="~/", "Root"="/")
  values$adapters <- read.delim(file = "../data/adapters/adapters.txt", header = FALSE, sep = "#", col.names = c("Adapter", "Description"))
  values$genomes <- read.delim(file = "../data/genomes/genomes.txt", header = FALSE, sep = ",", col.names = c("Description", "Version", "Dataset", "Status", "Interval Status", "Genome FASTA", "Bowtie Index", "Gene Sets"))

  observe({
    # Set the initial value of the input directory to the first entry in the volumes list
    values$input_dir <- values$input_volumes[[1]]
    values$output_dir <- values$output_volumes[[1]]
  })

  observe({
    # Only show the "Align dataset" button if a dataset has been selected
    toggleState(id = "align_dataset", condition = !is.null(input$input_file))
  })

  observe({
    # Only show the "Align dataset" button if a dataset has been selected
    toggleState(id = "view_dataset", condition = !is.null(input$processed_file))
  })

  output$input_dir <- renderText({
    # Print the dataset directory
    print("render directory name")
    return(values$input_dir)
  })

  output$output_dir <- renderText({
    # Print the dataset directory
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
      values$input_dir <- parseDirPath(values$volumes, input$find_input_dir)
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
      values$output_dir <- parseDirPath(values$output_volumes, input$find_output_dir)
    # Print the output directories
    print(values$output_dir)
  })

  observeEvent(values$input_dir, {
    # Update the file select box with the list of fastq files from the directory
    print("update input box")
    updateSelectInput(session = session,
                      inputId = "input_file",
                      label = NULL,
                      choices = load_directory(dir = values$input_dir, filetypes = c("fastq", "fq")),
                      selected = NULL)
  })

  observeEvent(values$output_dir, {
    # Update the file select box with the list of fastq files from the directory
    print("update output box")
    updateSelectInput(session = session,
                      inputId = "processed_file",
                      label = NULL,
                      choices = load_directory(dir = values$output_dir, filetypes = "", include_dirs = TRUE),
                      selected = NULL)
  })

  observeEvent(align_dataset_listener(), {
    print("button pressed")
    input_file <- input$input_file
    print(input_file)
    showModal(align_modal)
    output$align_modal_title <- renderText({paste("File: ", input$input_file, sep = "")})
  })

  ### Genome dialogue
  observeEvent(select_genome_listener(), {
    print("Select genome")
    showModal(genomes_modal)
    output$genome_table <- renderTable(expr =  values$genomes, rownames = TRUE)
    update_genome_index(session = session, genomes = values$genomes)
  })

  observeEvent(view_genome_listener(), {
    print("view genome listener")
    toggleElement(id = "genome_details", condition = TRUE)
    genome_row <- values$genomes[input$genome_index, ]
    print(genome_row)
    output$genome_interval_status <- renderUI({
      if(as.character(genome_row[1,"Interval.Status"])=="NA"){
        return(HTML("<font color='red'>Incomplete</font>"))
      } else {
        return(HTML("<font color='green'>Ready</font>"))
      }
    })
  })

  observeEvent(get_ensembl_genomes_listener(), {
    print("get_ensembl_genomes_listener")
    withProgress(message = "Getting genomes from ENSEMBL", {
      values$mart_info <- listMarts()[1,]
      incProgress(amount = .4)
      values$biomart <- useMart(biomart = values$mart_info[[1,1]])
      incProgress(amount = .8)
      values$ensembl_genome_index <- listDatasets(mart = values$biomart)
    })
    update_ensembl_genome_index(session = session, genomes = values$ensembl_genome_index)
  })

  observeEvent(add_genome_listener(), {
    print(input$ensembl_genome_index)
    if(values$ensembl_genome_index[input$ensembl_genome_index, "dataset"] %in% values$genomes$dataset) return(NULL)
    values$genomes <- rbind(values$genomes, data.frame(row.names = NULL, values$ensembl_genome_index[input$ensembl_genome_index, ], "Status"="Incomplete", "Interval Status"="NA", "Genome FASTA"="NA", "Bowtie Index"="NA", "Gene Sets"="NA"))
    print(values$genomes)
    update_genome_index(session = session, genomes = values$genomes)
  })

  ### Adapter modal
  observeEvent(view_adapters_listener(), {
    print("Viewing adapters modal")
    showModal(adapter_modal)
    output$adapter_table <- renderTable(expr =  values$adapters, rownames = TRUE)
    update_adapter_index(session = session, adapters = values$adapters)
    output$changes_saved <- renderText(expr = {return("")})
  })

  observeEvent(add_adapter_listener(), {
    print("Adding an adapter")
    values$adapters <- rbind(values$adapters, data.frame("Adapter"=input$add_adapter_sequence_input, "Description"=input$add_adapter_description_input))
    update_adapter_index(session = session, adapters = values$adapters)
  })

  observeEvent(remove_adapter_listener(), {
    values$adapters <- values$adapters[-c(as.numeric(input$adapter_index)), ]
    update_adapter_index(session = session, adapters = values$adapters)
    print("Removing an adapter")
  })

  observeEvent(save_adapters_listener(), {
    success <- save_adapters(x = values$adapters, path = "../data/adapters/adapters_new.txt")
    print(success)
    if(isTRUE(success)){
      output$changes_saved <- renderText(expr = {return("Saved")})
    }
  })

  ### Listeners
  align_dataset_listener <- reactive({input$align_dataset})
  view_results_listener <- reactive({input$view_results})
  begin_processing_listener <- reactive({input$begin_processing})

  select_genome_listener <- reactive({input$select_genome})
  load_genome_listener <- reactive({input$load_genome})
  view_genome_listener <- reactive({input$view_genome})
  get_ensembl_genomes_listener <- reactive({input$get_ensembl_genomes})
  add_genome_listener <- reactive({input$add_genome})

  view_adapters_listener <- reactive({input$view_adapters})
  add_adapter_listener <- reactive({input$add_adapter})
  remove_adapter_listener <- reactive({input$remove_adapter})
  save_adapters_listener <- reactive({input$save_adapters})




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

}
