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
  values$input_volumes <- c("Input"="../data/input","Home"="~/", "Root"="/")
  values$output_volumes <- c("Output"="../data/output","Home"="~/", "Root"="/")
  values$fasta_volumes <- c("Data"="../data/genomes","Home"="~/", "Root"="/")
  values$adapters <- read.delim(file = "../data/adapters/adapters.txt", header = FALSE, sep = "#", col.names = c("Adapter", "Description"))
  values$genomes <- read.delim(file = "../data/genomes/genomes.txt", header = FALSE, sep = ",", col.names = c("Description", "Version", "Dataset", "Status", "Interval Status", "Genome FASTA", "Bowtie Index", "Gene Sets"))
  values$genome_dir <- "../data/genomes/"

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
    print("genome dir WBcel999")
    print(getwd())
    #print(values$genome_dir)
    #print(paste(getwd(), values$genome_dir, sep = "/"))
    #print(get_ensembl_intervals(genome = "WBAAA", path = values$genome_dir))
    output$genome_table <- renderTable(expr =  {
      print(class(values$genomes))
      print(values$genomes)
      new_table <- values$genomes
      if(nrow(new_table)>0){
        new_table[,5:8][new_table[ ,5:8]!="NA"] <- "OK"
      }
      return(new_table)
      }, rownames = TRUE, spacing = "s")
    update_genome_index(session = session, genomes = values$genomes)
  })

  observeEvent(view_genome_listener(), {
    print("view the genome details")
    #toggleElement(id = "genome_details", condition = TRUE)
    toggleElement(id = "genome_details")
  })

  observeEvent(view_genome_listener(), {
    print("reload the genome view")
    genome_row <- values$genomes[input$genome_index, ]
    print(genome_row)
    output$genome_interval_status <- renderUI({
      if(as.character(genome_row[5])=="NA"){
        return(HTML("<font color='red'>Incomplete</font>"))
      } else {
        return(HTML("<font color='green'>Ready</font>"))
      }
    })

    output$gene_list_status <- renderUI({
      if(as.character(genome_row[["Gene.Sets"]])=="NA"){
        return(HTML("<font color='red'>None</font>"))
      } else {
        return(HTML(text = paste("<font color='green'>",
                    strsplit(x = genome_row[["Gene.Sets"]],
                             split = ";",
                             fixed = TRUE)[[1]],
                    "</font>", sep = "")))
      }
    })

    output$genome_fasta_location <- renderUI({
      if(as.character(genome_row[["Genome.FASTA"]])=="NA"){
        return(HTML("<font color='red'>Incomplete</font>"))
      } else {
        return(HTML("<font color='green'>",genome_row[["Genome.FASTA"]],"</font>"))
      }
    })

    shinyFileChoose(input = input,
                    session = session,
                    id = "genome_fasta_finder",
                    roots = isolate(values$fasta_volumes),
                    filetypes = c("fasta", "fa"))
    shinyFileChoose(input = input,
                    session = session,
                    id = "gene_list_finder",
                    roots = isolate(values$fasta_volumes),
                    filetypes = c("txt", "tsv"))

  })

  observe({
    print("getting genome fasta file")
    #print(input$genome_index)
    if(!is.null(input$genome_fasta_finder)){
      fasta_location <- parseFilePaths(roots = values$fasta_volumes, selection = input$genome_fasta_finder)
      print("fasta_location")
      print(fasta_location)
      print(as.character(fasta_location[["datapath"]]))
      values$genomes[isolate(input$genome_index), "Genome.FASTA"] <- as.character(fasta_location[["datapath"]])
      print(values$genomes)
      # updateTextInput(session, "genome_fasta_finder", value = NA)
    }
      #print(values$genomes[input$genome_index, "Genome.FASTA"])}
  })

  observe({
    print("getting gene lists")
    #print(input$genome_index)
    if(!is.null(input$gene_list_finder)){
      gene_list_location <- parseFilePaths(roots = values$fasta_volumes, selection = input$gene_list_finder)
      print("fasta_location")
      print(gene_list_location)
      print(as.character(gene_list_location[["datapath"]]))
      isolate(
        if(values$genomes[input$genome_index, "Gene.Sets"]=="NA"){
          values$genomes[[input$genome_index, "Gene.Sets"]] <- getAbsolutePath(gene_list_location[["datapath"]])
        } else {
          values$genomes[[input$genome_index, "Gene.Sets"]] <- paste(values$genomes[input$genome_index, "Gene.Sets"],
                                                                              getAbsolutePath(gene_list_location[["datapath"]]), sep = ";")
        }
      )
      #print(values$genomes)
      # updateTextInput(session, "genome_fasta_finder", value = NA)
    }
    #print(values$genomes[input$genome_index, "Genome.FASTA"])}
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
    values$genomes <- rbind(values$genomes, data.frame(stringsAsFactors = FALSE, row.names = NULL, values$ensembl_genome_index[input$ensembl_genome_index, ], "Status"="Incomplete", "Interval Status"="NA", "Genome FASTA"="NA", "Bowtie Index"="NA", "Gene Sets"="NA"))
    print(values$genomes)
    update_genome_index(session = session, genomes = values$genomes)
  })

  observeEvent(get_intervals_listener(), {
    print("Get gene intervals")
    genome <- values$genomes[isolate(input$genome_index), 3]
    print(genome)
    if(!check_for_intervals(path = "../data/genomes", genome = genome)){
      get_ensembl_intervals(path = "../data/genomes", genome = genome)
    }
    values$genomes[isolate(input$genome_index), 5] <- "OK"
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
  get_intervals_listener <- reactive({input$get_intervals})


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
