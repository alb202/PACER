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
  values$genome_volumes <- c("Data"="../data/genomes","Home"="~/", "Root"="/")
  values$adapters <- read.delim(file = "../data/adapters/adapters.txt", stringsAsFactors = FALSE, header = FALSE, sep = "#", col.names = c("Adapter", "Description"))
  values$genomes <- read.delim(file = "../data/genomes/genomes.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t", row.names= NULL, col.names = c("Description", "Version", "Dataset", "Status", "Interval.Status", "Genome.FASTA", "Bowtie.Index", "Gene.Sets"))
#  values$genomes <- read.delim(file = "../data/genomes/genomes_new.txt", stringsAsFactors = FALSE, header = FALSE, sep = ",", row.names= NULL, col.names = c("Description", "Version", "Dataset", "Status", "Interval.Status", "Genome.FASTA", "Bowtie.Index", "Gene.Sets"))
  values$genome_dir <- "../data/genomes/"
  #print("values$genomes")
  #print(values$genomes)
  observe({
    # Set the initial value of the input directory to the first entry in the volumes list
    values$input_dir <- values$input_volumes[[1]]
    values$output_dir <- values$output_volumes[[1]]
    print("values$genomes")
    print(values$genomes)
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
  observeEvent(label = "Select genome modal", select_genome_listener(),priority = 100, {
    print("Select genome")
    showModal(genomes_modal)
    print("genome dir WBcel999")
    print(getwd())
    print(isolate(nrow(values$genomes)))
    #print(values$genome_dir)
    #print(paste(getwd(), values$genome_dir, sep = "/"))
    #print(get_ensembl_intervals(genome = "WBAAA", path = values$genome_dir))
    output$genome_table <- renderTable(width = "100%" , expr =  {
      #print(class(values$genomes))
      print("printing the genome A")
      print(isolate(nrow(values$genomes)))
      print(values$genomes)
      new_table <- values$genomes
      print(nrow(new_table))
      if(nrow(new_table)>0){
        print("Changing values of new_table to OK")
        new_table[,5:8][new_table[,5:8]!="None"] <- "OK"
        print(new_table)
        print("new_table")
        print(new_table)
      }
      return(new_table)
    }, rownames = TRUE, spacing = "s")
    # update_genome_index(session = session, genomes = values$genomes)
    # print("print(isolate(input$genome_index))")
    # print(isolate(input$genome_index))
    print("printing B")
    update_genome_index(session = session, genomes = values$genomes)
  })

  # new_table <- data.frame(lapply(X = new_table,
  #                                FUN = function(x){
  #                                  if(x[5] == "OK")
  #                                    x[4] <- "Ready"
  #                                  return(x)
  #                                }))

  observeEvent(eventExpr = view_genome_listener(), priority = 0,{
    print("view the genome details")
    print(isolate(nrow(values$genomes)))
    #toggleElement(id = "genome_details", condition = TRUE)
    toggleElement(id = "genome_details")
    # update_genome_index(session = session, genomes = values$genomes)
    print("print(isolate(input$genome_index))")
    print(isolate(input$genome_index))
  })

  observeEvent(eventExpr = c(values$genomes[input$genome_index, 5:8],
                             view_genome_listener()),
               label = "update view genome", priority = -5, {

                 print("reload the genome view")
                 if(isolate(nrow(values$genomes)>0)){
                   #update_genome_index(session = session, genomes =  values$genomes)
                   print("genome index")
                   print(isolate(input$genome_index))
                   genome_row <- values$genomes[input$genome_index, ]
                   print("genomes B")
                   print(isolate(nrow(values$genomes)))
                   print("genome row")
                   print(genome_row)
                   print("genome row 5:7")
                   print(genome_row[5:7])
                   print("genome row NA")
                   # print(sum(genome_row[5:7] == "None"))
                   # print((genome_row[5] == "None" &
                   #          genome_row[6] == "None" &
                   #          genome_row[7] == "None"))

                   print("get genome row again")
                   #genome_row <- values$genomes[input$genome_index, ]
                   print(genome_row)
                   output$genome_interval_status <- renderUI({
                     if(as.character(genome_row[5])=="None"){
                       return(HTML("<font color='red'>Incomplete</font>"))
                     } else {
                       return(HTML("<font color='green'>Ready</font>"))
                     }
                   })
                   print("update gene list status")
                   print(genome_row)
                   # print(strsplit(x = genome_row[["Gene.Sets"]],
                   #                split = ";",
                   #                fixed = TRUE)[1])
                   # print(length(strsplit(x = genome_row[["Gene.Sets"]],
                   #                       split = ";",
                   #                       fixed = TRUE)[1]))
                   # print(length(strsplit(x = genome_row[["Gene.Sets"]],
                   #                       split = ";",
                   #                       fixed = TRUE)[[1]]))
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

                   print("output$genome_fasta_location ")
                   output$genome_fasta_location <- renderUI({
                     if(as.character(genome_row[["Genome.FASTA"]])=="None"){
                       return(HTML("<font color='red'>Incomplete</font>"))
                     } else {
                       return(HTML("<font color='green'>",genome_row[["Genome.FASTA"]],"</font>"))
                     }
                   })

                   print(" output$genome_index_location")
                   output$genome_index_location <- renderUI({
                     if(as.character(genome_row[["Bowtie.Index"]])=="None"){
                       return(HTML("<font color='red'>Incomplete</font>"))
                     } else {
                       return(HTML("<font color='green'>",genome_row[["Bowtie.Index"]],"</font>"))
                     }
                   })

                   print("shiny file choosers")
                   shinyFileChoose(input = input,
                                   session = session,
                                   id = "genome_fasta_finder",
                                   roots = isolate(values$genome_volumes),
                                   filetypes = c("fasta", "fa"))
                   print("genome_index_finder dir choose")
                   shinyDirChoose(input = input,
                                  session = session,
                                  id = "genome_index_finder",
                                  roots = isolate(values$genome_volumes))
                   print("gene_list_finder dir choose")
                   shinyFileChoose(input = input,
                                   session = session,
                                   id = "gene_list_finder",
                                   roots = isolate(values$genome_volumes),
                                   filetypes = c("txt", "tsv"))
                 }
                 # update_genome_index(session = session, genomes = values$genomes)
                 # print("print(isolate(input$genome_index))")
                 print(isolate(input$genome_index))
               })

  observeEvent(eventExpr = genome_fasta_finder_listener(), {
    print("getting Genome.FASTA file")
    #print(input$genome_index)
    print("input$genome_fasta_finder")
    print(isolate(input$genome_fasta_finder))
    if(!is.null(input$genome_fasta_finder)){
      fasta_location <- parseFilePaths(roots = values$genome_volumes, selection = input$genome_fasta_finder)
      print("fasta_location")
      print(fasta_location)
      print(as.character(fasta_location[["datapath"]]))
      print("update values$genomes")
      values$genomes[isolate(input$genome_index), "Genome.FASTA"] <- as.character(fasta_location[["datapath"]])
      print(values$genomes)
      print(isolate(nrow(values$genomes)))
      # updateTextInput(session, "genome_fasta_finder", value = NA)
    }
  })

  observeEvent(eventExpr = genome_index_finder_listener(), {
    print("getting genome index directory")
    print(isolate(input$genome_index_finder))
    print(isolate(nrow(values$genomes)))
    if(!is.null(input$genome_index_finder)){
      index_location <- parseDirPath(roots = values$genome_volumes, selection = input$genome_index_finder)
      updateTextInput(session = session, inputId = "genome_index_finder", value = NULL)
      print(isolate(input$genome_index_finder))
      print("index_location")
      print(index_location)
      #print(as.character(index_location))
      values$genomes[isolate(input$genome_index), "Bowtie.Index"] <- getAbsolutePath(as.character(index_location))
      print("genomes D")
      print(values$genomes)
      # updateTextInput(session, "genome_fasta_finder", value = NA)
      print(isolate(nrow(values$genomes)))
    }
  })

  observeEvent(label = "Checking complete genome",
               priority = -10,
               eventExpr = check_complete_genome_listener()
               # , suspended = TRUE
               , {
    print("Checking if genome is complete")
    if(nrow(values$genomes)>0){
      print("checking for complete genomes")
      values$genomes[,"Status"] <- check_for_complete_genome(genome = values$genomes)
    }
  })

  observeEvent(eventExpr = gene_list_finder_listener(), {
    print("getting gene lists")
    #print(input$genome_index)
    print(isolate(nrow(values$genomes)))
    if(!is.null(input$gene_list_finder)){
      gene_list_location <- parseFilePaths(roots = values$genome_volumes, selection = input$gene_list_finder)
      print("fasta_location")
      print(gene_list_location)
      print(as.character(gene_list_location[["datapath"]]))
      isolate(
        if(values$genomes[input$genome_index, "Gene.Sets"]=="None"){
          values$genomes[[input$genome_index, "Gene.Sets"]] <- getAbsolutePath(gene_list_location[["datapath"]])
        } else {
          values$genomes[[input$genome_index, "Gene.Sets"]] <- paste(values$genomes[input$genome_index, "Gene.Sets"],
                                                                     getAbsolutePath(gene_list_location[["datapath"]]),
                                                                     sep = ";")
        }
      )
    }
  })

  observeEvent(label = "remove gene", eventExpr = remove_gene_list_listener(), {
    print("remove gene list")
    print(isolate(nrow(values$genomes)))
    print(isolate(input$gene_list_status))
    isolate(genome_row <- values$genomes[input$genome_index, ])
    isolate(new_list <- paste(strsplit(x = genome_row[["Gene.Sets"]],
                                       split = ";",
                                       fixed = TRUE)[[1]][-(as.integer(input$gene_list_status))], collapse = ";"))
    print("updating the new gene list")
    if(new_list=="") new_list <- "None"
    isolate(values$genomes[input$genome_index, "Gene.Sets"] <- new_list)
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
    print("update_ensembl_genome_index(session = session, genomes = values$ensembl_genome_index)")
    update_ensembl_genome_index(session = session, genomes = values$ensembl_genome_index)
  })

  observeEvent(add_genome_listener(), {
    print("Adding a new genome")
    print(input$ensembl_genome_index)
    if(values$ensembl_genome_index[input$ensembl_genome_index, "dataset"] %in% values$genomes$dataset) return(NULL)
    print(isolate(nrow(values$genomes)))
    values$genomes <- add_genome(genomes = values$genomes, new_genome = values$ensembl_genome_index[input$ensembl_genome_index, ])
    # update_genome_index(session = session, genomes = values$genomes)
    # if(nrow(values$genomes)==0){
    #   print("create new genome list")
    #   values$genomes <- data.frame(stringsAsFactors = FALSE, row.names = NULL, values$ensembl_genome_index[input$ensembl_genome_index, ], "Status"="Incomplete", "Interval.Status"="None", "Genome.FASTA"="None", "Bowtie.Index"="None", "Gene.Sets"="None")
    # } else {
    #   print("append genome to list")
    #   values$genomes <- rbind(values$genomes, data.frame(stringsAsFactors = FALSE, row.names = NULL, values$ensembl_genome_index[input$ensembl_genome_index, ], "Status"="Incomplete", "Interval.Status"="None", "Genome.FASTA"="None", "Bowtie.Index"="None", "Gene.Sets"="None"))
    # }
    print(values$genomes)
    print("update_genome_index EEE")
    update_genome_index(session = session, genomes = values$genomes)
    print(isolate(input$genome_index))
    #print()
    print(isolate(nrow(values$genomes)))
  })

  observeEvent(remove_genome_listener(), {
    print("removing a genome")
    if(!is.na(input$genome_index)){
      values$genomes <- values$genomes[-as.integer(input$genome_index), ]
      # if(nrow(values$genomes)>0)
      #   row.names(values$genomes) <- 1:nrow(values$genomes)
    }
    update_genome_index(session = session, genomes = values$genomes)
  })


  observeEvent(eventExpr = get_intervals_listener(), {
    print("Get gene intervals")
    genome <- values$genomes[isolate(input$genome_index), 3]
    print(genome)
    if(!check_for_intervals(path = "../data/genomes", genome = genome)){
      print("getting ensembl intervals")
      get_ensembl_intervals(path = "../data/genomes", genome = genome)
    }
    print("update status of ")
    values$genomes[isolate(input$genome_index), 5] <- "OK"
    print(isolate(nrow(values$genomes)))
  })

  observeEvent(save_genomes_listener(), {
    success <- save_info(x = values$genomes, path = "../data/genomes/genomes.txt", delimeter = "\t", col_names = TRUE)
    print(success)
    if(isTRUE(success)){
      output$genome_changes_saved <- renderText(expr = {return("Saved")})
    }
  })

  observeEvent(load_genome_listener(), {
    print("loading selected genome")
    selected_genome <- values$genomes[input$genome_index, ]
  if(selected_genome["Status"] == "Ready"){
    print("ready - A")
    values$selected_genome <- selected_genome
    print("ready - B")
    output$selected_genome <- renderText({return(selected_genome[["Version"]])})
  } else {
    print("Incomplete - A")
    values$selected_genome <- NULL
    print("Incomplete - B")
    output$selected_genome <- renderText({"Select a complete genome"})
  }

  })

  # observeEvent(nrow_genomes_listener(), {
  #   print("Number of rows has changed")
  #   if(nrow(values$genomes)>0){
  #     row.names(values$genomes) <- 1:nrow(values$genomes)
  #     }
  #   update_genome_index(session = session, genomes = values$genomes)
  # })

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
    success <- save_info(x = values$adapters, path = "../data/adapters/adapters_new.txt", delimeter = " #", col_names = FALSE)
    print(success)
    if(isTRUE(success)){
      output$adapter_changes_saved <- renderText(expr = {return("Saved")})
    }
  })

  ### Listeners
  align_dataset_listener <- reactive({input$align_dataset})
  view_results_listener <- reactive({input$view_results})
  begin_processing_listener <- reactive({input$begin_processing})

  select_genome_listener <- reactive({input$select_genome})
  load_genome_listener <- reactive({input$load_genome})
  view_genome_listener <- reactive({input$view_genome})
  #update_view_genome_listener <- reactive({values$update_genome_viewer})
  get_ensembl_genomes_listener <- reactive({input$get_ensembl_genomes})
  add_genome_listener <- reactive({input$add_genome})
  remove_genome_listener <- reactive({input$remove_genome})
  get_intervals_listener <- reactive({input$get_intervals})
  check_complete_genome_listener <- reactive({values$genomes})
  remove_gene_list_listener <- reactive({input$remove_gene_list})
  update_genome_view_listener <- reactive({values$update_genome_view})
  genome_fasta_finder_listener <- reactive({ input$genome_fasta_finder})
  genome_index_finder_listener <- reactive({ input$genome_index_finder})
  gene_list_finder_listener <- reactive({input$gene_list_finder})
  save_genomes_listener <- reactive({input$save_genomes})
  load_genome_listener <- reactive({input$load_genome})
  nrow_genomes_listener <- reactive({values$genomes})
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
