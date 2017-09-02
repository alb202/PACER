#server.R
library(shiny)
library(shinyFiles)
library(shinyjs)
library(DT)

source("shiny_helper.R")

server <- function(input, output, session){

  values <- reactiveValues()
  values$input_volumes <- c("Input"="../data/input","Home"="~/", "Root"="/")
  values$output_volumes <- c("Output"="../data/output","Home"="~/", "Root"="/")
  values$adapters <- read.delim(file = "../data/adapters/adapters.txt", header = FALSE, sep = "#", col.names = c("Adapter", "Description"))

  observe({
    # Set the initial value of the input directory to the first entry in the volumes list
    values$input_dir <- values$input_volumes[[1]]
    values$output_dir <- values$output_volumes[[1]]
  })

  observe( {
    # Only show the "Align dataset" button if a dataset has been selected
    toggleState(id = "align_dataset", condition = !is.null(input$input_file))
  })

  observe( {
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
    print("update select box")
    updateSelectInput(session = session,
                      inputId = "input_file",
                      label = NULL,
                      choices = load_directory(dir = values$input_dir, filetypes = c("fastq", "fq")),
                      selected = NULL)
  })

  observeEvent(values$output_dir, {
    # Update the file select box with the list of fastq files from the directory
    print("update select box")
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
    showModal(modalDialog(
      title = paste("File: ", input$input_file, sep = ""),
      easyClose = TRUE,
      textInput("dataset_name", "Choose a name for this dataset (optional)",
                placeholder = ""
      ),
      footer = tagList(
        modalButton(label = "Cancel"),
        actionButton(inputId = "begin_processing", label = "Begin...")
      )
    ))
  })

  observeEvent(view_adapters_listener(), {
    #value$adapters <- read.delim(file = "../data/adapters/adapters.txt", header = FALSE, sep = "#", col.names = c("Adapter", "Description"))
    showModal(modalDialog(size = "l",
      title = "Adapters",
      easyClose = TRUE,
      renderTable(expr =  values$adapters, rownames = TRUE),
      div(tags$hr()),
      fluidRow(
        column(width = 4, style = "margin-top:0px; margin:0px;",
               textInput(inputId = "add_adapter_sequence_input", label = "Adapter")),
        column(width = 4, style = "margin-top:0px; margin:0px;",
               textInput(inputId = "add_adapter_description_input", label = "Description")),
        column(width = 4, style = "margin-top:25px",
               actionButton(inputId = "add_adapter", label = "Add Adapter"))),
      div(tags$hr()),
      fluidRow(
        column(width = 4, style = "margin-top:0px",
               selectInput(inputId = "adapter_index",label = "Select adapter to remove", multiple = FALSE, choices = 1:nrow(values$adapters))),
        column(width = 2, style = "margin-top:24px",
               actionButton(inputId = "remove_adapter", label = "Remove Adapter"))
        #, column(width = 6)
        ),
      footer = tagList(
        modalButton(label = "Exit"),
        actionButton(inputId = "save_adapters", label = "Save Changes"),
        textOutput(outputId = "changes_saved")
      )
    ))
    output$changes_saved <- renderText(expr = {return("")})
  })

  observeEvent(add_adapter_listener(), {
    values$adapters <- rbind(values$adapters, data.frame("Adapter"=input$add_adapter_sequence_input, "Description"=input$add_adapter_description_input))
    if(nrow(values$adapters)>0){
      row.names(values$adapters) <- 1:nrow(values$adapters)
      row_choices <- 1:nrow(values$adapters)
      toggleState(id = "remove_adapter", condition = TRUE)
    } else{
      row_choices <- c("")
      toggleState(id = "remove_adapter", condition = TRUE==FALSE)
    }
    #choice_list <- ifelse(test = nrow(values$adapters)>0, yes = 1:nrow(values$adapters), no = "")
    print(values$adapters)
    updateSelectInput(session = session,
                      inputId = "adapter_index",
                      label = NULL,
                      choices = row_choices, #1:nrow(values$adapters),
                      selected = NULL)
  })
  observeEvent(remove_adapter_listener(), {
    values$adapters <- values$adapters[-c(as.numeric(input$adapter_index)), ]
    if(nrow(values$adapters)>0){
      row.names(values$adapters) <- 1:nrow(values$adapters)
      row_choices <- 1:nrow(values$adapters)
      toggleState(id = "remove_adapter", condition = TRUE)
    } else{
      row_choices <- c("")
      toggleState(id = "remove_adapter", condition = TRUE==FALSE)
    }
    print(input$adapter_index)
    print(values$adapters)
    #choice_list <- ifelse(test = nrow(values$adapters)>0, yes = 1:nrow(values$adapters), no = "")
    updateSelectInput(session = session,
                      inputId = "adapter_index",
                      label = NULL,
                      choices = row_choices, #1:nrow(values$adapters),
                      selected = NULL)
  })
  observeEvent(save_adapters_listener(), {
    success <- save_adapters(x = values$adapters, path = "../data/adapters/adapters_new.txt")
    print(success)
    if(isTRUE(success)){
      output$changes_saved <- renderText(expr = {return("Saved")})
    }
  })

  observeEvent(begin_processing_listener(), {



  })
  ### Listeners
  align_dataset_listener <- reactive({input$align_dataset})
  view_results_listener <- reactive({input$view_results})
  begin_processing_listener <- reactive({input$begin_processing})
  view_adapters_listener <- reactive({input$view_adapters})
  add_adapter_listener <- reactive({input$add_adapter})
  remove_adapter_listener <- reactive({input$remove_adapter})
  save_adapters_listener <- reactive({input$save_adapters})




   # output$input_dir <- renderText({
  #   # If a save directory is selected, display the name of the save directory
  #   return(values$input_dir)
  # })


}

# # Get the current time
# timestamp <- strsplit(x = as.character(Sys.time()), split = " ")
#
# # Make the full output directory with timestamp
# values$output_dir <- paste(values$save_dir, "/", timestamp[[1]][1], "_", timestamp[[1]][2], sep = "")

