align_modal <- modalDialog(
  size = "m",
  title = textOutput(outputId = "align_modal_title"),
  easyClose = TRUE,
  fluidRow(column(width = 12,
                  textInput("dataset_name", "Choose a name for this dataset (optional)",
                            placeholder = ""
                  ))),
  footer = tagList(
    modalButton(label = "Cancel"),
    actionButton(inputId = "begin_processing", label = "Begin...")
  )
)

adapter_modal <- modalDialog(size = "l",
                             title = "Adapters",
                             easyClose = FALSE,
                             fluidRow(column(width = 12, tableOutput("adapter_table"))),
                             # renderTable(expr =  values$adapters,
                             #             rownames = TRUE),
                             div(tags$hr()),
                             fluidRow(
                               column(width = 4, style = "margin-top:0px; margin:0px;",
                                      textInput(inputId = "add_adapter_sequence_input",
                                                label = "Adapter")),
                               column(width = 4, style = "margin-top:0px; margin:0px;",
                                      textInput(inputId = "add_adapter_description_input",
                                                label = "Description")),
                               column(width = 4, style = "margin-top:25px",
                                      actionButton(inputId = "add_adapter",
                                                   label = "Add Adapter"))),
                             div(tags$hr()),
                             fluidRow(
                               column(width = 6, style = "margin-top:0px",
                                      selectInput(inputId = "adapter_index",
                                                  label = "Select adapter to remove",
                                                  multiple = FALSE,
                                                  choices = NULL)),
                               column(width = 2, style = "margin-top:24px",
                                      actionButton(inputId = "remove_adapter",
                                                   label = "Remove Adapter"))
                             ),
                             footer = tagList(
                               modalButton(label = "Exit"),
                               actionButton(inputId = "save_adapters",
                                            label = "Save Changes"),
                               textOutput(outputId = "changes_saved")
                             )
)

genomes_modal <- modalDialog(size = "l",
                             title = "Genomes",
                             easyClose = FALSE,
                             fluidRow(column(width = 12, tableOutput(outputId = "genome_table"))),
                               # renderTable(expr =  values$genomes,
                               #             rownames = TRUE)),
                             fluidRow(
                               column(width = 4,
                                      style = "margin-top:0px;",
                                      selectInput(inputId = "genome_index",
                                                  label = "Select genome",
                                                  multiple = FALSE,
                                                  choices = "")),
                               column(width = 2, style = "margin-top:23px;",
                                      disabled(actionButton(inputId = "load_genome",
                                                            label = "Load Genome"))),
                               column(width = 2, style = "margin-top:23px;",
                                      disabled(actionButton(inputId = "view_genome",
                                                            label = "View Genome Info")))),
                             div(tags$hr()),
                             hidden(
                               div(id = "genome_details", style = "border:5;border-color:grey;",
                                   wellPanel(
                                     fluidRow(
                                       column(width = 2, style = "margin-top:25px",
                                              actionButton(inputId = "get_intervals", label = "Get Intervals")),
                                       column(width = 2, style = "margin-top:0px; margin:0px;",
                                              "ENSEMBL Intervals"),
                                       column(width = 2, style = "margin-top:0px; margin:0px;",
                                              htmlOutput(outputId = "genome_interval_status"))
                                     ),
                                     fluidRow(
                                       column(width = 2, style = "margin-top:25px",
                                              shinyFilesButton(id = "gene_list_finder",
                                                               title = "Find gene lists",
                                                               label = "Find gene list(s) ...",
                                                               multiple = TRUE)),
                                       column(width = 2, style = "margin-top:0px; margin:0px;",
                                              "Gene lists"),
                                       column(width = 2, style = "margin-top:0px; margin:0px;",
                                              htmlOutput(outputId = "gene_list_status"))
                                     ),
                                     fluidRow(
                                       column(width = 2, style = "margin-top:25px",
                                              shinyFilesButton(id = "genome_fasta_finder",
                                                               title = "Find genome FASTA",
                                                               label = "Find file ...",
                                                               multiple = FALSE)),
                                       column(width = 2, style = "margin-top:0px; margin:0px;",
                                              "ENSEMBL Genome FASTA"),
                                       column(width = 2, style = "margin-top:0px; margin:0px;",
                                              htmlOutput(outputId = "genome_fasta_location"))
                                     )
                                   )),
                               div(tags$hr())),
                             fluidRow(
                               column(width = 6, style = "margin-top:0px",
                                      selectInput(inputId = "ensembl_genome_index", selectize=TRUE, width = "100%",
                                                  label = "ENSEMBL Genomes", multiple = FALSE, choices = "")),
                               column(width = 3, style = "margin-top:0px",
                                      actionButton(inputId = "get_ensembl_genomes", label = "View ENSEMBL genomes")),
                               column(width = 3, style = "margin-top:24px",
                                      disabled(actionButton(inputId = "add_genome", label = "Add genome")))),
                             footer = tagList(
                               modalButton(label = "Exit"),
                               actionButton(inputId = "save_adapters", label = "Save Changes"),
                               textOutput(outputId = "changes_saved")
                             )
)
