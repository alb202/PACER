#ui.R
library(shiny)
library(shinyFiles)
library(shinyjs)
library(biomaRt)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "lumen.css"),
    tags$style(HTML(".modal-lg {width:1100px;padding:0;}"))),
  fluidRow(h6('      PACER 2.0')),
  #titlePanel("PACER 2.0"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id="sidebar", type = "tabs",
                  tabPanel(title = "Align", value = "align", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel(tags$hr(style = "margin:0;padding:0;"),
                                     fluidRow(shinyFiles::shinyDirButton(id = "find_input_dir",
                                                                         label = "Select folder ...",
                                                                         title = "Select folder ...")),
                                     fluidRow(h5(textOutput("input_dir"))),
                                     fluidRow(selectInput(inputId = "input_file",
                                                          label = NULL,
                                                          choices = NULL,
                                                          selected = NULL,
                                                          selectize = FALSE,
                                                          multiple = FALSE,
                                                          size = 5)),
                                     fluidRow(
                                       disabled(actionButton(inputId = "use_dataset",
                                                             icon = icon("play",
                                                                         lib = "glyphicon"),
                                                             label="Use this dataset ..."))),
                                     fluidRow(tags$hr(style = "margin:1;padding:0;")),
                                     fluidRow(
                                       actionButton(inputId = "select_genome",
                                                    #icon = icon("check", lib = "glyphicon"),
                                                    label="Select genome")),
                                     fluidRow(h5(textOutput("selected_genome"), style = "margin:1;padding:0;")),
                                     fluidRow(tags$hr(style = "margin:1;padding:0;")))),

                  tabPanel(title = "Options", value = "options", icon = icon("cog", lib="glyphicon"),
                           wellPanel(
                             fluidRow(
                               actionButton(inputId = "view_adapters",
                                            label="View Adapters")
                             ), tags$br(), tags$br(),
                             fluidRow(
                               sliderInput(inputId = "get_range",
                                           label = "Range of read lengths",
                                           min = 10,
                                           max = 30,
                                           value = c(12, 30),
                                           step = 1,
                                           ticks = TRUE,
                                           width = "90%")
                             ), tags$br(), tags$br(),
                             fluidRow(
                               numericInput(inputId = "read_cutoff",
                                            label = "Cutoff for over-represented reads",
                                            value = .001,
                                            max = 1,
                                            width = "70%")),tags$br(), tags$br(),
                             fluidRow(
                               numericInput(inputId = "cores",
                                            label = "Processing cores",
                                            value = 2,
                                            min = 1,
                                            max = 64,
                                            step = 1,
                                            width = "50%"))
                           )),
                  tabPanel(title = "View", value = "view_results", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel(
                             fluidRow(shinyFiles::shinyDirButton(id = "find_output_dir",
                                                                 label = "Select folder ...",
                                                                 title = "Select folder ...")),
                             fluidRow(h5(textOutput("output_dir"))),
                             fluidRow(selectInput(inputId = "processed_file",
                                                  label = NULL,
                                                  choices = NULL,
                                                  selected = NULL,
                                                  selectize = FALSE,
                                                  multiple = FALSE,
                                                  size = 5)),
                             fluidRow(
                               disabled(actionButton(inputId = "view_dataset",
                                                     icon = icon("check", lib = "glyphicon"),
                                                     label="View dataset")))
                           )
                  ),
                  tabPanel(title = "Help", value = "help", icon = icon("help", lib="glyphicon")
                  )
      )),
    mainPanel(
      tabsetPanel(id="main", type = "tabs",
                  tabPanel(title = "Statistics", value = "stats", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel(
                             fluidRow(h5(HTML("Trimming log"))),
                             fluidRow(uiOutput("trim_log"))
                           ),
                           wellPanel(
                             fluidRow(h5(HTML("Alignment log"))),
                             fluidRow(uiOutput('align_log'))
                           )
                  ),
                  tabPanel(title = "Dataset1", value = "dataset1", icon = icon("folder-open", lib="glyphicon"),
                           tabsetPanel(id="datasets", type = "tabs",
                                       tabPanel(title = "5' Plots", value = "five_prime", icon = icon("folder-open", lib="glyphicon"),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - sense and antisense reads"))),
                                                  fluidRow(plotOutput("five_prime_plot__two_mm__both"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - sense reads"))),
                                                  fluidRow(plotOutput("five_prime_plot__two_mm__sense"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads"))),
                                                  fluidRow(plotOutput("five_prime_plot__two_mm__antisense"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - sense and antisense reads - 5' assigned to 22nt"))),
                                                  fluidRow(plotOutput("five_prime_plot__two_mm__both__22nt_5prime"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - sense reads - 5' assigned to 22nt"))),
                                                  fluidRow(plotOutput("five_prime_plot__two_mm__sense__22nt_5prime"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt"))),
                                                  fluidRow(plotOutput("five_prime_plot__two_mm__antisense__22nt_5prime"))
                                                )),
                                       tabPanel(title = "Scatters", value = "scatters", icon = icon("folder-open", lib="glyphicon"),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 15nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_15"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 16nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_16"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 17nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_17"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 18nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_18"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 19nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_19"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 20nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_20"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 21nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_21"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 23nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_23"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 24nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_24"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 25nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_25"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 26nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_26"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 27nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_27"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 28nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_28"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 29nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_29"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 30nt x 22nt"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22_30"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 22G x 22A"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22G_22A"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 22G x 22C"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22G_22C"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("Two mismatch - antisense reads - 5' assigned to 22nt - 22G x 22T"))),
                                                  fluidRow(plotOutput("scatter__two_mm__antisense__22nt_5prime__22G_22T"))
                                                )),
                                       tabPanel(title = "Logos", value = "logos", icon = icon("folder-open", lib="glyphicon"),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 20nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__20nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 21nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__21nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 22nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__22nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 23nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__23nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 24nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__24nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 25nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__25nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch in seed - antisense reads - 26nt sequence vs shuffled"))),
                                                  fluidRow(plotOutput("seq_logo__no_mm_in_seed__antisense__26nt"))
                                                )
                                       ),
                                       tabPanel(title = "Phasing", value = "phasing", icon = icon("folder-open", lib="glyphicon"),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch - 22nt"))),
                                                  fluidRow(plotOutput("phasing__no_mm__22nt"))
                                                ),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch - 26nt"))),
                                                  fluidRow(plotOutput("phasing__no_mm__26nt"))
                                                )
                                        ),
                                       tabPanel(title = "Offsets", value = "offsets", icon = icon("folder-open", lib="glyphicon"),
                                                wellPanel(
                                                  fluidRow(h5(HTML("No mismatch - Offset from 22nt"))),
                                                  fluidRow(plotOutput("offsets__no_mm__5prime_22nt__sense")),
                                                  fluidRow(plotOutput("offsets__no_mm__5prime_22nt__antisense"))
                                                )
                                       ),
                                       tabPanel(title = "Heatmaps", value = "heatmaps", icon = icon("folder-open", lib="glyphicon"))
                           )
                  )
      )
    )
  )
)
