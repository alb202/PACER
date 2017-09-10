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
  titlePanel("PACER"),
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
                                                   width = "100%")
                                     ), tags$br(), tags$br(),
                                     fluidRow(
                                       numericInput(inputId = "read_cutoff",
                                                   label = "Cutoff for over-represented reads",
                                                   value = .001,
                                                   max = 1,
                                                   width = "100%")
                                     )

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
                  tabPanel(title = "Align", value = "align", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel())))))

