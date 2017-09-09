#ui.R
library(shiny)
library(shinyFiles)
library(shinyjs)
library(biomaRt)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "lumen.css"),
    tags$style(HTML(".modal-lg {width:1050px;padding:0;}"))),
    #tags$style(HTML(".modal-adapter {width:800px;padding:0;}"))),
  titlePanel("PACER"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id="sidebar", type = "tabs",
                  tabPanel(title = "Align", value = "align", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel(tags$hr(style = "margin:0;padding:0;"),
                             #column(10,
                             fluidRow(shinyFiles::shinyDirButton(id = "find_input_dir",
                                                                 label = "Select folder ...",
                                                                 title = "Select folder ...")),
                             fluidRow(h5(textOutput("input_dir"))),
                             #fluidRow(tags$hr()),
                             fluidRow(selectInput(inputId = "input_file", label = NULL, choices = NULL, selected = NULL,selectize = FALSE, multiple = FALSE, size = 5)),
                             fluidRow(
                               disabled(actionButton(inputId = "align_dataset",
                                                     icon = icon("check", lib = "glyphicon"),
                                                     label="Align dataset"))),
                             fluidRow(tags$hr(style = "margin:1;padding:0;")),
                             fluidRow(
                               actionButton(inputId = "select_genome",
                                            icon = icon("check", lib = "glyphicon"),
                                            label="Select genome")),
                             fluidRow(h5(textOutput("selected_genome"), style = "margin:1;padding:0;")),
                             fluidRow(tags$hr(style = "margin:1;padding:0;")),
                             fluidRow(
                               actionButton(inputId = "view_adapters",
                                            icon = icon("check",
                                                        lib = "glyphicon"),
                                            label="View Adapters"))
                           )
                  ),
                  tabPanel(title = "View", value = "view_results", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel(
                               fluidRow(shinyFiles::shinyDirButton(id = "find_output_dir",
                                                                   label = "Select folder ...",
                                                                   title = "Select folder ...")),
                               fluidRow(h5(textOutput("output_dir"))),
                               # actionButton(inputId = "view_results",
                               #                       icon = icon("check", lib = "glyphicon"),
                               #                       label="Load results")),
                             fluidRow(selectInput(inputId = "processed_file", label = NULL, choices = NULL, selected = NULL,selectize = FALSE, multiple = FALSE, size = 5)),
                             fluidRow(
                               disabled(actionButton(inputId = "view_dataset",
                                                     icon = icon("check", lib = "glyphicon"),
                                                     label="View dataset")))
                           )
                  ),
                  tabPanel(title = "Help", value = "help", icon = icon("help", lib="glyphicon")
                  )
      )),
    #hr(),
    #fluidRow(column(3, verbatimTextOutput("value")))))),

    mainPanel(
      tabsetPanel(id="main", type = "tabs",
                  tabPanel(title = "Align", value = "align", icon = icon("folder-open", lib="glyphicon"),
                           wellPanel(

                           ))))))

