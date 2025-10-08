library(stringr)
library(ggplot2)
library(gggenes)
library(aplot)
library(ggeasy)
library(RColorBrewer)
library(labeling) #eventually not necessary
library(rebus)    #eventually not necessary
library(stringr)
library(shiny)
library(DT)
library(dplyr)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(colourpicker)
library(GEOquery)

source("Plotter.R")

###shiny
options(shiny.maxRequestSize=200*1024^2)

ui <-fluidPage(
  titlePanel("Microarray-Plotter"),
  theme = shinythemes::shinytheme("spacelab"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Main",
                 selectInput("spec", "Please select organism:", choices =l),
                 prettySwitch("own_map", label = "Do you want to use own Gene-map?", value= FALSE),
                 prettySwitch("own_anno", label = "Do you want to use own Microarray-Platform?", value= FALSE),
                 conditionalPanel(condition = "input.own_map",
                                  fileInput("file_map", "enter map-file [csv/gff/gff3]", multiple = FALSE,
                                            accept = c(".csv", ".gff", ".gff3"))),
                 conditionalPanel(condition = "input.own_anno",
                                  fileInput("file_anno", "enter Platform-file", multiple = FALSE,
                                            accept = c(".csv"))),
                 #prettySwitch("whichdata", label = "Do you want to add a fiel or a Geo-ID?", value= FALSE),
                 fileInput("file", "Choose Microarray GEO-file", multiple = FALSE,
                           accept = ".txt"),
                 textInput("GEOname", "Do you want to add Geo-data?"),
                 #prettySwitch("add_reads", label = "Do you want to add RNAread-data?", value= FALSE),
                 fileInput("Name", "enter names of Microarray conditions", multiple = FALSE, accept = c(".txt")),
                 uiOutput("whichReplicons"),
                 numericInput("start", "Please enter the start-position", 1, step = 1000),
                 numericInput("end", "Please enter the end-position", 1e4, step = 1000),
                 actionButton("do", "Plot"),
                 downloadButton('foo'),
                 sliderInput("width", "Width:", min = 3, max=30, value =11),
                 sliderInput("height", "Height:", min = 3, max=30, value =7),
                 selectInput("download_type", "Output Format:", 
                             choices=c("png", "pdf", "tex", "jpeg", "tiff", "bmp", "svg"),
                             selected = "png", width="40%")),
        tabPanel("optional",
                 numericInput("Gsize", "Graphsize:", 1.2, step = 0.1),
                 numericInput("min_micro", "min. y-axis:", NA, step = 0.1),
                 numericInput("max_micro", "max. y-axis:", NA, step = 0.1),
                 prettySwitch("subgenes", label = "display subgenes", value= TRUE),
                 prettySwitch("label", label = "change label size?", value= FALSE),
                 conditionalPanel(condition = "input.label",
                                  numericInput("msize", "label fontsize:", value=4)),
                 conditionalPanel(condition = "input.label",
                                  numericInput("ntlength", "only label genes longer than [nt]:", value=NULL)),
                 prettySwitch("line_visible", label = "line in maps?", value= TRUE),
                 uiOutput("whichmicros"),
                 uiOutput("ordermicros"),
                 textInput("incom_genes", "ending of incomplete genes:", 
                           value = "_partial", width = NULL, placeholder = NULL),
                 numericInput("graph_size", "map-graph ratio:", 3),
                 numericInput("arrow_body_height", "height of arrowbody:", 7),
                 numericInput("arrowhead_height", "height of arrow:", 10),
                 numericInput("arrowhead_width", "width of arrow:", 8)),
        tabPanel("color",
                 uiOutput("whichcolor"))
        
      )),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", uiOutput("whichplot")),
        tabPanel("Gene_map", DT::DTOutput('map_table'),
                 downloadButton('tabmap', label = "Download table")),
        tabPanel("Microarray_table", DT::DTOutput('microarray_table'),
                 downloadButton('tabgo', label = "Download table"))
      )
    )
  ),
  absolutePanel(
    top = 10, right = 10, draggable = TRUE,
    downloadButton("Tutorial", "Manual")
  )
)

