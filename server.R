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
library(shinycssloaders)
library(colourpicker)
library(shinyWidgets)
library(svglite)
library(sortable)
library(cowplot)
#library(grDevices)
library(tidyr)

source("Plotter.R")

#setwd("Y:/exchange/Marcus/M_pro/Microarray_plotter/")
l <-list.dirs("databases", full.names = FALSE, recursive = FALSE)


server <- function(input, output) {
  
  A <- reactive({
    if(input$own_anno){
      inFile <- input$file_anno$datapath
      if (is.null(inFile)) return(NULL)
    }else {
    lm <- list.files(str_c("databases/", input$spec), full.names = TRUE)
    inFile <- lm[which(str_detect(lm,"Plateform.csv"))]
    }
    read.csv(inFile, sep=";")
    })
  
  Gff <- reactive({
    if(input$own_map){
      inFile <- input$file_map$datapath
      if (is.null(inFile)) return(NULL)
    }else {
      lm <- list.files(str_c("databases/", input$spec), full.names = TRUE)
      inFile <- lm[which(str_detect(lm,"_map"%R%DOT))]
    }
      
      if(str_detect(inFile, ".csv"%R%END)){
        Gh <- load_csv(inFile)
      }else if(str_detect(inFile, ".gff3"%R%END)){
        Gh <- load_Gff3(inFile)
      }else if(str_detect(inFile, ".gff"%R%END)){
        Gh <- load_Gff(inFile)
      }
      
      Gh0 <- Gh
      Gh$start <- sapply(1:nrow(Gh), function(i)
        min(Gh0[i, which(colnames(Gh) %in%c("start", "end"))]))
      Gh$end <- sapply(1:nrow(Gh), function(i)
        max(Gh0[i, which(colnames(Gh) %in%c("start", "end"))]))
      Gh
  })
  
  
  Geo <- reactive({
    inFile <- input$file
    k <- str_extract_all(input$GEOname, one_or_more(WRD))[[1]]
    if (is.null(inFile) & length(k)==0) return(NULL)
    if(!is.null(inFile)){
      R <- read.delim(inFile$datapath)        #load geo-file
      R <- R[,which(colnames(R)!="X")]
    }
   if(length(k)!=0){
     for(i in k){
       g <- getGEO(i)
       
       if(str_sub(i,1,3)=="GSE"){
         Q <-g[[1]]@assayData$exprs
         colnames(Q) <-  g[[1]]@phenoData@data$title
       }else if(str_sub(i,1,3)=="GSM"){
         Q <-as.data.frame(g@dataTable@table$VALUE)
         colnames(Q) <-  g@header$title
       }
       if(i==k[1]){
         R1 <- Q
       }else{
         R1 <- cbind(R1,Q)
       }
     }
   }
    if(length(k)!=0 & !is.null(inFile)){
      R <- cbind(R, R1)
    }else if(length(k)!=0 & is.null(inFile)){
      R <- R1
    }
    
    if(!is.null(input$Name)){
      N <- str_remove(colnames(R), zero_or_more("[-_. ]")%R%one_or_more(DGT)%R%END)
      ex <- str_extract(colnames(R), zero_or_more("[-_. ]")%R%one_or_more(DGT)%R%END)
      Name <- readLines(input$Name$datapath)
      Name <- Name[which(Name!="" & Name!=" ")]
      
      N1 <- data.frame(old = unique(N[!is.na(N)]), new = Name)
      N <- sapply(N, function(i) N1$new[which(N1$old==i)])
      colnames(R) <- str_c(N,ex)
    }   
    
    R
  })
  
  output$whichReplicons <- renderUI({
    if(is.null(A())){ 
      p("please select a species")
    }else{
      U <- sort(unique(A()$Replicon[which(!is.na(A()$Replicon))]))
      tagList(selectInput("Rep", "Please select chromosome or plasmid:", 
                          choices =U, selected = "chromosome"))
    }
  })
  
  N <- reactive({unique(str_remove(colnames(Geo()), zero_or_more("[-_. ]")%R%one_or_more(DGT)%R%END))})
  
  output$whichmicros <- renderUI({
    if (is.null(N())){
      p("")
    } else{
      
      tagList(
        awesomeCheckboxGroup("filter_micro",  label = "Which Microarray should be displayed?", 
                             choices = N(),
                             status = "danger", selected = N()))
    }
  })
  
  output$ordermicros <- renderUI({
    if (is.null(N())){
      p("")
    } else{
      rank_list(
        text = "Choose the order of Microarrays in the plot:",
        labels = N(),
        input_id = "order_micro")
    }
  })
  
  
  output$whichcolor <- renderUI({
    if (is.null(N())){
      p("")
    } else{
      colo <- col1[1:length(N())]    #create colors for RNA-reads
      
      tagList(
        lapply(1:length(N()), function(i) {
          colourInput(str_c("P",i), label=N()[i], showColour = "both", value = colo[i])}))
    }
  })
  
  output$whichplot <- renderUI({
    withSpinner(plotOutput("plot", width = str_c(input$width,"in"), height = str_c(input$height,"in")), type=6, hide.ui = FALSE)
  })
  
  output$plot <- renderPlot(NULL)
  
  col <- reactive({
    if (is.null(input$P1)){
      co <- col1[1:length(N())]
    } else{
      co <- c(input$P1, input$P2, input$P3, input$P4, input$P5, input$P6, input$P7, input$P8,
              input$P9, input$P10, input$P11, input$P12, input$P13, input$P14, input$P15, input$P16,
              input$P17, input$P18, input$P19, input$P20, input$P21, input$P22, input$P23, input$P24)[1:length(N())]
    }
    names(co) <- N()
    co
  })
  
  observeEvent(input$do, {
    output$plot <- renderPlot({
      
      isolate({ Microarrayplot(Geo(), A(),Gff(), input$start, input$end, 
                               replicon = input$Rep,
                               #alpha = input$alpha,
                               graph_size = input$graph_size,
                               subgenes = isolate(input$subgenes),
                               line_visible = input$line_visible,
                               incom_genes = input$incom_genes,
                               color = col(),
                               min_micro = input$min_micro,
                               max_micro = input$max_micro,
                               arrow_body_height = input$arrow_body_height,
                               arrowhead_height = input$arrowhead_height,
                               arrowhead_width = input$arrowhead_width,
                               filter_micro = input$filter_micro, 
                               order_micro = input$order_micro, 
                               Gsize= input$Gsize,
                               msize= input$msize,
                               Mwidth = input$width,
                               ntlength = input$ntlength,
                               lineSize = input$Gsize) })
    })
    
  })
  
  
  output$map_table <- DT::renderDT({
    mapD2 <- Gff() %>% select(Replicon, gene, start, end, orientation)
    DT::datatable(mapD2)#, editable = TRUE)
  })
  
  output$microarray_table <- DT::renderDT({
    mapD2 <- A()  %>% select(Row, Col, ControlType, Name, Replicon, start, end, orientation)
    DT::datatable(mapD2)#, editable = TRUE)
  })
  
  output$foo = downloadHandler(
    filename = function() {
      paste("Microarrayplot.", input$download_type, sep="")
    },
    content = function(file) {
      data_to_save <- Microarrayplot(Geo(), A(),Gff(), input$start, input$end, 
                          replicon = input$Rep,
                          #alpha = input$alpha,
                          graph_size = input$graph_size,
                          subgenes = isolate(input$subgenes),
                          line_visible = input$line_visible,
                          incom_genes = input$incom_genes,
                          color = col(),
                          min_micro = input$min_micro,
                          max_micro = input$max_micro,
                          arrow_body_height = input$arrow_body_height,
                          arrowhead_height = input$arrowhead_height,
                          arrowhead_width = input$arrowhead_width,
                          filter_micro = input$filter_micro, 
                          order_micro = input$order_micro, 
                          Gsize= input$Gsize,
                          msize= input$msize,
                          Mwidth = input$width,
                          ntlength = input$ntlength,
                          lineSize = input$Gsize)
      if(input$download_type=="svg"){
        save_plot(file, plot = data_to_save, base_width = input$width*1.25, base_height = input$height*1.25, units = "in", 
                  device = input$download_type, fix_text_size = FALSE)
      }else{
        save_plot(file, plot = data_to_save, bg = "white", base_width = input$width*1.25, base_height = input$height*1.25, units = "in", 
                  device = input$download_type)
      }
    })
  output$tabmap <- downloadHandler(
    filename = function(){"Map.csv"}, 
    content = function(fname){
      write.csv2(Gff(), fname, row.names = FALSE )
    }
  )
  
  output$tabgo <- downloadHandler(
    filename = function(){"MicroarrayPlatform.csv"}, 
    content = function(fname){
      write.csv2(A(), fname, row.names = FALSE )
    }
  )
  
  
  output$Tutorial <- downloadHandler(
    filename = "MicroarrayPlotter_manual.pdf",
    content = function(file) {
      file.copy("MicroarrayPlotter_manual.pdf", file)
    }
  )
  
}

