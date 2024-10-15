library(shiny)
library(tidyverse)
library(DT)
library(shinyjs)
library(shinythemes)
library(reactable)
library(openxlsx)
library(writexl)
library(shinyalert)
library(renv)

source("logic/functions.R")

ui <- fluidPage(theme = shinytheme("flatly"),
  useShinyjs(),
  
  tags$head(
    tags$style(HTML(".shiny-output-error-validation {
                    color: red; }"))
  ),
  titlePanel("Fluoroprobe Analysis"),
  sidebarLayout(
    sidebarPanel(
      div(id = "File_Inputs",
          fileInput("fpdata_file", "Select FPdata file", accept = ".txt"),
          actionButton("deletefprows", "Delete FP Rows", class = "btn-info"),
          br(),
          br(),
          fileInput("labid_file", "Select LabID file", accept = ".csv"),
          actionButton("deleteIDrows", "Delete LabID Rows", class = "btn-info"),
          br(),
          br()
          ),
      textInput("description", "Description (optional)", value = ""),
      uiOutput("save_button"),
      span(textOutput("ready"), style = "color:green")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("FP Data",
          DT::dataTableOutput("fpdata_table")
        ),
        tabPanel("Lab ID",
          DT::dataTableOutput("labid_table")
        ),
        navbarMenu("QA/QC",
                   tabPanel("Duplicates",
                            reactableOutput("duplicates_table")),
                   tabPanel("Blanks",
                            reactableOutput("blanks_table"))
        ),
        tabPanel("FP Report",
                 reactableOutput("finalreport_table"))
      )
    )
  )
)

server <- function(input, output){

  #Create reactive value to store "raw" file contents
  origin <- reactiveValues(fpdata = NULL, labiddata = NULL)
  
  #Update ID reactive value upon file upload
  observeEvent(input$labid_file, {
    
    origin$labiddata <- read.csv(input$labid_file$datapath, header = FALSE, stringsAsFactors = FALSE, quote = "", skip = 1,
                           col.names = c("Number", "Site", "ID", "Type", "Start", "Svol", "MQvol", "Comments"),
                           colClasses = c("numeric", "factor", "character", "factor", "character", "numeric", "numeric", "character"))
    
  })
  
  #Update FP reactive value upon file upload
  observeEvent(input$fpdata_file, {
    req(input$fpdata_file)
    
    origin$fpdata <- read.table(input$fpdata_file$datapath, header = FALSE, sep = "\t", skip = 2, 
                    col.names = c("DateTime_Analyzed", "Greens", "Cyano", "Diatoms", "Crypto", 
                             "#5", "#6", "#7", "Yellow", "totChla", "Transmission", 
                             "Depth", "Temperature", "GreenCells", "CyanoCells", "DiatomCells", 
                             "CryptoCells", "#5cells", "#6cells", "#7cells", "Yellow2", 
                             "totCellCt", "T700", "LED3", "LED4", "LED5", "LED6", "LED7", 
                             "LED8", "Pressure", "TLED", "TSensor"),
                     colClasses = c("character", rep("numeric", 4), rep("NULL", 3), rep("numeric", 3), rep("NULL", 21)))
  })
  
  #Show raw FP data. User can delete errant rows (nrow() Always should be multiple of 10)
  output$fpdata_table <- renderDT({
    validate({
      req(input$fpdata_file)
      need(nrow(origin$fpdata)%%10 == 0, message = "FP data not a multiple of 10")
    })
    DT::datatable(origin$fpdata)
  })
  
  #Delete FP rows upon actionButton click
  observeEvent(input$deletefprows, {
    if(length(input$fpdata_table_rows_selected)==10){
      origin$fpdata <- origin$fpdata[as.numeric(-input$fpdata_table_rows_selected),]
    }
  })
  
  #Show raw ID data. Make "ID" column (3rd column) editable to fix names if needed
  output$labid_table <- renderDT({
    req(input$labid_file)
    
    DT::datatable(origin$labiddata,
                  options = list(
                    pageLength = -1,
                    lengthMenu = list(c(-1, 10), c("All", "10"))
                  ),
                  editable = list(target = "cell", disable = list(columns = c(1:2, 4:7))))
  })
  
  #Update index of corrected sample ID
  observeEvent(input$labid_table_cell_edit, {
    row  <- input$labid_table_cell_edit$row
    clmn <- input$labid_table_cell_edit$col
    origin$labIDdata[row, clmn] <- input$labid_table_cell_edit$value
  })
  
  #Delete ID rows upon actionButton click
  observeEvent(input$deleteIDrows, {
    if(!is.null(input$labid_table_rows_selected)){
      origin$labIDdata <- origin$labIDdata[as.numeric(-input$labid_table_rows_selected),]
    }
  })

  #Create fundamental dataframe for further calculations by merging FP and ID data
  merged_data <- reactive({
    merge_fp_id(origin$fpdata, origin$labiddata)
  })
  
  #Show sample Transmissions and RPD of duplicates
  output$duplicates_table <- renderReactable({
    validate({
      req(input$fpdata_file)
      req(input$labid_file)
      
      need(nrow(origin$fpdata)/nrow(origin$labiddata)== 10, message = "Data and LabID files not of proportional length")
    })
    process_duplicates(merged_data()) %>%
    reactable(
              pagination = F,
              columns = list(
                Transmission = colDef(style = function(value){
                  if (value < 90){
                    color <- "red"
                  }else {
                      color <- "#00CD00"
                    }
                    list(color = color)
                }),
                RPD_CV_Total = colDef(style = function(value){
                  if(is.na(value)==T){
                    color <- NULL
                  }
                  if(is.na(value) == F){
                    if(value > 15){
                    color <- "red"
                    } else {
                    color <- "#00CD00"
                    }
                  }
                  list(color = color)
                }),
                RPD_CV_Total_F = colDef(style = function(value){
                  if(is.na(value)==T){
                    color <- NULL
                  }
                  if(is.na(value) == F){
                    if(value > 15){
                      color <- "red"
                    } else {
                      color <- "#00CD00"
                    }
                  }
                  list(color = color)
                })
            )
    )
  })
  
  #Show blanks are below 0.2 ug/L total Chl-a
  output$blanks_table <- renderReactable({
    req(input$fpdata_file)
    req(input$labid_file)
    
    process_blanks(merged_data()) %>%
      reactable(pagination = F,
                columns = list(
                  totChla = colDef(style = function(value){
                    if (value > 0.2){
                      color <- "red"
                    } else{
                      color <- "#00CD00"
                    }
                    list(color = color)
                  }
                  )
                )
      )
  })
  
  #Show preview of processed output
  output$finalreport_table <- renderReactable({
    validate({
      req(input$fpdata_file)
      req(input$labid_file)
      need(nrow(origin$fpdata)/nrow(origin$labiddata)== 10, message = "Data and LabID files not of proportional length")
    })
    process_output_preview(merged_data()) %>%
    reactable(pagination = F)
  })

  #Display 'Save Data' button if FP and ID files match
  output$save_button <- renderUI({
    validate(
      need(str_sub(input$fpdata_file$name, 1, 6) == str_sub(input$labid_file$name, 1, 6), message = "Mistmatching Data and LabID files"),
      need(nrow(origin$fpdata)/nrow(origin$labiddata)== 10, message = "")
    )
      downloadButton("save_data", "Save Data", class = "btn-success")
  })
  
  #Save data upon clicking save button
  output$save_data <- downloadHandler(
    filename = function(){
      paste0(file_name(str_sub(input$fpdata_file$name, 1, 6), input$description), ".xlsx")
    },
    content = function(file){
      
      dataset <- list("Processed Data" = process_output_preview(merged_data()),
                      "QAQC" = process_duplicates(merged_data()) %>%
                        filter(!grepl("Blank", Type, ignore.case = T)),
                      "Summary of Blanks" = process_blanks(merged_data()), 
                      "Field Blanks" = process_field_blanks(merged_data()))
      dataset %>%
        write_xlsx(file)
      
      shinyalert(title = "Saved",
                 type = "success")
    })
  
  #Reset file params upon FP file uploads
  observeEvent(input$fpdata_file, {
    reset("File_params")
  })
  
  #Reset file params upon ID file uploads
  observeEvent(input$labid_file, {
    reset("File_params")
  })
}

shinyApp(ui = ui, server = server)