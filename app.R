library(shiny)
library(tidyverse)
library(DT)
library(shinyjs)
library(shinythemes)
library(reactable)

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
          actionButton("deletefprows", "Delete FP Rows", class = "btn-danger"),
          br(),
          br(),
          fileInput("labid_file", "Select LabID file", accept = ".csv"),
          actionButton("deleteIDrows", "Delete LabID Rows", class = "btn-danger"),
          br(),
          br()
          ),
      
      div(id = "File_params",
          dateInput("date_sel", "Select Date", value = NULL, format = "yy/mm/dd"),
          checkboxGroupInput("project_sel", "Select Relevant Project/Lake",
                             choices = c("Owasco" = "OW",
                                         "Honeoye" = "HN",
                                         "Seneca"= "SC",
                                         "KLA",
                                         "MC")),
          textInput("comment", "Comment", placeholder = "Include any descriptive text here")
          ),
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
        tabPanel("Raw Results",
                 reactableOutput("rawresults_table")
        ),
        tabPanel("FP Report",
                 DT::dataTableOutput("finalreport_table"))
      )
    )
  )
)

server <- function(input, output){

#-------------------Read FP data file and show in table-------------------------  
  #Create reactive value to store "raw" file contents
  origin <- reactiveValues(fpdata = NULL, labIDdata = NULL)
  
  #Update reactive values upon file uploads
  observeEvent(input$labid_file, {
    req(input$labid_file)
    
    IDs <- read.csv(input$labid_file$datapath, header = FALSE, stringsAsFactors = FALSE, quote = "", skip = 1,
                           col.names = c("Number", "Site", "ID", "Type", "Start", "Svol", "MQvol", "Comments"),
                           colClasses = c("numeric", "factor", "character", "factor", "character", "numeric", "numeric", "character"))
    
    origin$labIDdata <- IDs
  })
  
  observeEvent(input$fpdata_file, {
    req(input$fpdata_file)
    
    data <- read.table(input$fpdata_file$datapath, header = FALSE, sep = "\t", skip = 2, 
                    col.names = c("DateTime", "Greens", "Cyano", "Diatoms", "Crypto", 
                             "#5", "#6", "#7", "Yellow", "totChla", "Transmission", 
                             "Depth", "Temperature", "GreenCells", "CyanoCells", "DiatomCells", 
                             "CryptoCells", "#5cells", "#6cells", "#7cells", "Yellow2", 
                             "totCellCt", "T700", "LED3", "LED4", "LED5", "LED6", "LED7", 
                             "LED8", "Pressure", "TLED", "TSensor"),
                     colClasses = c("character", rep("numeric", 4), rep("NULL", 3), rep("numeric", 3), rep("NULL", 21)))
    
    origin$fpdata <- data
  })
  
  #Create a datatable with the capacity to delete 10 rows (if a run is bad or should not be included)
  output$fpdata_table <- renderDT({
    validate({
      req(input$fpdata_file)
      need(nrow(origin$fpdata)%%10 == 0, message = "FP data not a multiple of 10")
    })
    DT::datatable(origin$fpdata)
  })
  
  observeEvent(input$deletefprows, {
    if(length(input$fpdata_table_rows_selected)==10){
      origin$fpdata <- origin$fpdata[as.numeric(-input$fpdata_table_rows_selected),]
    }
  })
  
  #Create a datatable with the ability to edit Sample.ID names in case of mistyping of LabID sheet
  #and delete rows (if a run is bad or should not be included)
  output$labid_table <- renderDT(
    DT::datatable(origin$labIDdata,
                  editable = list(target = "cell", disable = list(columns = c(1:2, 4:8))))
  )
  #editable cells
  observeEvent(input$labid_table_cell_edit, {
    row  <- input$labid_table_cell_edit$row
    clmn <- input$labid_table_cell_edit$col
    origin$labIDdata[row, clmn] <- input$labid_table_cell_edit$value
  })
  
  #delete rows
  observeEvent(input$deleteIDrows, {
    if(!is.null(input$labid_table_rows_selected)){
      origin$labIDdata <- origin$labIDdata[as.numeric(-input$labid_table_rows_selected),]
    }
  })
#-------------------------------------------------------------------------------

#----------------Create "Raw Results" dataframe and show in table---------------
  Raw_Results <- reactive({
    
    #Process fpdata by removing empty rows, averaging every 10 data points, and converting 
    #DateTime back to legible format
    fp_file <- origin$fpdata
    
    fp_file[!apply(fp_file == "", 1, all),]
    
    fp_file$DateTime <- as.POSIXct(fp_file$DateTime, format = "%m/%d/%Y %H:%M", tz = "EST")
    
    fp_file <- fp_file %>%
      reframe(across(names(fp_file), ~ round(colMeans(matrix(.x, nrow =10)), digits=4), .names = "{.col}"))
    fp_file$Total <- round(fp_file$Greens + fp_file$Cyano + fp_file$Diatoms + fp_file$Crypto, digits = 4)
    
    fp_file$DateTime <- as.POSIXct(fp_file$DateTime, tz = "EST", origin = "1970-01-01")
    fp_file$DateTime <- format(fp_file$DateTime, "%m/%d/%Y, %H:%M:%S")
    
    
    #Process LabID sheet by removing empty rows and calculating dilution factor
    labid_file <- origin$labIDdata
    labid_file[!apply(labid_file == "", 1, all), ]
    
    labid_file <- labid_file %>%
      mutate(dilution_factor = (Svol+MQvol)/Svol)
    
    #Adjust FP data by multiplying by dilution factor for each sample
    multiply_by_df <- fp_file %>%
      reframe(across(-c(DateTime, Transmission), ~ .x*labid_file$dilution_factor, .names = "{.col}_f"))
    
    #Combine dataframes to produce "raw results" and move Transmission column to the front
    cbind(labid_file, fp_file, multiply_by_df) %>%
      select(DateTime, Site, ID, Type, Start, dilution_factor, Greens, Cyano, Diatoms, Crypto, Yellow, Total,
             Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, Total_f, Transmission, Comments) %>%
      relocate(Transmission, .before = Type)
    
  })
  
  #Display "raw results" in a table ONLY if the original FPdata file is 10x long as LabID file.
  #Conditional formatting of Transmission column to be red if < 90
  output$rawresults_table <- renderReactable({
    validate({
      req(input$fpdata_file)
      req(input$labid_file)
      
      need(nrow(origin$fpdata)/nrow(origin$labIDdata)== 10, message = "Data and LabID files not of proportional length")
    })
    reactable(Raw_Results(),
              columns = list(
                Transmission = colDef(style = function(value){
                  if (value < 90){
                    color <- "red"
                    
                    list(color = color)
                  }
                })
              )
            )
  })
#-------------------------------------------------------------------------------

#----------------Create "Final Report" dataframe and show in table--------------
  final_report <- reactive({
    
    #Create final dataframe by removing rows with Transmission < 90, group by Sample ID and Site,
    #average results rounded to 2 digits, and finally select desired columns to be included in final output
    fp_report <- Raw_Results() %>%
      subset(Transmission >= 90) %>%
      group_by(Site, ID, .add = TRUE) %>%
      summarise(count = n(), Green_Chl=round(mean(Greens_f), 2), Bluegreen_Chl = round(mean(Cyano_f), 2), Diatom_Chl=round(mean(Diatoms_f),2 ), 
                Cryptophyte_Chl=round(mean(Crypto_f), 2), Total_Chl=round(mean(Total_f),  2), Yellow_Sub=round(mean(Yellow_f), 2)) %>%
      ungroup() %>%
      select(Site, ID, count, Green_Chl, Bluegreen_Chl, Diatom_Chl, Cryptophyte_Chl, Total_Chl, Yellow_Sub)
  })
  
 #Display table ONLY if original FPdata file is 10x long as LabID file
  output$finalreport_table <- renderDT({
    validate({
      req(input$fpdata_file)
      req(input$labid_file)
      need(nrow(origin$fpdata)/nrow(origin$labIDdata)== 10, message = "Data and LabID files not of proportional length")
    })
    DT::datatable(final_report())
  })
#-------------------------------------------------------------------------------
  #Enable saving of data only if input parameters are selected, date of FP and LabID files match, and if
  #FP file is 10x long as LabID file
  
  output$save_button <- renderUI({
    validate(
      need(input$project_sel != "", message = "Select projects"),
      need(str_sub(input$fpdata_file$name, 1, 6) == str_sub(input$labid_file$name, 1, 6), message = "Mistmatching Data and LabID files"),
      need(nrow(origin$fpdata)/nrow(origin$labIDdata)== 10, message = "")
    )
    actionButton("save_data", "Save Data", class = "btn-success")
  })
  
  #Use input parameters to create file name for "QAQC" output
  raw_filename <- reactive({
    date <- strftime(input$date_sel, format = "%y%m%d")
    
    proj <- paste(input$project_sel, collapse = "_")
    
    if(isTruthy(input$comment)==TRUE){
      comment <- paste0(toupper(input$comment), "_")
    } else comment = input$comment
    
    paste0("QAQC","_", date, "_", comment, proj)
  })
  
  #Use input parameters to create file name for final report output
  fp_filename <- reactive({
    date <- strftime(input$date_sel, format = "%y%m%d")
    
    proj <- paste(input$project_sel, collapse = "_")
    
    if(isTruthy(input$comment)==TRUE){
      comment <- paste0(toupper(input$comment), "_")
    } else comment = input$comment
    
    paste0(date, "_", comment, proj)
  })
  
  #Save data upon clicking save button
  observeEvent(input$save_data, {
    req(input$labid_file, input$fpdata_file)
    
    write.csv(Raw_Results(), paste0("C:/Users/KELLEY/Documents/R Main Directory/Fluoroprobe/QAQC/",raw_filename(),".csv"), row.names = FALSE)
    write.csv(final_report(), paste0("C:/Users/KELLEY/Documents/R Main Directory/Fluoroprobe/final_reports/",fp_filename(),".csv"), row.names = FALSE)
  })
  
  #Reset file params upon file uploads
  observeEvent(input$fpdata_file, {
    reset("File_params")
  })
  
  observeEvent(input$labid_file, {
    reset("File_params")
  })
}

shinyApp(ui = ui, server = server)