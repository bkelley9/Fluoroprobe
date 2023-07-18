library(shiny)
library(tidyverse)
library(DT)
library(formattable)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  
  tags$head(
    tags$style(HTML(".shiny-output-error-validation {
                    color: red; }"))
  ),
  titlePanel("Fluoroprobe Analysis"),
  sidebarLayout(
    sidebarPanel(
      div(id = "File_Inputs",
          fileInput("fpdata", "Select FPdata file", accept = ".txt"),
          fileInput("labid", "Select LabID file", accept = ".csv")
          ),
      
      div(id = "File_params",
          dateInput("date_sel", "Select Date", value = NULL, format = "yy/mm/dd"),
          checkboxGroupInput("project_sel", "Select Relevant Project/Lake",
                             choices = c("Owasco" = "OW",
                                         "Honeoye" = "HN",
                                         "Seneca"= "SC",
                                         "KLA",
                                         "MCOW"))
          ),
      uiOutput("save_button"),
      span(textOutput("ready"), style = "color:green")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("FP Data",
          tableOutput("fpdata_table")
        ),
        tabPanel("Lab ID",
          tableOutput("labid_table")
        ),
        tabPanel("Raw Results",
                 DT::dataTableOutput("rawresults_table")
        ),
        tabPanel("FP Report",
                 DT::dataTableOutput("finalreport_table"))
      )
    )
  )
)

server <- function(input, output){

#-------------------Read FP data file and show in table-------------------------  
  fp_file <- reactive({
    req(input$fpdata)
    
    read.table(input$fpdata$datapath, header = FALSE, sep = "\t", skip = 2, 
               col.names = c("DateTime", "Greens", "Cyano", "Diatoms", "Crypto", 
                             "#5", "#6", "#7", "Yellow", "totChla", "Transmission", 
                             "Depth", "Temperature", "GreenCells", "CyanoCells", "DiatomCells", 
                             "CryptoCells", "#5cells", "#6cells", "#7cells", "Yellow2", 
                             "totCellCt", "T700", "LED3", "LED4", "LED5", "LED6", "LED7", 
                             "LED8", "Pressure", "TLED", "TSensor"),
               colClasses = c("character", rep("numeric", 4), rep("NULL", 3), rep("numeric", 3), rep("NULL", 21)))
  })
  
  output$fpdata_table <- renderTable(
    fp_file()
  )
#-------------------------------------------------------------------------------
  
#-----------------Read Lab ID file and show in table----------------------------
  labid_file <- reactive({
    req(input$labid)
    
    read.csv(input$labid$datapath, header = FALSE, stringsAsFactors = FALSE, quote = "", skip = 1,
              col.names = c("Number", "Site", "ID", "Type", "Start", "Svol", "MQvol", "Comments"),
             colClasses = c("numeric", "factor", "factor", "factor", "character", "numeric", "numeric", "character"))
  })

  output$labid_table <- renderTable(
    labid_file()
  )
#-------------------------------------------------------------------------------

#----------------Create "Raw Results" dataframe and show in table---------------
  Raw_Results <- reactive({
    fp_raw <- fp_file()
    
    fp_raw$DateTime <- as.POSIXct(fp_raw$DateTime, format = "%m/%d/%Y %H:%M")
    
    fp_raw <- fp_raw %>%
      reframe(across(names(fp_raw), ~ round(colMeans(matrix(.x, nrow =10)), digits=4), .names = "{.col}"))
    fp_raw$Total <- round(fp_raw$Greens + fp_raw$Cyano + fp_raw$Diatoms + fp_raw$Crypto, digits = 4)
    fp_raw$DateTime <- as.POSIXct(fp_raw$DateTime, tz = "EST", origin = "1970-01-01")
    
    labid_raw <- labid_file()
    labid_raw <- labid_raw %>%
      mutate(dilution_factor = (.$Svol+.$MQvol)/.$Svol)
    
    multiply_by_df <- fp_raw %>%
      reframe(across(-c(DateTime, Transmission), ~ .x*labid_raw$dilution_factor, .names = "{.col}_f"))
    multiply_by_df$DateTime <- fp_raw$DateTime
    multiply_by_df$Transmission <- fp_raw$Transmission
    
    Raw_Results_df <- cbind(labid_raw, fp_raw, subset(multiply_by_df, select = -c(DateTime, Transmission))) %>%
      select(DateTime, Site, ID, Type, Start, dilution_factor, Greens, Cyano, Diatoms, Crypto, Yellow, Total,
             Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, Total_f, Transmission, Comments)
    
  })
  
  output$rawresults_table <- renderDT(
    DT::datatable(Raw_Results())
  )
#-------------------------------------------------------------------------------

#----------------Create "Final Report" dataframe and show in table--------------
  final_report <- reactive({
    Raw_df <- Raw_Results()
    
    filter_dups <- Raw_df %>%
      subset(Transmission >= 90) %>%
      group_by(Site, ID, .add = TRUE) %>%
      summarise(Green_Chl=mean(Greens_f), Bluegreen_Chl = mean(Cyano_f), Diatom_Chl=mean(Diatoms_f), 
                Cryptophyte_Chl=mean(Crypto_f), Total_Chl=mean(Total_f), Yellow_Sub=mean(Yellow_f)) %>%
      ungroup()
    
    fp_report <- filter_dups %>%
      reframe(across(-c(Site, ID), ~ round(.x, digits = 2)))
    fp_report$Site <- filter_dups$Site
    fp_report$ID <- filter_dups$ID
    
    fp_report <- fp_report %>%
      select(Site, ID, Green_Chl, Bluegreen_Chl, Diatom_Chl, Cryptophyte_Chl, Total_Chl, Yellow_Sub)
  })
  
 
  output$finalreport_table <- renderDT(
    DT::datatable(final_report())
  )
#-------------------------------------------------------------------------------
  output$save_button <- renderUI({
    validate(
      need(input$project_sel != "", message = "Need to select projects")
    )
    actionButton("save_data", "Save Data", class = "btn-success")
  })
  
  raw_filename <- reactive({
    date <- strftime(input$date_sel, format = "%y%m%d")
    
    proj <- paste(input$project_sel, collapse = "_")
    
    paste0("Raw_Results","_", date, "_", proj)
  })
  
  fp_filename <- reactive({
    date <- strftime(input$date_sel, format = "%y%m%d")
    
    proj <- paste(input$project_sel, collapse = "_")
    
    paste0(date, "_", proj)
  })
  
  
  observeEvent(input$save_data, {
    req(input$labid, input$fpdata)
    
    write.csv(Raw_Results(), paste0("C:/Users/KELLEY/Documents/R Main Directory/Fluoroprobe/raw_results/",raw_filename(),".csv"), row.names = FALSE)
    write.csv(final_report(), paste0("C:/Users/KELLEY/Documents/R Main Directory/Fluoroprobe/final_reports/",fp_filename(),".csv"), row.names = FALSE)
  })
  
  observeEvent(input$fpdata, {
    reset("File_params")
  })
  
}

shinyApp(ui = ui, server = server)