merge_fp_id <- function(fp_data, id_data){
  fp_data %>%
    mutate(DateTime = as.POSIXct(DateTime, format = "%m/%d/%Y %H:%M", tz = "EST")) %>%
    reframe(across(1:8, ~ matrix(.x, nrow =10) %>% colMeans() %>% round(4), .names = "{.col}")) %>%
    mutate(DateTime = as.POSIXct(DateTime, tz = "EST", origin = "1970-01-01") %>% format("%m/%d/%Y, %H:%M:%S")) %>%
    cbind(id_data) %>%
    mutate(dilution_factor = (Svol+MQvol)/Svol %>% round(2),
           across(Greens:totChla, ~.x * dilution_factor, .names = "{.col}_f")) %>%
    select(DateTime, Site, ID, Transmission, Type, Start, dilution_factor, Greens, Cyano, Diatoms, Crypto, Yellow, totChla,
           Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, totChla_f, Comments) %>%
    filter(!grepl("NV", Comments))
}

process_duplicates <- function(data){
  data %>%
    select(DateTime, Site, ID, Transmission, Type, Start, dilution_factor, totChla, totChla_f, Comments) %>%
    filter(!grepl("Blank", Type, ignore.case = T)) %>%
    group_by(ID) %>%
    mutate(RPD_Total = ifelse(length(ID)==2 & grepl("Duplicate", Comments, ignore.case = T), round((max(totChla)-min(totChla))/mean(totChla)*100,2), NA)) %>%
    mutate(RPD_Total_F = ifelse(length(ID)==2 & grepl("Duplicate", Comments, ignore.case = T), round((max(totChla_f)-min(totChla_f))/mean(totChla_f)*100,2), NA)) %>%
    ungroup() %>%
    mutate(Date_Analyzed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Analyzed, DateTime, Site, Type, dilution_factor, Start, Transmission, totChla, totChla_f, Comments, RPD_Total, RPD_Total_F) %>%
    rename(Sample.ID = "ID")
}

process_blanks <- function(data){
  data %>%
    filter(grepl("Blank", Type, ignore.case = T)) %>%
    mutate(Date_Analyzed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Analyzed, DateTime, Site, Type, Start, dilution_factor, Transmission, Greens, Cyano, Diatoms, Crypto, Yellow, totChla, Comments) %>%
    rename(Sample.ID = "ID")
}

process_field_blanks <- function(data){
  data %>%
    filter(grepl("FB", ID)) %>%
    mutate(Date_Analyzed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Analyzed, DateTime, Site, Type, Start, dilution_factor, Transmission, Greens, Cyano, Diatoms, Crypto, Yellow, totChla, Comments) %>%
    rename(Sample.ID = "ID")
}

process_output_preview <- function(data){
  output <- data %>%
    filter(Transmission >= 90, !grepl("MQ", Site), !grepl("FB", ID)) %>%
    group_by(ID) %>%
    mutate(Count = n(), 
           across(Greens_f:totChla_f, ~ mean(.x) %>% round(2))) %>%
    ungroup() %>%
    mutate(Date_Analyzed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Analyzed, Site, Type, Count, Greens_f, Cyano_f, Diatoms_f, Crypto_f, totChla_f) %>%
    unique()
    
    names(output) <- c("Sample.ID", "Date_Analyzed", "Site", "Type", "Count", "Green_Chl", "Bluegreen_Chl", "Diatom_Chl", "Cryptophyte_Chl", "Total_Chl")

    return(output)
}

file_name <- function(date, description){
  if(nchar(description) < 1){
    name <- paste(date, "FP", sep = "_")
  } else{
    name <- paste(date, description, "FP", sep = "_")
  }
  return(name)
}