merge_fp_id <- function(fp_data, id_data){
  fp_data %>%
    mutate(DateTime_Analyzed = as.POSIXct(DateTime_Analyzed, format = "%m/%d/%Y %H:%M", tz = "EST")) %>%
    reframe(across(1:8, ~ matrix(.x, nrow =10) %>% colMeans() %>% round(4), .names = "{.col}")) %>%
    mutate(DateTime_Analyzed = as.POSIXct(DateTime_Analyzed, tz = "EST", origin = "1970-01-01") %>% format("%m/%d/%Y, %H:%M:%S")) %>%
    cbind(id_data) %>%
    mutate(dilution_factor = (Svol+MQvol)/Svol %>% round(2),
           across(Greens:totChla, ~.x * dilution_factor, .names = "{.col}_f")) %>%
    select(DateTime_Analyzed, Site, ID, Transmission, Type, Start, dilution_factor, Greens, Cyano, Diatoms, Crypto, Yellow, totChla,
           Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, totChla_f, Comments) %>%
    filter(!grepl("NV", Comments))
}

process_duplicates <- function(data){
  data %>%
    select(DateTime_Analyzed, Site, ID, Transmission, Type, Start, dilution_factor, totChla, totChla_f, Comments) %>%
    filter(!grepl("Blank", Type, ignore.case = T)) %>%
    group_by(ID) %>%
    mutate(RPD_CV_Total = ifelse(length(ID) > 2, round(sd(totChla)/mean(totChla) * 100, 2), ifelse(length(ID)==2 & grepl("Duplicate", Comments, ignore.case = T), round((max(totChla)-min(totChla))/mean(totChla)*100,2), NA))) %>%
    mutate(RPD_CV_Total_F = ifelse(length(ID) > 2, round(sd(totChla_f)/mean(totChla_f) * 100, 2), ifelse(length(ID)==2 & grepl("Duplicate", Comments, ignore.case = T), round((max(totChla_f)-min(totChla_f))/mean(totChla_f)*100,2), NA))) %>%
    ungroup() %>%
    mutate(Date_Processed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Processed, DateTime_Analyzed, Site, Type, dilution_factor, Start, Transmission, totChla, totChla_f, Comments, RPD_CV_Total, RPD_CV_Total_F) %>%
    rename(Sample.ID = "ID")
}

process_blanks <- function(data){
  data %>%
    filter(grepl("Blank", Type, ignore.case = T)) %>%
    mutate(Date_Processed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Processed, DateTime_Analyzed, Site, Type, Start, dilution_factor, Transmission, Greens, Cyano, Diatoms, Crypto, Yellow, totChla, Comments) %>%
    rename(Sample.ID = "ID")
}

process_field_blanks <- function(data){
  data %>%
    filter(grepl("FB", ID)) %>%
    mutate(Date_Processed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Processed, DateTime_Analyzed, Site, Type, Start, dilution_factor, Transmission, Greens, Cyano, Diatoms, Crypto, Yellow, totChla, Comments) %>%
    rename(Sample.ID = "ID")
}

process_output_preview <- function(data){
  output <- data %>%
    filter(Transmission >= 90, !grepl("MQ", Site), !grepl("FB", ID)) %>%
    group_by(ID) %>%
    mutate(Count = n(), 
           across(Greens_f:totChla_f, ~ mean(.x) %>% round(2))) %>%
    ungroup() %>%
    mutate(Date_Processed = Sys.Date() %>% format("%m/%d/%Y")) %>%
    select(ID, Date_Processed, DateTime_Analyzed, Site, Type, Count, Greens_f, Cyano_f, Diatoms_f, Crypto_f, totChla_f, Yellow_f) %>%
    distinct(ID, .keep_all = T)
    
    names(output) <- c("Sample.ID", "Date_Processed", "DateTime_Analyzed", "Site", "Type", "Count", "Green_Chl", "Bluegreen_Chl", "Diatom_Chl", "Cryptophyte_Chl", "Total_Chl", "Yellow_subs")

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