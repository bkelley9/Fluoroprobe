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
    mutate(RPD_Total_F = ifelse(length(ID)==2 & grepl("Duplicate", Comments, ignore.case = T), round((max(totChla_f)-min(totChla_f))/mean(totChla_f)*100,2), NA))
}

process_blanks <- function(data){
  data %>%
    filter(grepl("Blank", Type, ignore.case = T)) %>%
    select(-c(Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, totChla_f))
}

process_field_blanks <- function(data){
  data %>%
    filter(grepl("FB", ID)) %>%
    select(-c(Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, Total_f))
}

process_output_preview <- function(data){
  output <- data %>%
    filter(Transmission >= 90, !grepl("MQ", Site), !grepl("FB", ID)) %>%
    group_by(ID) %>%
    mutate(Count = n(), 
           across(Greens_f:totChla_f, ~ mean(.x) %>% round(2))) %>%
    select(Site, ID, Type, Count, Greens_f, Cyano_f, Diatoms_f, Crypto_f, totChla_f)
    
    names(output) <- c("Site", "ID", "Type", "Count", "Green_Chl", "Bluegreen_Chl", "Diatom_Chl", "Cryptophyte_Chl", "Total_Chl")

    return(output)
}

process_blanks <- function(data){
  output <- data %>%
    filter(grepl("Blank", Type, ignore.case = T)) %>%
    select(-c(Greens_f, Cyano_f, Diatoms_f, Crypto_f, Yellow_f, totChla_f))
  
  return(output)
}
