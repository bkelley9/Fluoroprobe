##2022 Script to read in processed lab, id files for HAB report

#load packages
library(dplyr)
library(readr)
library(data.table)
library(arsenal)
library(readr)
library(ggplot2)
library(lubridate)
library(ggfortify)
library(tidyverse)
library(data.table)
library(formattable)
library(stringr)

##set working directory; where you're getting files from
#getwd()
#will be user dependent
#setwd(dir = "C:/Users/CLECKNER/Documents/2022 HABs/RProcResults")

#fp_report2 <- fread("Combined/Results_221006.csv")
fp_report2 <- fread("Combined/Results22_221021.csv")
fp_report2 <- fread("Missing JH files Results.csv")
glimpse(fp_report2)

colnames(fp_report2) <- c("V1", "DateTime", "Greens", "Site", "ID", "Type", "Start", "Dil. Factor.", 
    "AnGreens", "AnCyanos", "AnDiatoms", "AnCryptos", "AnYellow", "AnTotal", 
    "FinalGreen", "FinalCyano", "FinalDiatom", "FinalCrypto", "FinalYellow", 
    "FinalTotal", "Transmission", "Comments")
glimpse(fp_report2)


#ensure transmission > 90, samples valid
valid <- subset(fp_report2, fp_report2$Transmission >= 90 & fp_report2$Comments != "NV") 
glimpse(valid)
str(valid)

valid$DateTime <- as.factor(valid$DateTime)
valid$Site <- as.factor(valid$Site)
valid$ID <- as.factor(valid$ID)
valid$Type <- as.factor(valid$Type)
valid$FinalGreen <- as.numeric(valid$FinalGreen)
valid$FinalCyano <- as.numeric(valid$FinalCyano)
valid$FinalDiatom <- as.numeric(valid$FinalDiatom)
valid$FinalCrypto <- as.numeric(valid$FinalCrypto)

#screen out blanks
validNB <- subset(valid, Type != "Blank")

#need to add date into file
fp_valid <- validNB %>% 
  select(DateTime, Site, ID, Type, FinalGreen, FinalCyano, FinalDiatom, FinalCrypto, FinalTotal, FinalYellow) %>%
  group_by(Site, ID, Type, .add = FALSE) %>%
  summarize(Green_Chl=mean(FinalGreen), Bluegreen_Chl = mean(FinalCyano), Diatom_Chl=mean(FinalDiatom), 
            Cryptophyte_Chl=mean(FinalCrypto), Total_Chl=mean(FinalTotal), Yellow_Sub=mean(FinalYellow))
glimpse(fp_valid)
str(fp_valid)

write.csv(fp_valid, file ="Combined/MissingJH221016.csv", row.names=TRUE)

HAB_report <- subset(fp_valid, fp_valid$Type == "HAB")

#drop YS from HAB report
HAB_report <- subset(HAB_report, select = -9)

#Round chl and yellow_sub values to 2 decimal places
HAB_report$Green_Chl <- round(HAB_report$Green_Chl, digits=2)
HAB_report$Bluegreen_Chl <- round(HAB_report$Bluegreen_Chl, digits=2)
HAB_report$Diatom_Chl <- round(HAB_report$Diatom_Chl, digits=2)
HAB_report$Cryptophyte_Chl <- round(HAB_report$Cryptophyte_Chl, digits=2)
HAB_report$Total_Chl <- round(HAB_report$Total_Chl, digits=2)
View(HAB_report)

colnames(HAB_report) <- c("Site", "ID", "Type", "GreenChl (µg/L)", "CyanoChl (µg/L)",
                          "DiatomChl (µg/L)", "CryptophyteChl (µg/L)", "Total Chlorophyll (µg/L)")

KE_chl <- subset(HAB_report, Site=="KE")
write.csv(KE_chl, file ="LakeReports/220817KEchl.csv", row.names=TRUE)

CN_chl <- subset(HAB_report, Site=="CN")
write.csv(CN_chl, file ="LakeReports/220823CNchl.csv", row.names=TRUE)


#filter by lake


#save csv HAB
write.csv(HAB_report, file ="LakeReports/KE_HAB_220816.csv", row.names=TRUE)



#library(plyr)

#combine results files
setwd("C:/Users/CLECKNER/Documents/2022 HABs/RProcResults/merged files")
getwd()
mypath = "C:/Users/CLECKNER/Documents/2022 HABs/RProcResults/merged files"

all_res <- list.files(mypath, pattern = '*.csv', full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()

#delete first column
all_results <- subset(all_res, select=-1)
#write code to rename 1st column as Sample Number or something
write.csv(all_results, file = "C:/Users/CLECKNER/Documents/2022 HABs/RProcResults/Combined/Results22_221021.csv", row.names=TRUE)





#make fp_report data frame with valid data
#average duplicates (rows with same ID) and include site, sample type as well)
fp_sample <- fp_report %>% 
  group_by(ID,Site,Type,.add=TRUE) %>%
  summarize(Green_Chl=mean(Green_Chl), Bluegreen_Chl = mean(Bluegreen_Chl), Diatom_Chl=mean(Diatom_Chl), 
            Cryptophyte_Chl=mean(Cryptophyte_Chl), Total_Chl=mean(FinalTotal), Yellow_Sub=mean(FinalYellow))
View(fp_sample)

HAB_report <- subset(fp_sample, fp_sample$Type == "HAB")


#Round chl and yellow_sub values to 2 decimal places
HAB_report$Green_Chl <- round(fp_report$Green_Chl, digits=2)
HAB_report$Bluegreen_Chl <- round(fp_report$Bluegreen_Chl, digits=2)
HAB_report$Diatom_Chl <- round(fp_report$Diatom_Chl, digits=2)
HAB_report$Cryptophyte_Chl <- round(fp_report$Cryptophyte_Chl, digits=2)
HAB_report$Total_Chl <- round(fp_report$Total_Chl, digits=2)
View(HAB_report)






##fp_data______________________________________________
#reading in raw FP text file, deleting first two rows, adding column names, subsetting variables/columns relevant to analysis 
#converting variables to working format (e.g. time, numeric), averaging 10 measurements to get one value per sample

#need to change name of input FP file
fp <- read.table("220720_22_fpdata.txt", header = FALSE, sep = "\t", skip = 2)
str(fp)

#add column names to FP file
colnames(fp) <- c("DateTime", "Greens", "Cyano", "Diatoms", "Crypto", "#5", "#6", "#7", "Yellow", "totChla", "Transmission", "Depth", "Temperature", "GreenCells", "CyanoCells", "DiatomCells", "CryptoCells", "#5cells", "#6cells", "#7cells", "Yellow2", "totCellCt", "T700", "LED3", "LED4", "LED5", "LED6", "LED7", "LED8", "Pressure", "TLED", "TSensor")

fp


#Subset data from relevant columns
calcdata <- fp[ ,c(1:5,9:11)]
str(calcdata)


#inspect format of DateTime variable
tail(calcdata$DateTime)
#change DateTime written format to "YYYY-MM-DD HH:MM:SS EDT", but stays as character
tail(as.POSIXct(calcdata$DateTime, format = "%m/%d/%Y %H:%M"))
#create Timestamp column, with date and time in POSIXct format (R recognizes it as actual date and time, rather than reading as characters)
calcdata$Timestamp <- as.POSIXct(calcdata$DateTime, format = "%m/%d/%Y %H:%M")
str(calcdata)

#Convert characters to numeric so the data can be averaged below
calcdata$Greens <- as.numeric(calcdata$Greens)
calcdata$Cyano <- as.numeric(calcdata$Cyano)
calcdata$Diatoms <- as.numeric(calcdata$Diatoms)
calcdata$Crypto <- as.numeric(calcdata$Crypto)
calcdata$Yellow <- as.numeric(calcdata$Yellow)
calcdata$totChla <- as.numeric(calcdata$totChla)
calcdata$Transmission <- as.numeric(calcdata$Transmission)

#check calcdata structure
str(calcdata)

#Average data for every 10 rows
timeavg <- colMeans(matrix(calcdata$Timestamp, nrow=10))
timeavg
greensavg <- colMeans(matrix(calcdata$Greens, nrow=10))
greensavg
cyanoavg <- colMeans(matrix(calcdata$Cyano, nrow=10))
diatomsavg <- colMeans(matrix(calcdata$Diatoms, nrow=10))
cryptoavg <- colMeans(matrix(calcdata$Crypto, nrow=10))
totalavg <- greensavg + cyanoavg + diatomsavg + cryptoavg
yelSubavg <- colMeans(matrix(calcdata$Yellow, nrow=10))
transavg <- colMeans(matrix(calcdata$Transmission, nrow=10))

#make new dataframe with averaged values
datfp <- data.frame(timeavg, greensavg, cyanoavg, diatomsavg, cryptoavg, totalavg, yelSubavg, transavg)
head(datfp)
str(datfp)

#need to convert averaged date/time back to POSIXct DateTime format 
time.cv <- as.POSIXct(datfp$timeavg, tz = "", origin = "1970-01-01")  
tail(time.cv)
str(time.cv)

#make dataframe out of data frame with averaged data and POSIXct date/time
fpcalc <- cbind(time.cv, datfp)
View(fpcalc)



##labid____________________________________________________
#read in typed Lab ID sheet, add correct column names,convert some variables to numeric
#calculate dilution factor and multiply analyzed concentrations by it to get final concentrations

#read in csv file of labid 
#need to update name of input file
id <- read.csv("220720_22_labid.csv", header = FALSE,stringsAsFactors = F, quote = "")
str(id)

##delete first row of labid (original column names in csv file)
id = id[-1,]

#add column names to labid 
colnames(id) <- c("Number", "Site", "ID", "Type", "Start", "Svol", "MQvol", "Comments")

View(id)

#delete empty rows
id <- id[!apply(id == "", 1, all),]

#convert sample volume (Svol) and milli-q volume (MQvol) to numeric 
id$Svol <- as.numeric(id$Svol)
id$MQvol <- as.numeric(id$MQvol)
str(id)

#calculate dilution factor based on mL sample and mL MQ Water
dilf <- 1/(id$Svol/(id$Svol + id$MQvol))

#calculate dilution factor (should be 1 if sample not diluted)
x3 <- cbind(id, dilf, fpcalc)

#get final concentrations by multiplying dilution factor * average concentration
greenfinal <- x3$dilf * x3$greensavg
cyanofinal <- x3$dilf * x3$cyanoavg
diatomfinal <- x3$dilf * x3$diatomsavg
cryptofinal <- x3$dilf * x3$cryptoavg
yellowfinal <- x3$dilf * x3$yelSubavg
totalfinal <- x3$dilf * x3$totalavg

#make new dataframe
fpfinal <- data.frame(greenfinal, cyanofinal, diatomfinal, cryptofinal, yellowfinal, totalfinal)

#make new dataframe combining labid info, dilution factor, calculated averages, and final concentrations after multiplying by dilution factor
x4 <- cbind(id,dilf, fpcalc, fpfinal)


##Write raw_results_______________________________________
#raw results include the analyzed and calculated values of each parameter,
#duplicates not averaged and not all samples are valid (transmission may be <90 percent for some samples)

#compile Raw Results with desired columns
Raw_Results <- data.frame(x4$time.cv, x4$Site, x4$ID, x4$Type, x4$Start, x4$dilf, x4$greensavg, x4$cyanoavg, x4$diatomsavg, x4$cryptoavg, x4$yelSubavg, x4$totalavg, x4$greenfinal, x4$cyanofinal, x4$diatomfinal, x4$cryptofinal, x4$yellowfinal, x4$totalfinal, x4$transavg, x4$Comments)  
#add column names to report
colnames(Raw_Results) <- c("DateTime", "Site", "ID", "Type", "Start", "Dil. Factor", "AnGreens", "AnCyanos", "AnDiatoms", "AnCryptos", "AnYellow","AnTotal",  
                           "FinalGreen", "FinalCyano", "FinalDiatom", "FinalCrypto", "FinalYellow","FinalTotal", "Transmission",  "Comments")
View(Raw_Results)

#Write and export raw_results .csv (duplicates not averaged and transmission not filtered out) CSV in R
#need to update name of Results file
write.csv(Raw_Results, file = "Results_220720_22.csv", row.names=TRUE)

#stop here for daily reports

##Write fp_report__________________________________________
#fp_report includes: the mean values for each ID/site/type of: 
#FinalGreen, FinalCyano, FinalDiatom, FinalCrypto, FinalTotal, and FinalYellow.

#filter out rows without transmission >= 90 to get valid data
trans_filter <- subset(Raw_Results, x4$transavg >= 90) 
View(trans_filter)


#make fp_report data frame with valid data
#average duplicates (rows with same ID) and include site, sample type as well)
fp_report <- trans_filter %>% 
  group_by(ID,Site,.add=TRUE) %>%
  summarize(Green_Chl=mean(FinalGreen), Bluegreen_Chl = mean(FinalCyano), Diatom_Chl=mean(FinalDiatom), 
            Cryptophyte_Chl=mean(FinalCrypto), Total_Chl=mean(FinalTotal), Yellow_Sub=mean(FinalYellow))
View(fp_report)

#Round chl and yellow_sub values to 2 decimal places
fp_report$Green_Chl <- round(fp_report$Green_Chl, digits=2)
fp_report$Bluegreen_Chl <- round(fp_report$Bluegreen_Chl, digits=2)
fp_report$Diatom_Chl <- round(fp_report$Diatom_Chl, digits=2)
fp_report$Cryptophyte_Chl <- round(fp_report$Cryptophyte_Chl, digits=2)
fp_report$Total_Chl <- round(fp_report$Total_Chl, digits=2)
fp_report$Yellow_Sub <- round(fp_report$Yellow_Sub, digits=2)
View(fp_report)

#write and export fp_report .csv in R
write.csv(fp_report, file = "fp_Report_220712.csv", row.names=TRUE)

#read in processed data from Box


#ADD TO CSV REPOSITORIES____________________________________
#code from AW to combine all csvs in a folder
library(plyr)

#labid
mypath = "C:/Users/johns/Documents/FLI/Fluoroprobe/labid"
all_labid <- list.files(mypath, pattern = '*.csv', full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
#write code to rename 1st column as Sample Number or something
write.csv(all_labid, file = "Repositories/labid_2021.csv", row.names=TRUE)


#fpdata
mypath = "C:/Users/johns/Documents/FLI/Fluoroprobe/fpdata"
all_fpdata <- list.files(mypath, pattern = '*.csv', full.names = TRUE) %>%
  lapply(read.csv) %>% #Need to find a way to delete the first row of each data file as it's read in
  bind_rows()
#delete first row
all_fpdata = all_fpdata[-1,]
#rename columns
#add column names to FP file
colnames(all_fpdata) <- c("DateTime", "Greens", "Cyano", "Diatoms", "Crypto", "#5", "#6", "#7", "Yellow", "totChla", "Transmission", "Depth", "Temperature", "GreenCells", "CyanoCells", "DiatomCells", "CryptoCells", "#5cells", "#6cells", "#7cells", "Yellow2", "totCellCt", "T700", "LED3", "LED4", "LED5", "LED6", "LED7", "LED8", "Pressure", "TLED", "TSensor")
write.csv(all_fpdata, file = "Repositories/fpdata_2021.csv", row.names=TRUE)


#raw results
mypath = "C:/Users/johns/Documents/FLI/Fluoroprobe/fpresults/raw_results"
all_raw_results <- list.files(mypath, pattern = '*.csv', full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
write.csv(all_raw_results, file = "Repositories/fp_results_2021.csv", row.names=TRUE)


#finished reports
mypath = "C:/Users/johns/Documents/FLI/Fluoroprobe/fpresults/fp_reports"
all_fin_reports <- list.files(mypath, pattern = '*.csv', full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
write.csv(all_fin_reports, file = "Repositories/fp_report_2021.csv", row.names=TRUE)



#______________________________________________________________
#AW edits making report from raw results 08/09/2021
#making report from results csv that Lisa sent for Dock-a-palooza
#Results does not have yellow substances data in it
Raw_Results <- read.csv("Results_210805dp.csv")

#Can only skip over trans_filter step if all transmission is 100 
#averaging duplicates
fp_report <- Raw_Results %>% 
  group_by(ID,Site,Type,.add=TRUE) %>%
  summarize(Green_Chl=mean(FinalGreen), Bluegreen_Chl = mean(FinalCyano), Diatom_Chl=mean(FinalDiatom), Cryptophyte_Chl=mean(FinalCrypto), Total_Chl=mean(FinalTotal))
View(fp_report)

#Round chl and yellow_sub values to 2 decimal places
fp_report$Green_Chl <- round(fp_report$Green_Chl, digits=2)
fp_report$Bluegreen_Chl <- round(fp_report$Bluegreen_Chl, digits=2)
fp_report$Diatom_Chl <- round(fp_report$Diatom_Chl, digits=2)
fp_report$Cryptophyte_Chl <- round(fp_report$Cryptophyte_Chl, digits=2)
fp_report$Total_Chl <- round(fp_report$Total_Chl, digits=2)

#write and export fp_report .csv in R
write.csv(fp_report, file = "C:\\Users\\Owner\\Desktop\\Razavi Lab\\Abby\\ESF PhD\\FLI summer 2021\\fpdata\\fp_Report_210811.csv", row.names=TRUE)
