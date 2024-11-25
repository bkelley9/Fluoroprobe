# Fluoroprobe Data Processing App

## Before the Run
Each run on the benchtop Fluoroprobe (FP) must be accompanied by a bench sheet with relevant metadata for the run. These benchsheets must contain the following columns in this order:
- Number
- Site
- ID (Sample ID)
- Type (Raw Water, HAB, etc.)
- Start (time of the first record in the FP data file for the respective sample)
- Svol (volume of sample in cuvette)
- MQvol (volume of MQ water in cuvette)
- Comments
  - For sample duplicates, enter 'duplicate' in **both** rows where the duplicate occurs
  - For samples that are out of Transmission range (ie. Transmission < 90), enter 'NV' for 'Not Valid'
  
Additionally, the FP data file should always be titled 'YYMMDD_fpdata.txt' where the date corresponds to the day the FP data was collected, and the bench sheet file should be 'YYMMDD_labid.csv' with an identical date as the FP data file.

## Processing Data
The app is expecting an FP data file that is **exactly 10x** the length of the benchsheet file (i.e. if the benchsheet is 6 rows long, the FP data file should be 60). This is because we take an average of 10 measurements for each sample run on the benchtop Fluoroprobe.
- **NOTE**: occasionally excel files will save empty rows at the bottom of the file. If the app isn't working, this is likely why, and you should erase the empty rows using the 'erase' tool in excel or make a new file.

If the FP data and bench sheet files are formatted correctly, tables will generate once both files have been uploaded. QAQC will be performed on all samples by checking if the average Transmission is > 90, and if the relative percent difference (RPD) between duplicates is < 15%. Samples with Transmission < 90 will be removed from the final processed output.