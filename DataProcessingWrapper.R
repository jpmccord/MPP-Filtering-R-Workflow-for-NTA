#James McCord
#November 26, 2017
 
require("tidyverse") || install.packages("tidyverse")
require("stringr") || install.packages("stringr")
require("sqldf") || install.packages("sqldf")
require("readxl") || install.packages("readxl")
require("data.table") || install.packages("data.table")
require("tcltk") || install.packages("tcltk")
library("data.table")
library("tidyverse")
library("stringr")
library("sqldf")
library("readxl")
library("tcltk")

####### Create global functions for use later ######
SampleInput <- # Generic data input from multiple file types
  function(filename){
    inputFilename <- c(filename)
    
    if (grepl("\\.(tsv|txt)$",inputFilename)) {
      Data <- read_tsv(inputFilename, comment = "#", na = c("NA"," ", "-"), guess_max=100000)
    }
    
    if (grepl("\\.(csv)$",inputFilename)) {
      Data <- read_csv(inputFilename, comment = "#", na = c("NA"," ", "-"), guess_max=100000)
    }
    
    if (grepl("\\.(xls|xlsx)$",inputFilename) | grepl(".xls$",inputFilename)) {
      Data <- read_excel(inputFilename, col_types = NULL, na = c("NA"," ", "-", ""), guess_max=100000)
    }
    
    return(Data)
  }

adductSearchSQL <- # ID RT and Mass clustering with SQL
  function(input,
           adductMass,
           masserror,
           timetolerance,
           ppm = TRUE) {
    if (ppm == FALSE) {
      upper_limit = adductMass + (masserror/1000)
      lower_limit = adductMass - (masserror/1000)
    }
    if (ppm == TRUE) {
      new_center = max(input$Mass)
      upper_limit = adductMass + (new_center * masserror * 10 ^ -6)
      lower_limit = adductMass - (new_center * masserror * 10 ^ -6)
    }
    input_filtered <- input %>%
      select(Mass, RetentionTime, FeatureNumber, IonMode)
    matchSubset <-
      sqldf(
        paste0(
          "SELECT a.*, b.FeatureNumber as adduct_FeatureNumber,
          b.IonMode as adduct_IonMode,
          b.Mass as adduct_Mass,
          b.RetentionTime as adduct_RetentionTime
          FROM input_filtered as a, input_filtered as b
          WHERE a.FeatureNumber != b.FeatureNumber
          and abs(a.RetentionTime-b.RetentionTime)<=",timetolerance,
          " and (b.Mass-a.Mass)>=",lower_limit,
          " and (b.Mass-a.Mass)<=",upper_limit,
          " order by a.RetentionTime, b.RetentionTime"
        )
        ) 
    
    if (ppm == TRUE) {
      matchSubset <- matchSubset %>%
      mutate(ppm_error = as.numeric((abs(adduct_Mass-Mass)-adductMass)/Mass *  10 ^ 6)) %>%
      filter(ppm_error <= masserror)
    }
    
    output <-
      mutate(
        input_filtered,
        is_adduct = (FeatureNumber %in% matchSubset$adduct_FeatureNumber),
        has_adduct = (FeatureNumber %in% matchSubset$FeatureNumber)
      )
    return(output)
  } 


UniqueCompounds <- #Parse unique compounds from output files for Dashboard Search
  function(filename) {
    
    input <- SampleInput(filename) 
    
    Compound_List <- input %>% select(Compound)
    
    if (grepl(filename,"$")) {
      var <- str_extract(filename,"(.*)(?=\\$_)") 
    } else {
      var <- str_extract(filename,"(?<=^)(.*)(?=_[0-9])")
    }
    
    Compound_working <- Compound_List %>%
      select(Compound, -contains("Flag")) %>%
      mutate(Compound_Flag = grepl("[a-z]",Compound)) %>%
      filter(Compound_Flag == TRUE) %>%
      mutate(Compound = sub("(\\s[0-9].*)","",Compound)) %>%
      unique()
    
    write_csv(Compound_working, paste0(var,"_UniqueCompounds.csv"))
    
  }

####### Specify Global Variables From File #######
source("config.R")

####### Specify Per-File Processing Workflow #######
source("ExperimentProcessing.R")

####### Specify Cross-Mode Combination Method (If Needed) #######
source("CrossExperimentMatching.R")

####### Process Files #######
processed_samples <- lapply(Filelist, NTA_Process)

####### Join Processed Files #######
if (PairedProcessing == TRUE) {
joined_samples <- NTA_Join(processed_samples)
}

### Placeholder compound list generation ###
if (PairedProcessing == TRUE) {
  lapply(joined_samples$rawfiles, UniqueCompounds)
} else {
  lapply(processed_samples, UniqueCompounds)
}

####### Eventually Run Dahsboard Search via API #######

####### Eventually Generate Toxpi Output In-Line #######
source("ToxPiScores.R")
 







