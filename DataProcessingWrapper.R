#James McCord
#November 26, 2017

##### Check for, install, and load necessary packages ####

check.packages <- function(package){
  new.package <- package[!(package %in% installed.packages()[, "Package"])]
  if (length(new.package)) 
    install.packages(new.package, dependencies = TRUE)
  sapply(package, require, character.only = TRUE)
}

packages<-c("tidyverse", "stringr", "sqldf", "readxl", "data.table", "tcltk2")

check.packages(packages)

####### Create global functions for use later ######
SampleInput <- # Generic data input from multiple file types
  function(filename){
    inputFilename <- c(filename)
    
    if (grepl("\\.(tsv|txt)$",inputFilename)) {
      Data <- read_tsv(inputFilename, comment = "#", na = c("NA"," ", "-", "0"), guess_max=100000)
    }
    
    if (grepl("\\.(csv)$",inputFilename)) {
      Data <- read_csv(inputFilename, comment = "#", na = c("NA"," ", "-", "0"), guess_max=100000)
    }
    
    if (grepl("\\.(xls|xlsx)$",inputFilename) | grepl(".xls$",inputFilename)) {
      Data <- read_excel(inputFilename, col_types = NULL, na = c("NA"," ", "-", ""), guess_max=100000)
    }
    
    return(Data)
  }

adductSearchSQL <- # RT and Mass matching with SQL
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
      mutate(Compound_Flag = grepl("[A-Z]",Compound)) %>%
      filter(Compound_Flag == TRUE) %>%
      mutate(Compound = sub("(\\s[0-9].*)","",Compound)) %>%
      unique()
    
    write_csv(Compound_working, paste0(input,"_UniqueCompounds.csv"))
    
  }

####### Specify Global Variables From File #######
source("config.R")

####### Eventually Specify Sample Input Method #####
source("SampleInput_MPP.R")

####### Specify Per-File Processing Workflow #######
source("ExperimentProcessing.R")

####### Specify Cross-Mode Combination Method (If Needed) #######
source("CrossExperimentMatching.R")

####### Batch Process Files #######
processed_samples <- lapply(Filelist, NTA_Process)

####### Join Processed Files in Pos/Neg Pairs #######
if (PairedProcessing == TRUE) {
joined_samples <- NTA_Join(processed_samples)
}

### Compound list generation for Dashboad Search ###
if (PairedProcessing == TRUE) {
  lapply(joined_samples$rawfiles, UniqueCompounds)
} else {
  lapply(processed_samples, UniqueCompounds)
}



####### Eventually Run Dahsboard Search via API #######


if (autoDashboard == TRUE){

  myList <- UniqueCompounds("AHHS Soil_$Esi-$_Filtered_2017-05-24.csv") %>%
    mutate(Compounds = paste0(Compound, "\n"))
  
  writeClipboard(head(myList$Compound))
  
check.packages("RSelenium")

selscroll <- function() webElem$sendKeysToElement(list(key = "page_down"))

rD <- rsDriver(remoteServerAddr = "localhost",
               browser = "chrome",
               check = FALSE,
               verbose = FALSE,
               )
remDr <- rD[["client"]]
remDr$setWindowSize(800, 800)
remDr$navigate("https://comptox.epa.gov/dashboard/dsstoxdb/batch_search")

webElem <- remDr$findElement("css", "body")

selscroll()

remDr$findElement('xpath', "//*[@name = 'input_types[]' and @value = 'compounds.mol_formula']")$clickElement()

textbox <- remDr$findElement(using= 'id', value = "list-search-text-box")
textbox$sendKeysToElement(list(key = 'control',"v", key = 'control'))

selscroll()

remDr$findElement('xpath', "//*[@href = '#batch-search-panel']")$clickElement()

selscroll()

#iframe <- remDr$findElement(using='id', value="list-search-select-box")
#remDr$switchToFrame(iframe)

option <- remDr$findElement(using = 'xpath', "//*/option[@value = 'tsv']")
option$clickElement()

CASRN <- remDr$findElement('xpath', "//*[@value = 'generic_substances.casrn']")
CASRN$clickElement()
remDr$findElement('xpath', "//*[@value = 'compounds.acd_iupac_name']")$clickElement()

molform <- remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'compounds.mol_formula']")
molform$clickElement()

remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'compounds.monoisotopic_mass']")$clickElement()
remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'opera_predictions']")$clickElement()
remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'test_predictions']")$clickElement()
remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'count(distinct b.fk_chemical_list_id) as data_sources']")$clickElement()
remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'hit_counts.active_assay_count as percent_active_calls, hit_counts.total_assay_count as number_active_assays_vs_total']")$clickElement()
remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'expocast_models.casrn as expocast, expocast_models.Total_median as expocast_median_exposure_prediction, expocast_models.inNHANES as nhanes']")$clickElement()
remDr$findElement('xpath', "//*[@name = 'columns[]' and  @value = 'toxval.toxval_type as toxval_data']")$clickElement()

webElem$sendKeysToElement(list(key = "end"))
remDr$findElement('xpath', "//*[@id = 'list-search-submit']")$clickElement()

remDr$close()
rD[["server"]]$stop()

}

####### Eventually Generate Toxpi Output In-Line #######
source("ToxPiScores.R")
 







