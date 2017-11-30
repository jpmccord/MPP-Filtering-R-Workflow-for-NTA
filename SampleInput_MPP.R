#Sample Input MPP

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

PrepareFile <- function(inputFilename) {
Data_Input <- SampleInput(inputFilename) %>%
  arrange(Mass, `Retention Time`) %>%
  mutate(FeatureNumber = as.numeric(seq.int(nrow(.)))) %>%
  gather(contains("_"), key = "filename", value = "volume") %>%
  separate(
    "filename", # Split filenames at "_" to gain metadata
    into = c("ExpID", "IonMode", "ModeDetail", "SampleDetail", "Rep"),
    sep = "_"
  ) %>%
  mutate(
    volume = ifelse(volume == 1, NA, volume), #Replace missing value export with missing
    SampleDetail = ifelse(SampleDetail %in% Blank_Sample_Names, "Blank", SampleDetail) #Rename Blank sample ID internally
  ) %>%
  mutate(Score = as.numeric(str_extract(
    Annotations,"(?<=overall=)(.*)(?=,)"))) %>%
  separate(Compound, c("Compound", "Delete"), sep = " Esi") %>%
  select(-Delete) %>%
  rename("RetentionTime" = `Retention Time`)

return(Data_Input)
}
