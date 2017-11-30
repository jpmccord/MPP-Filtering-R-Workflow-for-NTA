#SampleInput No Formatting

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

PrepareFile() <- function(inputFilename) {
SampleInput(inputFilename)
}