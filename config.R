# Spectral Processing Source Parameters
# Specify input file(s) from MPP output: (txt, tsv, xls, or xlsx) separated by commas (i.e. c("File1.txt", "File2.txt") )

Filelist <- c(
    "data/AHHS Soil Negative Database Matched 05022017.txt",
    "data/AHHS Soil Positive Database Matched 05022017.txt"
  ) 

# One or more SampleDetail IDs to be used as reference blanks
Blank_Sample_Names <- c(
  "MB01",
  "MB02",
  "MB03"
)

# Output Variables
OutputFormat_wide = TRUE #Output format long or wide
OutputRaw = FALSE # Output unfiltered data in addition to filtered data
PairedProcessing = TRUE #If True, pair Pos and Neg outputs into one file for Dashboard searching
Collapse_Reps = TRUE # Collapse replicate injections into a single avg/med abundance

#Processing Workflow variables
ReplicateFeature_Threshold = 0.8  #Fraction of replicates where feature is observed to classify as "real"
ReplicateCV_Threshold = FALSE #False if no CV filtering, otherwise minium replicate CV to classify as "real"

SampleFeature_Threshold = 1 # Number of Unique Sample IDs hits to count as "real", NOT replicate injections
SampleThresholdFraction = FALSE #SampleFeature_Threshold as a fraction of samples analyzed rather than a raw number

BlankRatio_Threshold = 3 #Fold-change by which sample must exceed blank to be a sample feature
SignificantSample_Threshold = 1 # Number of samples where feature must be more abundant than the blank to classify as "real"

adduct_ppmsearch = TRUE #Units for adduct_masserror TRUE = ppm, FALSE = mDa
adduct_masserror = 20 #Mass error allowed for adducts
adduct_timetolerance = 0.05 #min, time difference allowed for adducts

duplicate_ppmsearch = TRUE #Units for duplicate_masserror TRUE = ppm, FALSE = mDa
duplicate_masserror = 10 #Mass error allowed for duplicates, mDa by default
duplicate_timetolerance =  0.05 #min, time difference allowed for features to be pooled as duplicates instead of unique isomers
#Duplicate settings also apply to joining hits across modes

LowScore_Threshold = 90 #Score threshold for formulae match to be kept
WeakScore_Threshold = 50 # Score threshold for formulas predicted to match poorly (i.e. contain F,Br,Cl,P)



 