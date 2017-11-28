NTA_Process <- function(inputFilename)
{
  Data_Input <- SampleInput(inputFilename) %>%
  arrange(Mass, `Retention Time`) %>%
  mutate(FeatureNumber = as.numeric(seq.int(nrow(.))))

Abundance_Data <- Data_Input %>%
  gather(contains("_"), key = "filename", value = "volume") %>%
  separate(
    "filename", # Split filenames at "_" to gain metadata
    into = c("ExpID", "IonMode", "ModeDetail", "SampleDetail", "Rep"),
    sep = "_"
  ) %>%
  mutate(
    volume = ifelse(volume == 1 | volume == 0, NA, volume), #Replace missing value export with missing
    SampleDetail = ifelse(SampleDetail %in% Blank_Sample_Names, "Blank", SampleDetail) #Rename Blank sample ID internally
  ) %>%
  select(FeatureNumber,ExpID,IonMode,ModeDetail,SampleDetail,Rep,volume) #Drop Feature Data from this table for speed

ExpID <- Abundance_Data$ExpID[1]                # Store ExpID and list of SampleIDs for later
theMode <- Abundance_Data$IonMode[1]
SampleIDs <- unique(Abundance_Data$SampleDetail)

print(paste0("Generating Feature Date for ", ExpID,"_",theMode))

Feature_Data <- Data_Input %>%
  select(-contains(ExpID)) %>% #Drop Abundance Data from this matrix for speed 
  mutate(Score = as.numeric(str_extract(
    Annotations,"(?<=overall=)(.*)(?=,)"))) %>%
  separate(Compound, c("Compound", "Delete"), sep = " Esi") %>%
  rename("RetentionTime" = `Retention Time`) %>%
  select(FeatureNumber,Mass,RetentionTime,Compound,Score,Frequency) %>%
  left_join(distinct(select(Abundance_Data,IonMode,FeatureNumber))) %>%
  distinct(.keep_all=TRUE)


####### Summarize data. ####

print(paste0("Calculating Summary Statistics for ", ExpID,"_",theMode))

Abundance_Data.table <- data.table(Abundance_Data) # Performing Calculations in data.table package for efficiency

setkey(Abundance_Data.table, FeatureNumber, SampleDetail)

Graph_Table.table <- Abundance_Data.table[,list(NAbun = sum(!is.na(volume)),
                                                MedAbun = round(median(volume, na.rm =TRUE),0),
                                                MeanAbun = round(mean(volume, na.rm = TRUE),0),
                                                StdAbun = round(sd(volume, na.rm = TRUE),2),
                                                Count = length(volume)),
                                          by=list(FeatureNumber, SampleDetail)]

Graph_Table <- tbl_df(Graph_Table.table[, `:=` (CV = round((StdAbun / MeanAbun), 2),
                                                ReplicateThresholdFlag = ((NAbun/Count) >= ReplicateFeature_Threshold),
                                                CVThresholdFlag = 
                                                ((is.numeric(ReplicateCV_Threshold) &
                                                (StdAbun/MeanAbun < ReplicateCV_Threshold)) |
                                                ReplicateCV_Threshold == FALSE))
                                        ]) %>%
  mutate(SampleFeature = ifelse(ReplicateThresholdFlag & CVThresholdFlag, 1, 0))

BlankSub_Data <- Graph_Table %>% 
  filter(SampleDetail == "Blank") %>%
  mutate(Blank_avgvol = ifelse(is.nan(MeanAbun), 0, MeanAbun),
         Blank_medvol = ifelse(is.na(MedAbun), 0 , MedAbun)) %>%
  select(FeatureNumber, ends_with("vol")) %>%
  full_join(filter(Graph_Table,SampleDetail != "Blank")) %>%
  mutate(subMeanAbun = ifelse(MeanAbun - Blank_avgvol < 0, 0, MeanAbun - Blank_avgvol),
         subMedAbun = ifelse(MedAbun - Blank_medvol < 0, 0, MedAbun - Blank_medvol)) 

# Create flag for features that are not considered sample features.
Remove_by_Sample_Feature <- BlankSub_Data %>%
  group_by(FeatureNumber) %>%
  summarise(PresentInSamples = sum(SampleFeature, na.rm = TRUE)) %>%
  mutate(staticFilter = (PresentInSamples < SampleFeature_Threshold),
         thresholdFilter = (PresentInSamples/length(SampleIDs) <= SampleFeature_Threshold))

if (SampleThresholdFraction == TRUE) {
  Remove_by_Sample_Feature <- mutate(Remove_by_Sample_Feature, FilterMe = thresholdFilter)  %>%
    filter(FilterMe == TRUE)
} else {
  Remove_by_Sample_Feature <- mutate(Remove_by_Sample_Feature, FilterMe = staticFilter) %>%
    filter(FilterMe == TRUE)
}

SampleFeature_Flag <- Feature_Data %>%
  select(FeatureNumber) %>%
  mutate(
    FalseSampleFeature_Flag = FeatureNumber %in% Remove_by_Sample_Feature$FeatureNumber
  )

# Calculate blank ratio and create flag for features that are below the spike threshold.
SpikeThreshold_Flag <- BlankSub_Data %>%
mutate(BlankRatio = MedAbun / (Blank_medvol + 1),
       AboveSpikeThreshold = (BlankRatio >= BlankRatio_Threshold)) %>%
  select(FeatureNumber, MedAbun, AboveSpikeThreshold) %>%
  group_by(FeatureNumber) %>%
  summarise(SumSpikeThreshold = sum(AboveSpikeThreshold, na.rm = TRUE)) %>%
  select(FeatureNumber, SumSpikeThreshold) %>%
  arrange(FeatureNumber) %>%
  ungroup() %>%
  mutate(BelowSpike_Flag = (SumSpikeThreshold < SignificantSample_Threshold)) %>%
  select(FeatureNumber, BelowSpike_Flag)

#Create negative mass defect variable, and flags for MissingScore, LowScore, and expected WeakMatches
Score_Mass_Flags <- Feature_Data %>%
  select(Score, Mass, FeatureNumber, Compound) %>%
  mutate(
    RoundedMass = round(Mass, 0),
    MassDefect = Mass - RoundedMass,
    NegMassDefect_Flag = (MassDefect < 0.05),
    MissingScore_Flag = (is.na(Score)),
    LowScore_Flag = (Score < LowScore_Threshold | is.na(Score)),
    WeakMatch_Flag = ifelse((NegMassDefect_Flag ==TRUE &
                               (grepl("F[0-9]", Compound) |
                                  grepl("Br[0-9]", Compound) |
                                  grepl("Cl[0-9]", Compound) |
                                  grepl("P[0-9]", Compound))),
                            (Score > WeakScore_Threshold), FALSE)  
    # Allow poor matches for compounds predicted to have poor algorithm matches 
  ) %>%
  select(FeatureNumber, Compound, Score, NegMassDefect_Flag, MissingScore_Flag, LowScore_Flag, WeakMatch_Flag)



####### Create Flags for potential adducts ######
adductFlag <- function(input,adductName,adductMass,masserror=adduct_masserror,timetolerance = adduct_timetolerance){
  
  name_is_adduct <- (paste0("is_",adductName,"_adduct"))
  name_has_adduct <- (paste0("has_",adductName,"_adduct"))
  
  df <- adductSearchSQL(input,
                        adductMass = adductMass,
                        masserror = masserror,
                        timetolerance = timetolerance,
                        ppm = adduct_ppmsearch) %>% 
    rename_(.dots = setNames("is_adduct", name_is_adduct))%>%
    rename_(.dots = setNames("has_adduct", name_has_adduct)) %>%
    left_join(input)
  return(df)
}

print(paste0("Searching Adducts for ", ExpID,"_",theMode))

if (theMode == "Negative") {
  Adduct_Flags <- Feature_Data %>%
    select(FeatureNumber,Mass,RetentionTime,IonMode) %>%
    adductFlag(adductName = "Formate", adductMass = 44.0409) %>%
    select(FeatureNumber,contains("adduct"))
}

if (theMode == "Positive") {
  Adduct_Flags <- Feature_Data %>%
    select(FeatureNumber,Mass,RetentionTime,IonMode) %>%
    adductFlag(adductName = "Na", adductMass = 21.98194) %>%
    adductFlag(adductName = "Ammonium", adductMass = 17.02655) %>%
    select(FeatureNumber,contains("adduct"))
}

####### Flag and/or remove potential duplicates #####
print(paste0("Finding Duplicates in ", ExpID,"_",theMode))

duplicateSQL <- adductSearchSQL(
  Feature_Data,
  adductMass = 0,
  masserror = duplicate_masserror,
  timetolerance = duplicate_timetolerance,
  ppm = duplicate_ppmsearch) %>%
  filter(is_adduct==TRUE | has_adduct == TRUE) %>%
  arrange(Mass)

UniquenessFilter <- Feature_Data %>%
  mutate(
    Duplicate_Flag = FeatureNumber %in% duplicateSQL$FeatureNumber
  ) %>%
  filter(Duplicate_Flag == TRUE) %>%
  separate(Compound, c("Compound", "Delete"), sep = "@") %>%
  select(-Delete) %>%
  separate(Compound, c("Compound", "Delete"), sep = "\\s") %>%
  select(-contains("adduct"), -Delete) %>%
  mutate(RetentionTime = round(RetentionTime,1)) %>% # Use RT grouping at 0.1 min resolution, I'm Sorry
  mutate(
    holdme = ifelse(grepl("^[0-9]", Compound), round(as.numeric(Compound),1), 0),
    Roughname = ifelse(holdme == 0, Compound, holdme),
    newScore = ifelse(is.na(Score), 0, Score)
  ) %>%
  mutate(Roughname = ifelse(holdme == 0, Compound, holdme)) %>%
  group_by(Roughname, RetentionTime) %>%
  top_n(1, newScore) %>%
  top_n(1, Frequency) %>%
  top_n(-1, FeatureNumber) %>% # Select a single Feature from clusters prioritized based on Score, # of Observations, then lowest FeatureNumber
  ungroup() %>%
  select(FeatureNumber)

Duplicate_Flags <- Feature_Data %>%
  select(FeatureNumber) %>%
  mutate(
    Duplicate_Flag = FeatureNumber %in% duplicateSQL$FeatureNumber,
    Duplicate_Keep_Flag = FeatureNumber %in% UniquenessFilter$FeatureNumber
  )

All_Flags <- inner_join(SampleFeature_Flag,Adduct_Flags) %>%
  inner_join(SpikeThreshold_Flag)%>%
  inner_join(Score_Mass_Flags)%>%
  inner_join(Duplicate_Flags)

All_Data <- left_join(left_join(BlankSub_Data, Feature_Data), All_Flags) %>%
  select(FeatureNumber, 
         SampleDetail,
         MedAbun,
         MeanAbun,
         subMeanAbun,
         subMedAbun,
         Mass,
         RetentionTime,
         Compound,
         Score,
         IonMode,
         contains("Flag"),
         contains("adduct")) %>%
  mutate(
    Compound = ifelse(
      (LowScore_Flag == FALSE | (LowScore_Flag == TRUE & WeakMatch_Flag == TRUE)),
      Compound,
      paste0(Mass,"@",RetentionTime)
    )
  )

if (Collapse_Reps == FALSE) {
  All_Data <- Abundance_Data %>%
    left_join(unique(select(BlankSub_Data,FeatureNumber,Blank_avgvol,Blank_medvol))) %>%
    filter(SampleDetail != "Blank") %>%
    mutate(SampleDetail = gsub(" ",".",SampleDetail))%>%
    mutate(subMeanAbun = volume - Blank_avgvol,
           subMedAbun = volume - Blank_medvol,
           MeanAbun = volume,
           MedAbun = volume,
           Rep = str_extract(Rep,"[0-9]*"),
           SampleDetail = paste0(SampleDetail,"-",Rep)) 
}

Raw_Flagged <- All_Data %>%
  mutate(MeanAbun = ifelse(is.nan(MeanAbun), NA, MeanAbun)) %>%
  left_join(Feature_Data) %>%
  left_join(All_Flags)

if (OutputFormat_wide == TRUE | PairedProcessing == TRUE) {
  Wide_Table <- Raw_Flagged %>%
    unite(new,
          MeanAbun,
          MedAbun,
          subMeanAbun,
          subMedAbun,
          sep = "_") %>%
    select(FeatureNumber,new,SampleDetail)
  
  Wide_Table$SampleDetail <-
    paste0("joined_", Wide_Table$SampleDetail) #Append a label for later
  
  Spread_Table <- arrange(distinct(Wide_Table), FeatureNumber) %>%
    spread(key = SampleDetail, value = new) #Spread into columns
  
  DummyFrame <-
    Spread_Table %>% select(FeatureNumber, starts_with("joined")) #Grab joined column subset from before
  
  SampleColumns <-
    colnames(select(DummyFrame, FeatureNumber, starts_with("joined_"))) #Make a list of joined sample names
  
  ColumnSplit <-
    function(incol, indata) {
      # Define a custom function to apply separate_ across variables for defined I/O formats
      babel <- sub(".*_", "", incol)        # Grab each sample ID
      output <- select(separate_(
        indata,
        incol,
        c(
          #Prep columns for each variable explicitly
          paste0(babel, "_MeanAbun"),
          paste0(babel, "_MedAbun"),
          paste0(babel, "_subMeanAbun"),
          paste0(babel, "_subMedAbun")),
        sep = "_"
      ),
      starts_with(babel))
      
      return(output)
    }
  
  DummyFrame_Split <-
    data.frame(lapply(SampleColumns, FUN = ColumnSplit, indata = DummyFrame)) %>% # This separates sample columns with new names
    rename(FeatureNumber = FeatureNumber_MeanAbun) %>%
    select(-starts_with("FeatureNumber_"))
  
  #Coerce variables back to proper types
  DummyFrame_Split[DummyFrame_Split == "NA"] <- NA
  DummyFrame_Split[DummyFrame_Split == "NaN"] <- NA
  DummyFrame_Split <- mutate_all(DummyFrame_Split, funs(as.numeric))
  
  
  #Join Features with split new abundance values
  Sample_Raws <- inner_join(
    select(Feature_Data,
           FeatureNumber,
           Compound,
           IonMode,
           Mass,
           RetentionTime,
           Score),
           DummyFrame_Split)
  
  Raw_Flagged_Wide <- left_join(Sample_Raws, All_Flags) %>%
    mutate(
      Compound = ifelse(
        (LowScore_Flag == FALSE | (LowScore_Flag == TRUE & WeakMatch_Flag == TRUE)),
        Compound,
        paste0(Mass,"@",RetentionTime)
      )
    )
  
  if (OutputFormat_wide == TRUE) {
    Raw_Flagged <- Raw_Flagged_Wide
  }
}

Filtered_Data <- Raw_Flagged %>%
  filter(
    FalseSampleFeature_Flag == FALSE,
    BelowSpike_Flag == FALSE,
    Duplicate_Flag == FALSE | (Duplicate_Flag == TRUE & Duplicate_Keep_Flag == TRUE)
  )  %>%
  mutate(
    Compound = ifelse(
      (LowScore_Flag == FALSE | (LowScore_Flag == TRUE & WeakMatch_Flag == TRUE)),
      Compound,
      paste0(Mass,"@",RetentionTime)
    )
  )
#Export file that includes all features and all flags, no additional filtering

outputName <- paste0(ExpID,"_$",Raw_Flagged$IonMode[1],"$_")

print("Writing Output Files")

returnName <- paste0(outputName,"Filtered_", Sys.Date(),".csv")

if (OutputRaw == TRUE) {
write_csv(Raw_Flagged, paste0(outputName,"RawFlagged_", Sys.Date(),".csv"), na = "")
}

#Export file that includes features where the blanks have been subtracted and some features have been filtered out along with the flag columns.
write_csv(Filtered_Data, paste0(outputName,"Filtered_", Sys.Date(),".csv"), na = "")

if (PairedProcessing == TRUE & OutputFormat_wide == FALSE) {
  Filtered_Data <- Raw_Flagged_Wide %>%
    filter(
      FalseSampleFeature_Flag == FALSE,
      BelowSpike_Flag == FALSE,
      Duplicate_Flag == FALSE | (Duplicate_Flag == TRUE & Duplicate_Keep_Flag == TRUE)
    )
  write_csv(Filtered_Data, paste0(outputName,"W+Filtered_", Sys.Date(),".csv"), na = "")
  returnName <- paste0(outputName,"W+Filtered_", Sys.Date(),".csv")
}

print(paste0("Files Exported Successfully"))

return(returnName)
}