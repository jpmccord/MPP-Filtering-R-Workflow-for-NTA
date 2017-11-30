NTA_Join <- function(processed_samples) {
  
options(stringsAsFactors = FALSE)

var <- data.frame(unlist(processed_samples)) %>%
  transmute(filenames = unlist.processed_samples.) %>%
  mutate(ModeType = str_extract(filenames,"(?<=\\$)(.*)(?=\\$)")) %>%
  filter(!is.na(ModeType)) %>%
  mutate(rawfiles = sub("_\\$(.*)\\$_","_",filenames))

pairedvar <- sqldf(
  paste0(
    "select a.*, b.filenames as match_filenames,
    b.rawfiles as match_rawfiles,
    b.ModeType as match_ModeType
    from var as a, var as b  
    where a.rawfiles = b.rawfiles
    and a.ModeType != b.ModeType"
  )
  ) %>% filter(ModeType %in% c("Positive", "Pos", "Esi+"))

for (n in 1:(length(pairedvar$filenames))) {
  
  Positive_File <- SampleInput(pairedvar$filenames[1])
  
  Negative_File <- SampleInput(pairedvar$match_filenames[1])
  
  Pos1 <- Positive_File %>% 
    select(FeatureNumber,
           Compound,
           IonMode,
           Score,
           Mass,
           RetentionTime)
  
  Neg1 <- Negative_File %>% 
    select(FeatureNumber,
           Compound,
           IonMode,
           Score,
           Mass,
           RetentionTime)
  
  Cross_Mode <- 
    sqldf(
      paste0(
        "select a.*, b.FeatureNumber as match_FeatureNumber,
        b.IonMode as match_IonMode,
        b.Compound as match_Compound,
        b.Score as match_Score, 
        b.Mass as match_Mass,
        b.RetentionTime as match_RetentionTime
        from Pos1 as a, Neg1 as b  
        where a.Compound = b.Compound
        and abs(a.Mass-b.Mass)<=0.005 	
        and abs(a.RetentionTime-b.RetentionTime)<=1 
        order by a.Mass, b.Mass;"
      )
      ) %>%
    mutate(ppm_error = as.numeric((abs(match_Mass-Mass))/Mass *  10 ^ 6)) %>%
    filter(ppm_error <= duplicate_masserror)
  
  Pos_Unique_Flag <- Pos1 %>%
    select(FeatureNumber) %>%
    mutate(Pos_Unique_Flag = !(FeatureNumber %in% Cross_Mode$FeatureNumber))
  
  Neg_Unique_Flag <- Neg1 %>%
    select(FeatureNumber) %>%
    mutate(Neg_Unique_Flag = !(FeatureNumber %in% Cross_Mode$match_FeatureNumber))
  
  Positive_File <- inner_join(Positive_File,Pos_Unique_Flag)%>%
    select(-FeatureNumber)
  
  Negative_File <- inner_join(Negative_File,Neg_Unique_Flag)%>%
    select(-FeatureNumber)
  
  Full_Joined <- full_join(Positive_File,Negative_File) %>%
    mutate(BothModes_Flag = (Neg_Unique_Flag == FALSE | Pos_Unique_Flag == FALSE)) %>%
    mutate(BothModes_Flag = !is.na(BothModes_Flag)) %>%
    select(
      Compound,
      IonMode,
      Mass,
      RetentionTime,
      Score,
      contains("MeanAbun"),
      contains("Med"),
      contains("StdAbun"),
      contains("CV"),
      contains("adduct"),
      BothModes_Flag) %>%
    arrange(Compound,RetentionTime)
  
  Compound_Hits <- Full_Joined %>%
    select(Compound,IonMode,BothModes_Flag) %>%
    group_by(IonMode,BothModes_Flag) %>%
    count(Compound) %>%
    rename(N_Hits = n) %>%
    ungroup()%>%
    filter(BothModes_Flag == FALSE |
          (BothModes_Flag == TRUE & (IonMode %in% c("Positive", "Pos", "Esi+")))) %>%
    group_by(Compound)%>%
    summarize(N_CompoundHits = sum(N_Hits))
  
  Full_Joined <- full_join(Full_Joined,Compound_Hits)
  
  Full_Joined[Full_Joined == "0"] <- NA
  
  write_csv(Full_Joined, paste0("Joined_",pairedvar$rawfiles[n]))
}

joined_filelist <- pairedvar %>%
  mutate(rawfiles = paste0("Joined_",rawfiles))
          
return(joined_filelist)
}