
formatToxPi <- function(inputFile){
Dashboard_Tidy <- read_tsv(inputFile, na="-") %>%
  setnames(old = colnames(.), new = gsub("\\.","_",make.names(colnames(.)))) %>%
  mutate(Category = ifelse(is.na(TOXCAST___ACTIVE) | is.na(EXPOCAST_MEDIAN_EXPOSURE_PREDICTION_MG_KG_BW_DAY), "B", "A")) %>%
  arrange(INPUT,desc(DATA_SOURCES)) %>%
  group_by(INPUT) %>%
  mutate(max_Data_Sources = max(as.numeric(DATA_SOURCES))) %>%
  ungroup() %>%
  mutate(DataSourceRatio = (as.numeric(DATA_SOURCES) / max_Data_Sources),
         Category_Number = ifelse((DataSourceRatio == 1), 1, 2),
         Compound_Classification = paste0(Category,Category_Number)) %>%
  select(-max_Data_Sources) %>%
  mutate(ExposureCategory = (9 + ceiling(log10(EXPOCAST_MEDIAN_EXPOSURE_PREDICTION_MG_KG_BW_DAY)))) %>%
  mutate(Tox21_Activities = TOXCAST___ACTIVE/100) 
 return(Dashboard_Tidy)
}

Bio_Weight = 2
Abun_Weight = 1
Expo_Weight = 1
Freq_Weight = 2

Dashboard_inputFile <- c("data/Dashboard Search Output.tsv")
Abundance_inputFile <- c("AHHS Soil_$Negative$_Filtered_2017-11-30.csv") 

Dashboard_Tidy <- formatToxPi(Dashboard_inputFile) %>%
  rename(Compound = INPUT)

Sample_Data <- SampleInput(Abundance_inputFile)

N_abuns <- Sample_Data %>%
  select(
    IonMode,
    Compound,
    Mass,
    RetentionTime,
    Score,
    contains("Med")) %>%
    filter(grepl("[A-Z]",Compound)) %>%
    gather(contains("subMed"), key = File, value = MedAbun) %>%
    mutate(present = ifelse(is.na(MedAbun),0,1)) %>%
  group_by(Compound,Mass) %>%
  summarize(SampleN = sum(present, na.rm = TRUE))
  
Sample_ToxPi_Data <- Sample_Data %>%
  select(
    IonMode,
    Compound,
    Mass,
    RetentionTime,
    Score,
    contains("subMed")) %>%
  filter(grepl("[A-Z]",Compound)) %>%
  left_join(N_abuns) %>%
  arrange(Compound,RetentionTime)%>%
  ungroup()

Combined_File <- full_join(Sample_ToxPi_Data,Dashboard_Tidy) %>%
  mutate(RowID = seq.int(nrow(.))) %>%
  gather(contains("subMed"), key="SampleID", value = "volume") %>%
  mutate(SampleID = paste0("log10_",SampleID),
         log_vol = log10(volume)) %>%
  select(-volume) %>%
  filter(log_vol != log10(0)) %>%
  spread(SampleID,log_vol) %>%
  mutate(log10_avgAbun = rowMeans(select(.,contains("Med")), na.rm = TRUE)) %>%
  select(-contains("Med")) %>%
  mutate(Chemical=str_c(Compound,Preferred_Name, sep=": "),
         Pi_Bioactivity = (Tox21_Activities-min(Tox21_Activities, na.rm=TRUE))/(max(Tox21_Activities, na.rm=TRUE)-min(Tox21_Activities, na.rm=TRUE)),
         Pi_Exposure=(ExposureCategory-min(ExposureCategory,na.rm=TRUE))/(max(ExposureCategory, na.rm=TRUE)-min(ExposureCategory, na.rm=TRUE)),
         Pi_Abundance=(log10_avgAbun-min(log10_avgAbun, na.rm=TRUE))/(max(log10_avgAbun, na.rm=TRUE)-min(log10_avgAbun, na.rm=TRUE)),
         Pi_Frequency=(SampleN-min(SampleN, na.rm=TRUE))/(max(SampleN, na.rm=TRUE)-min(SampleN, na.rm=TRUE)),
         ToxPiScore=(Bio_Weight*Pi_Bioactivity+Expo_Weight*Pi_Exposure+Abun_Weight*Pi_Abundance+Freq_Weight*Pi_Frequency)) %>%
  select(IonMode,
         Compound,
         Chemical,
         contains("Name"),
         contains("Mass"),
         DTXSID,
         contains("Pi"),
         Category,
         Compound_Classification) %>%
  arrange(Chemical) %>%
  group_by(Chemical) %>%
  filter(ToxPiScore == max(ToxPiScore)) %>%
  ungroup() %>% group_by(Compound) %>%
  mutate(alpha_beta = ifelse(ToxPiScore == max(ToxPiScore), "\U03B1","\U03B2")) %>%
  mutate(Class = paste0(Compound_Classification, alpha_beta)) 

write.csv(Combined_File,"ToxPiTable.csv")


Bio_Weight = 2
Abun_Weight = 1
Expo_Weight = 1
Freq_Weight = 2

Weight_tbl <- data.frame( sub_Pi = c("Pi_Bioactivity", "Pi_Abundance", "Pi_Exposure", "Pi_Frequency"),
                          weight = c(Bio_Weight, Abun_Weight, Expo_Weight, Freq_Weight), stringsAsFactors = FALSE)

Graph_Weights <- Weight_tbl

Graph_Weights$w <- cumsum(Graph_Weights$weight)
Graph_Weights$wm <- Graph_Weights$w - Graph_Weights$weight
Graph_Weights$wt <- with(Graph_Weights, wm + (w - wm)/2)

Pi_Test <- Combined_File

Graph_Table <- Pi_Test %>%
  gather(contains("Pi_"), key = sub_Pi, value = height) %>%
  full_join(Weight_tbl) %>%
  group_by(Chemical) %>%
  summarize(ToxPi = sum(height*weight)) %>%
  right_join(Pi_Test) %>%
  gather(contains("Pi_"), key = sub_Pi, value = height) %>%
  left_join(Graph_Weights) %>%
  ungroup() %>%
  mutate(Chemical=str_wrap(Chemical, width=20))%>%
  arrange(desc(ToxPiScore))

Chem_List <- unique(Graph_Table$Chemical)

paginate = 9
chunks <- split(Chem_List, ceiling(seq_along(Chem_List)/paginate))

pdf("ToxPiGraphs.pdf")
for (i in 1:length(chunks)){
  graphable <- filter(Graph_Table,
                      Chemical %in% chunks[[i]])
  
  p <- ggplot(graphable, aes(ymin = 0))
  
  p1 <- p + geom_rect(aes(xmin = wm, xmax = w, ymax = height, fill = sub_Pi)) +
    coord_polar() +
    theme_bw() +
    theme(axis.text.x = element_blank())+
    ylim(0,1) +
    facet_wrap(~Chemical)
  
  print(p1)
}
dev.off()

#Random Subsample
#p  <- ggplot(filter(Graph_Table, Chemical %in% Graph_Table$Chemical[sample(1:nrow(Graph_Table), 6)]), aes(ymin = 0))
#p + geom_rect(aes(xmin = wm, xmax = w, ymax = height, fill = sub_Pi)) +
#  coord_polar() +
#  theme_bw() +
#  theme(axis.text.x = element_blank())+
#  ylim(0,1) +
#  facet_wrap(~Chemical)


# For Generically scaled graphics
#p + geom_rect(aes(xmin = wm, xmax = w, ymax = height, fill = sub_Pi)) +
#  coord_polar() +
#  theme_bw() +
#  theme(axis.text.x = element_blank())+
#  facet_wrap(~Chemical)
  




