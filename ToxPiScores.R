
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

Dashboard_inputFile <- c("ChemistryDashboard-AdvancedSearch_2017-08-14_13-15-59.tsv")
Abundance_inputFile <- c("Joined_AHHS Soil_Filtered_2017-08-14.csv")

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
    gather(contains("Med"), key = File, value = MedAbun) %>%
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
    contains("Med")) %>%
  filter(grepl("[A-Z]",Compound)) %>%
  mutate(SampleN = rowSums(select(.,contains("_Med")) > 0, na.rm = TRUE )) %>%
  arrange(Compound,RetentionTime)

Combined_File <- full_join(Sample_ToxPi_Data,Dashboard_Tidy, by="Compound") %>%
  mutate(RowID = as.numeric(seq.int(nrow(.)))) %>%
  gather(contains("Med"), key="SampleID", value = "volume") %>%
  mutate(SampleID = paste0("log10_",SampleID),
         log_vol = log10(volume)) %>%
  select(-volume) %>%
  filter(log_vol != log10(0)) %>%
  spread(SampleID,log_vol) %>%
  mutate(log10_avgAbun = rowMeans(select(.,contains("Med")), na.rm = TRUE)) %>%
  select(-contains("Med")) %>%
  mutate(Chemical=str_c(Compound,PREFERRED_NAME, sep=": "),
         Chemical=str_wrap(Chemical, width=20),
         Pi_Bioactivity = (Tox21_Activities-min(Tox21_Activities, na.rm=TRUE))/(max(Tox21_Activities, na.rm=TRUE)-min(Tox21_Activities, na.rm=TRUE)),
         Pi_Exposure=(ExposureCategory-min(ExposureCategory,na.rm=TRUE))/(max(ExposureCategory, na.rm=TRUE)-min(ExposureCategory, na.rm=TRUE)),
         Pi_Abundance=(log10_avgAbun-min(log10_avgAbun, na.rm=TRUE))/(max(log10_avgAbun, na.rm=TRUE)-min(log10_avgAbun, na.rm=TRUE)),
         Pi_Frequency=(SampleN-min(SampleN, na.rm=TRUE))/(max(SampleN, na.rm=TRUE)-min(SampleN, na.rm=TRUE)),
         ToxPiScore=(Bio_Weight*Pi_Bioactivity+Expo_Weight*Pi_Exposure+Abun_Weight*Pi_Abundance+Freq_Weight*Pi_Frequency)) %>%
  select(IonMode,
         Chemical,
         contains("Name"),
         contains("Mass"),
         DTXSID,
         contains("Pi")) %>%
  arrange(Chemical)
  
Graph_File <- Combined_File %>%
    filter(ToxPiScore > 3) %>%
  select(Chemical,contains("Pi_")) %>%
  gather(contains("Pi_"), key = Pi_Categories, value = Score) %>%
  group_by(Chemical,Pi_Categories) %>%
  summarise_each(funs(mean)) %>%
  arrange(Chemical,Pi_Categories)
  #mutate(Pi_Categories = factor(Pi_Categories, levels = colnames(select(Combined_File,contains("Pi_"))))) 

New_Approach <- Combined_File %>%
  filter(ToxPiScore >3) %>%
  select(Chemical,contains("Pi"),-contains("mass"))

bar_w <- c(Abun_Weight,Bio_Weight,Expo_Weight,Freq_Weight)

#bar_w <- c(4/6,8/6,4/6,8/6)
         
bar_pos <- 0.5 * (cumsum(bar_w) + cumsum(c(0, bar_w[-length(bar_w)])))

ggplot(Graph_File[1:4,]) +
  geom_col(mapping=aes(x=bar_pos, y=Score, fill=Pi_Categories),
           width=bar_w,
           position=position_nudge()  ) +
  scale_x_continuous(breaks = bar_pos) +
  coord_polar() +
  theme_void() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values=c('dodgerblue', 'orange', 'green', 'magenta')) +
  ggtitle(Graph_File$Chemical) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-.075, 1) +
  theme(strip.text.x=element_text(size=12)) +
  labs(title=" \n Chemicals in Tire Crumbs with ToxPi Scores > 3")




