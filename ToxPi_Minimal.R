
Bio_Weight = 2
Abun_Weight = 1
Expo_Weight = 1
Freq_Weight = 2

Weight_tbl <- data.frame( sub_Pi = c("Bio_Pi", "Abun_Pi", "Expo_Pi", "Freq_Pi"),
                          weight = c(Bio_Weight, Abun_Weight, Expo_Weight, Freq_Weight), stringsAsFactors = FALSE)

Graph_Weights <- Weight_tbl

Graph_Weights$w <- cumsum(Graph_Weights$weight)
Graph_Weights$wm <- Graph_Weights$w - Graph_Weights$weight
Graph_Weights$wt <- with(Graph_Weights, wm + (w - wm)/2)

Pi_Test <- SampleInput("ToxPi_TestArea.xlsx")
  
Graph_Table <-Pi_Test %>%
         gather(contains("_Pi"), key = sub_Pi, value = height) %>%
  full_join(Weight_tbl) %>%
  group_by(Compound) %>%
  summarize(ToxPi = sum(height*weight)) %>%
  right_join(Pi_Test) %>%
  gather(contains("_Pi"), key = sub_Pi, value = height) %>%
  left_join(Graph_Weights)


library(ggplot2)
p  <- ggplot(Graph_Table, aes(ymin = 0))
p + geom_rect(aes(xmin = wm, xmax = w, ymax = height, fill = sub_Pi)) +
  coord_polar() +
  theme_bw() +
  facet_wrap(~Compound)