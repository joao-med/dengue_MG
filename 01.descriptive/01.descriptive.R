###### Descriptive Analysis ###############

# removing sci notation
options(scipen=999)
# requiring packages
library(tidyverse)
library(raster)
library(geofacet)
require(tmap)
library(geobr)
library(ggmap)
library(RColorBrewer)
library(xlsx)
library(patchwork)
library(gridExtra)
theme_set(theme_bw())
setwd("final_codes/final_final_codes")
# loading data
data <- read_rds("00.data/data.rds")
grid <- read.csv("mg_grid.csv") %>% dplyr::select(-X)
data <- data %>% mutate(code = code %>% as.numeric()) %>% left_join(grid) 
# Filter for microregions with more than 10.000 cases in the span of time analysed
filter <- data %>% group_by(code) %>%
  summarise(sum = sum(cases)) %>%
  filter (sum >= 10000)
# Map minas ---------------------------------------------------------------
## Minas gerais geometry
MG <- read_micro_region() %>% filter(abbrev_state %in% "MG")

MG <- MG %>% mutate(cat = ifelse(code_micro %in% filter$code, "Included", "Excluded"),
                    id = 1:66)
p1 <- MG %>% 
  ggplot()+
  geom_sf(aes(fill = cat))+
  geom_sf_text(aes(label = id),color = 'white')+ 
  scale_fill_manual(values = c("grey", "black"), name = NULL)+
  scale_color_discrete(breaks = gtools::mixedsort(MG$name_state))+
  theme(legend.text = element_text(margin = margin(0, 0, 0, -24)))+
  xlab("")+
  ylab("")+
  theme_minimal()
p1+   tableGrob (cbind(MG$id[1:22],
                 MG$name_micro[1:22],
                 MG$id[23:44],
                 MG$name_micro[23:44],
                 MG$id[45:66],
                 MG$name_micro[45:66]))
ggsave("03.figs/fig01.png", height = 7, width = 16)
# Dengue ------------------------------------------------------------------
# sum of cases
summary(data)
data$cases %>% sum
# daily distribution by microregion
tab_cases <- data %>%
  group_by(microregion_name) %>% 
  summarize(Min=min(cases),
            Q1=quantile(cases,0.25),
            Q2=quantile(cases,0.5),
            Q3=quantile(cases,0.75),
            Max=round(max(cases),2),
            Mean=round(mean(cases),2),
            DP=round(sd(cases),2),
            Sum_cases = sum(cases),
            Mean_pop=round(mean(pop),0),
            Rate = round(Sum_cases/Mean_pop*10^5,2),
            annual_cases_mean = round(Sum_cases/10))
#Saving table
write.csv(tab_cases, "01.descriptive/tab_cases.csv")

# Exploring all cases -----------------------------------------------------
# time series
p1 <- data %>% 
  ggplot(aes(x = date, y = cases))+
  geom_line()+
  ggtitle('A')
# by month
p2 <- ggplot(data=data,aes(x=as.factor(month), y=cases)) +
  geom_bar(stat = "identity")+
  xlab ("month")+
  ggtitle('B')
# by day of the week
p3 <- ggplot(data=data,aes(x=as.factor(weekday), y=cases)) +
  geom_bar(stat = "identity") +
  xlab ("weekday") +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  ggtitle('C')
#Sunday and Saturday are the days with the least cases
# Heatmap of rate by month and year
p4 <- data %>% mutate(year = year %>% as.factor) %>% 
  group_by(year, month, microregion_name) %>%
  summarise(cases = sum (cases),
            pop = max(pop)) %>% 
  summarise(pop = sum(pop),
            cases =sum(cases)) %>%
  mutate (rate = cases/pop*10^5) %>%
  ggplot(aes(month, year, fill= rate)) +
  geom_raster()+
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Rate", colours=brewer.pal(9, "Greys"), trans = "log1p",
                       breaks = c(0, 10, 100, 300, 1000), labels = c(0, 10, 100, 300, 1000)) +
  scale_x_discrete(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  ggtitle('D')

(p1 | p2 )/ (p3 | p4)
ggsave("03.figs/fig02.png", dpi = 200, width = 8, height = 8)
# heat map based on year, month and rate by microregion  ---------------------------------------------
data %>% mutate(year = year %>% as.factor) %>% 
  group_by(year, month, name, pop, code) %>% 
  summarise(cases = sum (cases),
            pop = max(pop)) %>% 
  mutate (rate= cases/pop*10^5)  %>%
  ggplot(aes(month, year, fill= rate)) +
  geom_raster()+
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Rate", colours=brewer.pal(9, "Greys"), trans = "log1p",
                       breaks = c(0, 10, 100, 300, 1000), labels = c(0, 10, 100, 300, 1000)) +
  scale_x_discrete(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  facet_geo(~name, grid = grid)
ggsave("03.figs/fig03.png", dpi = 200, width = 12, height = 13)

# Temperature -------------------------------------------------------------

## Exploring weather data Minas gerais - ERA5
# matching classes of code
data <- data %>% mutate(code = code %>% as.integer())
# Selecting mean variables by year
yearly_ERA_mg <- data %>%
  group_by(microregion_name,code) %>% 
  summarise(t_min = t_min %>% mean(),
            t_max = t_max %>% mean(),
            t_mean = t_mean %>% mean())

# merging polygons with weather
yearly_ERA_mg <- MG %>% 
  dplyr::select(code = code_micro, geom) %>% 
  left_join(yearly_ERA_mg)

# Minimum Temperature daily distribution by microregion
big_table <- tibble()
for (i in data$microregion_name  %>% unique) {
  temp_data <- data %>% filter(microregion_name  == i)
  temp <- data$t_min
  cities <- as.vector(data$microregion_name  %>% unique)
  
  # temperature quantis
  summary_temp <- round(summary(temp_data$t_min),1) %>%
    as.vector() %>%
    as.tibble %>%
    t()
  colnames(summary_temp) <-  
    c("Min", "Percentile_25", "Median", "Mean", "Percentile_75", "Max")
  summary_temp <- summary_temp %>% as.tibble
  
  P0 <-  round(quantile(temp,probs=0),1) %>% tibble("0%" = .)
  P1 <-  round(quantile(temp,probs=0.01),1) %>% tibble("1%" = .)
  P2.5 = round(quantile(temp,probs=0.025),1) %>% tibble("2.5%" = .)
  P10 = round(quantile(temp,probs=0.1),1) %>% tibble("10%" = .)
  P25 = round(quantile(temp,probs=0.25),1) %>% tibble("25%" = .)
  P50 = round(quantile(temp,probs=0.5),1) %>% tibble("50%" = .)
  P75 = round(quantile(temp,probs=0.75),1) %>% tibble("75%" = .)
  P90 = round(quantile(temp,probs=0.9),1) %>% tibble("90%" = .)
  P97.5 = round(quantile(temp,probs=0.975),1) %>% tibble("97.5%" = .)
  P99 = round(quantile(temp,probs=0.99),1) %>% tibble("99%" = .)
  P100 = round(quantile(temp,probs=1),1) %>% tibble("100%" = .)
  
  lil_table <- tibble(i,
                      Minimum = summary_temp$Min,
                      P2.5,P10,P25,
                      Median = summary_temp$Median,
                      Mean = summary_temp$Mean,
                      P50,P75,P90,P97.5,
                      Max = summary_temp$Max)
  
  big_table <- bind_rows(lil_table, big_table)
}
write.xlsx(big_table, "01.descriptive/table_temp.xlsx")
####  maps of average min temperature by microregion
# Yearly min temp
p1 <- yearly_ERA_mg %>% ggplot(aes(fill = t_min)) +
  geom_sf ()+
  scale_fill_gradientn(name = "TºC min", colours=brewer.pal(9, "Greys"),
                       breaks = c(14, 16, 18, 20), labels = c('14ºC', '16ºC', '18ºC', '20ºC')) +
  ggtitle("A")
# create a list with regional datasets
regions <- data$microregion_name %>% unique()
data_region <- lapply(regions,function(x) data[data$microregion_name==x,])
names(data_region) <- regions
n <- length(regions)
# temperature ranges for each region
ranges_data <- t(sapply(data_region, function(x) range(x$t_min,na.rm=T)))
ranges_data <- data.frame(ranges_data)
ranges_data$indice <- 1:length(regions)
ranges_data$region <- row.names(ranges_data)
ranges_data$temperature <- (ranges_data$X1+ranges_data$X2)/2 
ranges_data <-  ranges_data %>% as.tibble()

p2 <- ranges_data %>% 
  ggplot(aes(x = reorder(region,temperature), y = temperature))+
  geom_point()+
  geom_errorbar(aes(ymin= X1, ymax = X2),
                width = 0.8)+
  ylab ("Temperature (ºC)")+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  ggtitle("B")

p1/p2
ggsave('03.figs/fig04.png', height = 11)

