###### Descriptive Analysis ###############

# removing sci notation
options(scipen=999)
# requiring packages
library(tidyverse)
library(geofacet)
require(tmap)
library(geobr)
library(ggmap)
library(RColorBrewer)
library(xlsx)
library(patchwork)
library(gridExtra)
library(ggsn)
theme_set(theme_bw())
# loading data
data <- read_rds("00.data/data.rds")
grid <- read.csv("mg_grid.csv") %>% dplyr::select(-X)
data <- data %>% mutate(code = code %>% as.numeric()) %>% left_join(grid) 
# Filter for microregions with more than 10.000 cases in
filter <- data %>% group_by(code) %>%
  summarise(sum = sum(cases),code_meso) %>%
  filter (sum >= 10000) %>% arrange(code_meso) %>% unique
# Map minas ---------------------------------------------------------------
MG <- read_micro_region() %>% filter(abbrev_state %in% "MG")
MG_meso <- read_meso_region() %>% filter(abbrev_state %in% "MG")
MG <- MG %>% mutate(cat = ifelse(code_micro %in% filter$code, "Included", "Excluded")) %>%
  left_join(data[,c(1,28,29)] %>% unique, by  = c('code_micro'= 'code')) %>% arrange(code_meso) %>% 
  mutate(id = 1:66)
# Dengue ------------------------------------------------------------------
# sum of cases
summary(data)
data$cases %>% sum
# daily distribution by microregion
tab_cases <- data %>%
  group_by(microregion_name) %>% 
  summarize(
    Mean_pop=round(mean(pop),0),
    Sum_cases = sum(cases),
    annual_cases_mean = round(Sum_cases/10),
    Rate_annual = round(annual_cases_mean/Mean_pop*10^5,2),
    Min=min(cases),
    Percentil_25=quantile(cases,0.25),
    Mediana=quantile(cases,0.5),
    Mean=round(mean(cases),2),
    Percentil_75=quantile(cases,0.75),
    Max=round(max(cases),2),
    DP=round(sd(cases),2))
#Saving table
write.xlsx(tab_cases, "01.descriptive/tab_cases.xlsx")
# Exploring all cases -----------------------------------------------------
# time series
p1 <- data %>% 
  ggplot(aes(x = date, y = cases))+
  geom_line()+
  xlab('Data')+
  ylab("Casos")
  ggtitle('A')
# by month
p2 <- ggplot(data=data,aes(x=as.factor(month), y=cases)) +
  geom_bar(stat = "identity")+
  xlab ("Mês")+
  ylab('Casos')+
  ggtitle('B')
# by day of the week
p3 <- ggplot(data=data,aes(x=as.factor(weekday), y=cases)) +
  geom_bar(stat = "identity", fill = "black") +
  xlab ("Dia da Semana") +
  ylab("Casos")+
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
  ylab("Ano") + 
  xlab("Mês") +
  scale_fill_gradientn(name = "Taxa", colours=brewer.pal(9, "Greys"), trans = "log1p",
                       breaks = c(0, 10, 100, 300, 1000), labels = c(0, 10, 100, 300, 1000)) +
  scale_x_discrete(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  ggtitle('D')
# heat map based on year, month and rate by microregion  -------------------------
p5 <- data %>% mutate(year = year %>% as.factor) %>% 
  group_by(year, month, name, pop, code) %>% 
  summarise(cases = sum (cases),
            pop = max(pop)) %>% 
  mutate (rate= cases/pop*10^5)  %>%
  ggplot(aes(month, year, fill= rate)) +
  geom_raster()+
  ylab("Ano") + 
  xlab("Mês") +
  scale_fill_gradientn(name = "Taxa", colours=brewer.pal(9, "Greys"), trans = "log1p",
                       breaks = c(0, 10, 100, 300, 1000), labels = c(0, 10, 100, 300, 1000)) +
  scale_x_discrete(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  facet_geo(~name, grid = grid)+
  ggtitle('E')
# Layout for plot
layout <- "
AABBEEEE
CCDDEEEE
"
p1 + p2 + p3 + p4 + p5 +
  plot_layout(design = layout)
ggsave("03.figs/fig01.png", dpi = 600, width = 24, height = 15)

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
  temp <- temp_data$t_min
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
  
  big_table <- bind_rows(lil_table, big_table) %>% arrange(i)
}
write.xlsx(big_table, "01.descriptive/table_temp.xlsx")

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
