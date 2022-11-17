## Second stage -- Meta-Analysis

## Removing sci notation
options(scipen=999)
# Loading libraries
library(dlnm)
library(mvmeta) 
library(splines)
library(tidyverse)
library(png)

# loading data
data_complete <- read_rds("00.data/data.rds")
# Filtering microregions with populations more than 10.000 cases in the span of time analysed
filter <- data_complete %>% group_by(microregion_name) %>% 
  summarise(sum = sum(cases)) %>%
  filter (sum >= 10000)

data_complete <- data_complete %>% filter(microregion_name %in% filter$microregion_name)
# Creating cross-basis
cb_t_min <- crossbasis(data_complete$t_min, lag = 22, 
                       argvar=list(fun="ns",df=3),
                       arglag=list(fun="poly",degree=3))
# create vector with regions names
regions <- as.character(unique(data_complete$microregion_name)) 

# create a list with regional datasets
data <- lapply(regions,function(x) data_complete[data_complete$microregion_name==x,])
names(data) <- regions
n <- length(regions)

# TEMPERATURE DISTRIBUTION
# percentiles 10 e 90 of temperature distribution
P2.5 <- round(quantile(data_complete$t_min, probs = 0.025),1)
P2.5
P10 <- round(quantile(data_complete$t_min, probs = 0.1),1)
P10
P90 <- round(quantile(data_complete$t_min, probs = 0.9),1)
P90
P97.5 <- round(quantile(data_complete$t_min, probs = 0.975),1)
P97.5

# Temperature range for all microregions
bound <- round(colMeans(ranges_data[,c(1,2)]),1)
bound

load("02.models/output/Sall_meta.RData")
load("02.models/output/yall_meta.RData")
# BASES OF TEMPERATURE AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)
bvar <- do.call("onebasis",c(list(x=xvar),attr(cb_t_min,"argvar")))
bvar_data = bvar
xlag <- 0:140/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb_t_min,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall,1,function(x) exp(bvar%*%x))

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall~1,Sall,method="reml")
summary(mvall)

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall_data <- crosspred(bvar,coef=coef(mvall),
                        vcov=vcov(mvall),model.link="log",
                        by=0.1,from=bound[1],to=bound[2],cen=16)
summary(cpall_data)
# PLOT OVERALL CUMULLATIVE
png(filename = paste0("03.figs/fig05.png"), 
    height = 5, width = 7, 
    units = 'in', res = 300)

par(mar=c(4,4,4,4),mfrow=c(1,1))
plot(cpall_data,type="l",ylab="RR",xlab="TºC",main="")
abline(v=16)
dev.off()

# OVERALL EFFECTS AT PREDICTOR LEVELS (P2.5, P10, P90, P97.5)
quantis <- round(quantile(data_complete$t_min, probs = c(0.025,0.1,0.9,0.975)),1) %>% tibble(Temperature = .)
p1 <- round(with(cpall_data,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P2.5),]),2) %>% t() %>% as.tibble()
p2 <- round(with(cpall_data,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P10),]),2)%>% t() %>% as.tibble()
p3 <- round(with(cpall_data,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P90),]),2)%>% t() %>% as.tibble()
p4 <- round(with(cpall_data,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P97.5),]),2)%>% t() %>% as.tibble()
# creating table with the data
tab_meta <- bind_cols(quantis, quantis = c('P2.5','P10','P90','P97.5'), bind_rows(p1,p2,p3,p4))
write.csv(tab_meta, "02.models/tab_meta.csv")


