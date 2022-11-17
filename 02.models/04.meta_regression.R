## MULTIVARIATE META-REGRESSION
library(tidyverse)
library(mgcv)
library(mvmeta) 
library(dlnm)
options(scipen=999)
# Loading data for meta regression
data_mt <- read_rds("00.data/data.rds") 
filter <- data_complete %>% group_by(microregion_name) %>% 
  summarise(sum = sum(cases)) %>%
  filter (sum >= 10000)
data_mt <- data_mt %>% filter(microregion_name %in% filter$microregion_name)
data_mt <- lapply(regions,function(x) data_mt[data_mt$microregion_name==x,])
names(data_mt) <- regions
# Wald test
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))}

# pop
meanpop <- sapply(data_mt, function(x) mean(x$pop,na.rm=T))
round(quantile(meanpop,na.rm=T,c(10,90)/100))
mvpop <- update(mvall,.~meanpop)
wald_pop<- round(fwald(mvpop,"pop"),3)%>% tibble()
summary(mvpop)

# GDP
meanGDP <- sapply(data_mt, function(x) mean(x$PIB,na.rm=T))
round(quantile(meanGDP,na.rm=T,c(10,90)/100),1)
mvGDP <- update(mvall,.~meanGDP)
wald_GDP <- round(fwald(mvGDP,"meanGDP"),3)%>% tibble()
summary(mvGDP)

# pcGDP
meanpcGDP <- sapply(data_mt, function(x) mean(x$pibpc,na.rm=T))
round(quantile(meanpcGDP,na.rm=T,c(10,90)/100))
mvpcGDP <- update(mvall,.~meanpcGDP)
wald_pcGDP <- round(fwald(mvpcGDP,"meanpcGDP"),3) %>% tibble()
summary(mvpcGDP)

# elevation
meanelev <- sapply(data_mt, function(x) mean(x$elevation,na.rm=T))
round(quantile(meanelev,na.rm=T,c(10,90)/100),1)
mvelev <- update(mvall,.~meanelev)
wald_elev <- round(fwald(mvelev,"meanelev"),3)%>% tibble()
summary(mvelev)

# Water access
meanwater <- sapply(data_mt, function(x) mean(x$pagua2010,na.rm=T))
round(quantile(meanwater,na.rm=T,c(10,90)/100),1)
mvwater <- update(mvall,.~meanwater)
wald_water <- round(fwald(mvwater,"meanwater"),3)%>% tibble()
summary(mvwater)

# Waste management access 
meanwaste <- sapply(data_mt, function(x) mean(x$plixo2010,na.rm=T))
round(quantile(meanwaste,na.rm=T,c(10,90)/100),1)
mvwaste <- update(mvall,.~meanwaste)
wald_waste <- round(fwald(mvwaste,"meanwaste"),3)%>% tibble()
summary(mvwaste)

# Urban population
meanurban <- sapply(data_mt, function(x) mean(x$purb2010,na.rm=T))
round(quantile(meanurban,na.rm=T,c(10,90)/100),1)
mvurban <- update(mvall,.~meanurban)
wald_urban <- round(fwald(mvurban,"meanurban"),3)%>% tibble()
summary(mvurban)

# latitude
meanlat <- sapply(data_mt, function(x) mean(x$lat,na.rm=T))
round(quantile(meanlat,na.rm=T,c(10,90)/100),2)
mvlat <- update(mvall,.~meanlat)
wald_lat <- round(fwald(mvlat,"meanlat"),3)%>% tibble()
summary(mvlat)

# longitude
meanlon <- sapply(data_mt, function(x) mean(x$lon,na.rm=T))
round(quantile(meanlon,na.rm=T,c(10,90)/100),2)
mvlong <- update(mvall,.~meanlon)
wald_long <- round(fwald(mvlong,"meanlon"),3)%>% tibble()
summary(mvlong)

variables <- c('Population', 'GDP', 'Per capita GDP', 'Elevation', 'Treated Water Access',
               'Waste Management', 'Urban Population', 'Latitute', 'Longitude')
wald_values <- bind_rows(wald_pop,wald_GDP,wald_pcGDP,wald_elev,wald_water,
                         wald_waste,wald_urban,wald_lat,wald_long) 
colnames(wald_values) <- c('Wald_test')
wald_tab <- bind_cols(Variables = variables, wald_values)

write.csv(wald_tab, "02.models/wald_tab.csv" )

## PREDICTION FROM META-REGRESSION
# elevation
val.elev <- round(quantile(meanelev,na.rm=T,c(10,90)/100),1)
val.elev
pred.elev <- predict(mvelev,data.frame(val.elev),vcov=T)
cpall_data_elev10 <- crosspred(bvar,coef=pred.elev[[1]]$fit,vcov=pred.elev[[1]]$vcov,model.link="log",by=0.1,cen=16)
cpall_data_elev90 <- crosspred(bvar,coef=pred.elev[[2]]$fit,vcov=pred.elev[[2]]$vcov,model.link="log",by=0.1,cen=16)
# saving image
png(filename = paste0("03.figs/fig06.png"), 
    height = 5, width = 7, 
    units = 'in', res = 300)

plot(cpall_data,type="l",ylab="RR",xlab="TºC",main="")
lines(cpall_data_elev10,type="l",ylab="RR",xlab="TºC",main="P10",col="blue",ci="lines")
lines(cpall_data_elev90,type="l",ylab="RR",xlab="TºC",main="P90",col="red",ci="lines")
dev.off()

# Estimates for p10 of elevation
quantis <- round(quantile(data_complete$t_min, probs = c(0.025,0.1,0.9,0.975)),1) %>% tibble(Temperature = .)
p1 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P2.5),]),2) %>% t() %>% as.tibble()
p2 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P10),]),2) %>% t() %>% as.tibble()
p3 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P90),]),2) %>% t() %>% as.tibble()
p4 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P97.5),]),2) %>% t() %>% as.tibble()
tab_meta_reg_p10 <- bind_cols(quantis, quantis = c('P2.5','P10','P90','P97.5'), bind_rows(p1,p2,p3,p4))

# Estimates for p90 of elevation
quantis <- round(quantile(data_complete$t_min, probs = c(0.025,0.1,0.9,0.975)),1) %>% tibble(Temperature = .)
p1 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P2.5),]),2) %>% t() %>% as.tibble()
p2 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P10),]),2) %>% t() %>% as.tibble()
p3 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P90),]),2) %>% t() %>% as.tibble()
p4 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P97.5),]),2) %>% t() %>% as.tibble()
tab_meta_reg_p90 <- bind_cols(quantis, quantis = c('P2.5','P10','P90','P97.5'), bind_rows(p1,p2,p3,p4))


write.csv(tab_meta_reg_p10, "02.models/tab_meta_reg_p10.csv")
write.csv(tab_meta_reg_p90, "02.models/tab_meta_reg_p90.csv")
