## MULTIVARIATE META-REGRESSION
library(tidyverse)
library(mgcv)
library(mvmeta) 
library(dlnm)
library(patchwork)
options(scipen=999)
theme_set(theme_bw())

# Loading data for meta regression
data_mt <- read_rds("00.data/data.rds") 
filter <- data_mt %>% group_by(microregion_name) %>% 
  summarise(sum = sum(cases)) %>%
  filter (sum >= 10000)
data_mt <- data_mt %>% filter(microregion_name %in% filter$microregion_name)
cen  <-  round(data_mt$t_min %>% mean,0)
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

variables <- c('População', 'PIB', 'PIB per capta', 'Altitude', 'Acesso à Tratamento de Água',
               'Acesso à Coleta de Lixo', 'Porcetagem de População Urbana', 'Latitute', 'Longitude')
wald_values <- bind_rows(wald_pop,wald_GDP,wald_pcGDP,wald_elev,wald_water,
                         wald_waste,wald_urban,wald_lat,wald_long) 
colnames(wald_values) <- c('Wald_test')
wald_tab <- bind_cols(Variables = variables, wald_values)
xlsx::write.xlsx(wald_tab, "02.models/wald_tab.xlsx" )

## PREDICTION FROM META-REGRESSION
# elevation
val.elev <- round(quantile(meanelev,na.rm=T,c(10,90)/100),1)
val.elev
pred.elev <- predict(mvelev,data.frame(val.elev),vcov=T)
cpall_data_elev10 <- crosspred(bvar,coef=pred.elev[[1]]$fit,vcov=pred.elev[[1]]$vcov,model.link="log",by=0.1,cen=cen)
cpall_data_elev90 <- crosspred(bvar,coef=pred.elev[[2]]$fit,vcov=pred.elev[[2]]$vcov,model.link="log",by=0.1,cen=cen)
# plotting image

plot(cpall_data,type="l",ylab="RR",xlab="TºC",main="")
lines(cpall_data_elev10,type="l",ylab="RR",xlab="TºC",main="P10",col="blue",ci="lines")
lines(cpall_data_elev90,type="l",ylab="RR",xlab="TºC",main="P90",col="red",ci="lines")
text(10,1.54,'I² = 58.1%')

# Estimates for p10 of elevation
quantis <- round(quantile(data_complete$t_min, probs = c(0.025,0.1,0.9,0.975)),1) %>% tibble(Temperature = .)
p1 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P2.5),]),2) %>% t() %>% as.tibble()
p2 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P10),]),2) %>% t() %>% as.tibble()
p3 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P90),]),2) %>% t() %>% as.tibble()
p4 <- round(with(cpall_data_elev10,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P97.5),]),2) %>% t() %>% as.tibble()
tab_meta_reg_p10 <- 
  bind_cols(quantis, quantis = c('P2.5','P10','P90','P97.5'), 
            bind_rows(p1,p2,p3,p4)) %>% 
  mutate(RR = paste0(allRRfit," (",allRRlow,"—",allRRhigh,")")) %>% select(c(1,2,6)) 

# Estimates for p90 of elevation
quantis <- round(quantile(data_complete$t_min, probs = c(0.025,0.1,0.9,0.975)),1) %>% tibble(Temperature = .)
p1 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P2.5),]),2) %>% t() %>% as.tibble()
p2 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P10),]),2) %>% t() %>% as.tibble()
p3 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P90),]),2) %>% t() %>% as.tibble()
p4 <- round(with(cpall_data_elev90,cbind(allRRfit,allRRlow,allRRhigh)[as.character(P97.5),]),2) %>% t() %>% as.tibble()
tab_meta_reg_p90 <- 
  bind_cols(quantis, quantis = c('P2.5','P10','P90','P97.5'),
            bind_rows(p1,p2,p3,p4))%>% 
  mutate(RR = paste0(allRRfit," (",allRRlow,"—",allRRhigh,")")) %>% select(c(1,2,6)) 


write.csv(tab_meta_reg_p10, "02.models/tab_meta_reg_p10.csv")
write.csv(tab_meta_reg_p90, "02.models/tab_meta_reg_p90.csv")

# Plotting combined image of baseline and elevation models
p1 <- tibble(fit = cpall_data$allRRfit, low = cpall_data$allRRlow,
             high = cpall_data$allRRhigh, temp = cpall_data$predvar) %>% 
  ggplot(aes(x =temp , y = fit))+
  geom_line()+
  geom_ribbon(aes(temp, ymin = low, ymax = high),
              fill = "grey", alpha = 0.5)+
  geom_vline(xintercept = cen)+
  geom_hline(yintercept = 1)+
  xlab("TºC")+
  ylab("RR")+
  geom_text(aes(x= 10, y = 1.4, label = 'I² = 60,0%'))+
  scale_y_continuous(limits=c(0.4,1.4))

p2 <- tibble(fit = cpall_data$allRRfit, 
             low = cpall_data$allRRlow,
             high = cpall_data$allRRhigh, 
             temp = cpall_data$predvar,
             fit10 = cpall_data_elev10$allRRfit, 
             low10 = cpall_data_elev10$allRRlow,
             high10 = cpall_data_elev10$allRRhigh,
             fit90 = cpall_data_elev90$allRRfit,
             low90 = cpall_data_elev90$allRRlow,
             high90 = cpall_data_elev90$allRRhigh) %>% 
  ggplot()+
  geom_line(aes(x = temp , y = fit))+
  geom_ribbon(aes(temp, ymin = low, ymax = high),
              fill = "grey", alpha = 0.5)+
  geom_line(aes(x =temp , y = fit10, color = "p10"))+
  geom_line(aes(x =temp , y = low10, color = "p10"), linetype = "dashed")+
  geom_line(aes(x =temp , y = high10,color = "p10"), linetype = "dashed")+
  geom_line(aes(x =temp , y = fit90, color = "p90"))+
  geom_line(aes(x =temp , y = low90, color = "p90"), linetype = "dashed")+
  geom_line(aes(x =temp , y = high90,color = "p90"), linetype = "dashed")+
  geom_vline(xintercept = cen)+
  geom_hline(yintercept = 1)+ labs(color='Elevation')+
  xlab("TºC")+
  ylab("RR")+
  geom_text(aes(x= 10, y = 1.4, label = 'I² = 58,1%'))+
  scale_y_continuous(limits=c(0.4,1.4))
# layout for plot
layout <- "
AAABBB
AAABBB
CCDDEE
"
p1+
  p2+
  tableGrob(mutate_all(tab_meta, ~ str_replace_all(.,"\\.",",")), rows = c("","","",""),cols = c("Temperatura","Percentis","RR"))+
  tableGrob(mutate_all(wald_tab,~ str_replace_all(.,"\\.",",")), rows = c("","","","","","","","",""), cols= c("Variáveis", "Teste de Wald"))+
  tableGrob(bind_cols(tab_meta_reg_p10,
                      tab_meta_reg_p90[,3]) %>% mutate_all(.,~ str_replace_all(.,"\\.",",")),
            cols = c("Temperatura","Percentis","RR (p10)", "RR (p90)"),
            rows = c("","","",""))+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')

ggsave("03.figs/fig02.png", height = 10, width = 16)

