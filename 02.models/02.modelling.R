# First Stage -- GAM and DLMN
options(scipen=999)
library(tidyverse)
library(mgcv)
library(splines)
library(dlnm)
library(patchwork)
library(png)
theme_set(theme_bw())
# Preparing data ----------------------------------------------------------
# Loading data
data <- read_rds("00.data/data.rds")
# Filtering microregions with populations more than 10.000 cases in the span of time analysed
filter <- data %>% group_by(microregion_name) %>% 
  summarise(sum = sum(cases)) %>%
  filter (sum >= 10000)
data <- data %>% filter(microregion_name %in% filter$microregion_name)
# Tables and matrices for meta-analysis -----------------------------------
# OVERALL CUMULATIVE SUMMARIES
yall <- matrix(NA,length(data$microregion_name %>% unique),
               3,dimnames=list(data$microregion_name %>% 
                                 unique,paste("b",seq(3),sep="")))
# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data$microregion_name %>% unique))
names(Sall) <- data$microregion_name %>% unique
for (i in unique(data$microregion_name)){
  print(i)
  # selecting city ----------------------------------------------------------
  data_city <- data %>% filter(microregion_name == i)
  data_city$cases %>% sum
  # span of time in days
  data_city$time <- 1:3652
  # -------------------------------------------------------------------------
  ## Creating cross-basis of min temperature (lags up to 22 days)
  temp <- data_city$t_min
  cb_t_min <- crossbasis(temp, lag = 22, 
                         argvar=list(fun="ns",df=3),
                         arglag=list(fun="poly",degree=3))
  # Applying GAM with variables selected beforehand --------------------------
  mod <- gam(cases ~ ns(data_city$time,10*8) + cb_t_min + weekday,
             family = nb(), data = data_city)
  print(paste(i, " done, saving model"))
  # saving model
  save(mod, file = paste0("02.models/output/",
                          paste0('mod_',i),
                          ".RData"))
  
  # Reducing models ---------------------------------------------------------
  crall <- crossreduce(cb_t_min,mod,cen=16)
  yall[i,] <- coef(crall)
  Sall[[i]] <- vcov(crall)
 }
save(yall,file = "02.models/output/yall_meta.RData")
save(Sall,file = "02.models/output/Sall_meta.RData")

# list of models ----------------------------------------------------------

list_mod <- list.files("02.models/output/", pattern = "mod_", 
                       full.names = T)
# Applying DLNM -----------------------------------------------------------

for ( i in 1:length(list_mod)){
  #loading model
  load(list_mod[i])
  #getting name
  name <- list_mod[i] %>% 
    str_remove("02.models/output/mod_") %>% 
    str_remove(".RData") %>% 
    str_remove("/")
  print(name)
  # dlnm 
  par(mar=c(4,4,4,4))
  ## Prediction model with as reference
  pred <- crosspred(basis=cb_t_min,model=mod,
                    by=0.1,cumul=TRUE, cen = 16)
  # heat map
  png(filename = paste0("02.models/output/dlnm_plots/heat_map_weekday_",name,".png"), 
      height = 5, width = 6, 
      units = 'in', res = 200)
  plot(pred, "contour",cumul = T,lag = 22, xlab=name,ylab="RR (lag 0-22)")
  dev.off()
  # Accumulated RR
  png(filename = paste0("02.models/output/dlnm_plots/RR_acumul_", name,".png"), 
      height = 5, width = 6, 
      units = 'in', res = 200)
  plot(pred, "slices",cumul = T,lag = 22, xlab=name,ylab="RR (lag 0-22)")
  abline(v=16)
  dev.off()
  # Plot RR at specific temperature for lag range (non-cumulative)
  png(filename = paste0("02.models/output/dlnm_plots/RR_percentis_n_cumulative", name,".png"), 
      height = 5, width = 7, 
      units = 'in', res = 200)
  par(mar=c(4,4,4,4),mfrow=c(2,2))
  
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.025),1)),lty=1,cumul=F,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.025),1),name,"(P2.5)"))
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.1),1)),lty=1,cumul=F,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.1),1),name,"(P10)"))
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.9),1)),lty=1,cumul=F,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.9),1),name,"(P90)"))
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.975),1)),lty=1,cumul=F,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.975),1),name,"(P97.5)"))
  dev.off()
  ## Plot RR at specific temperature for lag range (cumulative)
  png(filename = paste0("02.models/output/dlnm_plots/RR_percentis_cumulative", name,".png"), 
      height = 5, width = 7, 
      units = 'in', res = 200)
  par(mar=c(4,4,4,4),mfrow=c(2,2))
  
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.025),1)),lty=1,cumul = T,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.025),1),name,"(P2.5)"))
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.1),1)),lty=1,cumul = T,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.1),1),name,"(P10)"))
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.9),1)),lty=1,cumul = T,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.9),1),name,"(P90)"))
  plot(pred,"slices",
       var=as.character(round(quantile(temp,probs=0.975),1)),lty=1,cumul = T,ylab="RR", 
       main=paste("Lag-response curve at ",
                  round(quantile(temp,probs=0.975),1),name,"(P97.5)"))
  
  dev.off()
}

# data regarding RR of each temperature span analyzed 


table_RR <- tibble()

mod_files <- list.files("02.models/output/", pattern = "mod", full.names = T)


for (i in 1:38){
  
  # Loading models
  load(mod_files[i])
  
  ## Prediction model with 16 degrees Celsius
  
  pred <- crosspred(basis=cb_t_min,model=mod,
                    by=0.1,cumul=TRUE, cen = 16)
  
  ## RR accumulated total for the entire period for percentiles 2.5/10/90/97.5
  percentile_2.5 <- round(cbind(pred$allRRfit,
                                pred$allRRlow,
                                pred$allRRhigh)[round(quantile(temp,probs=0.025),1) %>% as.numeric(),],2) %>% 
    tibble(percentile_2.5 = .)
  percentile_10 <- round(cbind(pred$allRRfit,
                               pred$allRRlow,
                               pred$allRRhigh)[round(quantile(temp,probs=0.1),1) %>% as.numeric(),],2) %>%
    tibble(percentile_10 = .)
  percentile_90 <- round(cbind(pred$allRRfit,
                               pred$allRRlow,
                               pred$allRRhigh)[round(quantile(temp,probs=0.9),1) %>% as.numeric(),],2) %>%
    tibble(percentile_90 = .)
  percentile_97.5 <- round(cbind(pred$allRRfit,
                                 pred$allRRlow,
                                 pred$allRRhigh)[round(quantile(temp,probs=0.975),1) %>% as.numeric(),],2)%>%
    tibble(percentile_97.5 = .)
  
  
  lil_RR_n_CI <- bind_cols(variables = c("RR", "IC min", "IC max"),
                           percentile_2.5,
                           percentile_10,
                           percentile_90,
                           percentile_97.5,
                           Region = list_mod[i] %>% 
                             str_remove("06.models/output_cities/mod_") %>% 
                             str_remove(".RData"))
  
  RR_n_CI <- lil_RR_n_CI %>% pivot_longer(cols = c(percentile_2.5, 
                                                   percentile_10,
                                                   percentile_90,
                                                   percentile_97.5)) %>% 
    pivot_wider(names_from = variables, values_from = value)
  table_RR <- bind_rows(table_RR, RR_n_CI)
}

xlsx::write.xlsx(table_RR, "02.models/table_RR.xlsx")
