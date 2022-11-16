
# Selecting model ---------------------------------------------------------
options(scipen=999)
library(tidyverse)
library(splines)
library(dlnm)
library(mgcv)
library(splines)
theme_set(theme_bw())
# Preparing data ----------------------------------------------------------
# Loading data
data <- read_rds("00.data/data.rds")
# Filtering microregions with populations more than 10.000 cases in the span of time analysed
filter <- data %>% group_by(microregion_name) %>%
  summarise(sum = sum(cases)) %>%
  filter (sum >= 10000)
# Selecting DF of time spline
data <- data %>% filter(microregion_name %in% filter$microregion_name)
for (a in 4:8){
  sum_AIC <- c(0)
  
  for (i in unique(data$microregion_name)){
    print(i)
    # selecting city ----------------------------------------------------------
    data_city <- data %>% filter(microregion_name == i)
    data_city$cases %>% sum %>% print
    # span of time in days
    data_city$time <- 1:3652
    # -------------------------------------------------------------------------
    ## Creating cross-basis of min temperature (lags up to 21 days)
    temp <- data_city$t_min
    cb_t_min <- crossbasis(temp, lag = 21, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="poly",degree=3))
    # Applying GAM with variables selected beforehand --------------------------
    mod <- gam(cases ~ ns(data_city$time,10*a) + cb_t_min + weekday,
               family = nb(), data = data_city)
    print(paste(i, " done, saving model"))
    # AIC
    sum_AIC <- mod$aic + sum_AIC
  }
  assign( paste0("sum_AIC",a),sum_AIC)
}
sum_AIC4
sum_AIC5
sum_AIC6
sum_AIC7
sum_AIC8
# the model with the lowest AIC was with 8 df per year in the time spline

# Next, we seek to see if the addition of population as an offset (#4)
# the use a lag span based on Incubation periods leas to a better result of the AIC sum.

# Previous analysis used a lag based on the an exploratory analysis (21 - days)
# In order to apply biological features we'll apply the span based
# on the sum of the described previously incubation periods
# 8-12 days for EIP and 3-10 days for IIP, then: 11, 15  and 22 days 
# (Chan & Johansson,2012)

sum_AIC_1 <- c(0)
sum_AIC_2 <- c(0)
sum_AIC_3 <- c(0)

for (i in unique(data$microregion_name)){
  # selecting city ----------------------------------------------------------
  print(i)
  data_city <- data %>% filter(microregion_name == i)
  data_city$cases %>% sum %>% print
  # span of time in days
  data_city$time <- 1:3652
  temp <- data_city$t_min
  # ------------------------------------------------------------------------
  # Model 1 -----------------------------------------------------------------
  # Longer lag : 22 days   # 1
  cb_t_min_long <- crossbasis(temp, lag = 22, 
                              argvar=list(fun="ns",df=3),
                              arglag=list(fun="poly",degree=3))
  mod <- gam(cases ~ ns(data_city$time,10*8) + cb_t_min_long + weekday,
             family = nb(), data = data_city)
  
  # AIC Sum ---------------------------------------------------------
  sum_AIC_1 <- mod$aic + sum_AIC_1
  print(paste(i, " done model 1"))
  
  # Model 2 -----------------------------------------------------------------
  # Longer lag : 15 days   # 2
  cb_t_min_mid <- crossbasis(temp, lag = 15, 
                             argvar=list(fun="ns",df=3),
                             arglag=list(fun="poly",degree=3))
  mod <- gam(cases ~ ns(data_city$time,10*8) + cb_t_min_mid + weekday,
             family = nb(), data = data_city)
  # AIC Sum ---------------------------------------------------------
  sum_AIC_2 <- mod$aic + sum_AIC_2
  print(paste(i, " done model 2"))
  # Model 7 -----------------------------------------------------------------
  # Longer lag : 11 days   # 7
  cb_t_min_short <- crossbasis(temp, lag = 11, 
                               argvar=list(fun="ns",df=3),
                               arglag=list(fun="poly",degree=3))
  mod <- gam(cases ~ ns(data_city$time,10*8) + cb_t_min_short + weekday,
             family = nb(), data = data_city)
  # AIC Sum ---------------------------------------------------------
  sum_AIC_3 <- mod$aic + sum_AIC_3
  print(paste(i, " done model 3"))
}

sum_AIC_1
sum_AIC_2
sum_AIC_3

# -------------------------------------------------------------------------
# The model with the lowest AIC was the one with the longer span of lag, in 
# similar way as the exploratory analysis

# Introducing population as an offset. We expect that the AIC would not change 
# much, since the time spline already controls changes in a long span of time

sum_AIC_4 <- c(0)

for (i in unique(data$microregion_name)){
  # selecting city ----------------------------------------------------------
  print(i)
  data_city <- data %>% filter(microregion_name == i)
  data_city$cases %>% sum %>% print
  # span of time in days
  data_city$time <- 1:3652
  temp <- data_city$t_min
  # ------------------------------------------------------------------------
  # Model 4 -----------------------------------------------------------------
  # population as an offset   # 4
  cb_t_min_long <- crossbasis(temp, lag = 22, 
                              argvar=list(fun="ns",df=3),
                              arglag=list(fun="poly",degree=3))
  mod <- gam(cases ~ ns(data_city$time,10*8) + cb_t_min_long + weekday + offset(log(pop)),
             family = nb(), data = data_city)
  
  # AIC Sum ---------------------------------------------------------
  sum_AIC_4 <- mod$aic + sum_AIC_4
  print(paste(i, " done model 4"))
}

sum_AIC_4

# The value of the AIC changed slightly with the population as an offset7
AIC_table <- tibble(name = c("4 DF","5 DF","6 DF","7 DF","8 DF", 
                "1-22 lag", "1-15 lag", "1-11 lag", "Population as offset"),
       AIC_sum = c(sum_AIC4,sum_AIC5,sum_AIC6,sum_AIC7,sum_AIC8,
                   sum_AIC_1,sum_AIC_2,sum_AIC_3,sum_AIC_4))

xlsx::write.xlsx(AIC_table, "02.models//AIC_table.xlsx" )
