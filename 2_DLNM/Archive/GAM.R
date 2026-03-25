### ------------------------------------------------------------------------ ### 
###  DNLM - uni_test GAM model for univariate/multivariate modelling             ###
### ------------------------------------------------------------------------ ### 

set.seed(123)

# LOAD THE PACKAGE
# install.packages(c('xgboost', 'SHAPforxgboost', 'data.table', 'dplyr', 'caret'))
library(mgcv)
library(patchwork)
library(ggplotify)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

df <- df %>% filter(!location %in% c("Bilaspur", "Reckong Peo"))
df$season  <- relevel(df$season, ref = "Pre-Monsoon (Mar-May)")
df$weekday <- relevel(df$weekday, ref = "Sun")


# 2. Fit GAMM
model <- gamm(
  total_patients ~ 
    s(location, bs = "re") +          # random intercept for location
    elev + popu + season + weekday +  # fixed effects
    s(time, k = 10 * length(unique(df$year))) +  # long-term trend
    s(t2m) + s(tp) + s(fg10) + s(rh) + s(pm25),  # smooth expos.
  family = poisson(link = "log"),
  # correlation = corAR1(form = ~ time | location)   # AR1 within each location
  data = df
)
summary(model)
plot(model$gam)


# DLNM Univariate ####

uni_test <- function(df, variable, x){
  
  variable <<- variable
  x <<- x
  
  # 1. Create crossbasis
  cb <<- crossbasis(df[[variable]], lag=21,
                   argvar=list(fun="ns", df=4),
                   arglag=list(fun="ns", df=3))

  # 2. Fit GAMM
  model <<- gamm(total_patients ~ cb + s(time, k = 10*length(unique(df$year))) +
                   season + weekday + elev + popu,
                random = list(location = ~1),  # random intercept
                family = poisson(link="log"),
                # correlation = corAR1(form = ~ time | location),
                data = df)
  
  pred <<- crosspred(cb, model$gam,
                    cen = mean(df[[variable]], na.rm=TRUE),
                    bylag = 1)  
  
  # 3. Plot univariate outputs
  level <- as.ggplot(~plot(pred, "overall", xlab=variable, ylab="Relative Risk",
                           main=paste0("Cumulative ", variable, " Effect")))
  
  lags <- as.ggplot(~plot(pred, "slices", var=x, main=paste0("RR over lags at ", variable, " = ", x)))
  
  all <- as.ggplot(~plot(pred, "contour", xlab=variable, ylab="Lag (days)",
                         key.title=title("RR"),
                         main=paste0(variable, "-Lag-Response")))
  
  res <- level + lags + all
  
  return(list(res=res, level=level, lags=lags, all=all))
}

variables <- list('t2m'=0,
                  'tp'=3,
                  'pm25'=100,
                  'fg10'=50,
                  'rh'=80)
uni_res <- list()
for (v in names(variables)){
  out <- uni_test(df, variable=v, x= variables[[v]])
  uni_res[[v]] <- out
}

uni_res$t2m$all
uni_res$tp$all
uni_res$pm25$all
uni_res$fg10$all
uni_res$rh$all


# Multivariate ####


uni_test <- function(df, variable, x){
  
  variable <<- variable
  x <<- x
  
  # 1. Create crossbasis
  cb.t2m <<- crossbasis(df[['t2m']], lag=21, argvar=list(fun="ns", df=4), arglag=list(fun="ns", df=3))
  
  # 2. Fit GAMM
  model <<- gamm(total_patients ~ cb + s(time, k = 10*length(unique(df$year))) +
                   season + weekday + elev + popu,
                 random = list(location = ~1),  # random intercept
                 family = poisson(link="log"),
                 # correlation = corAR1(form = ~ time | location),
                 data = df)
  
  pred <<- crosspred(cb, model$gam,
                     cen = mean(df[[variable]], na.rm=TRUE),
                     bylag = 1)  
  
  # 3. Plot univariate outputs
  level <- as.ggplot(~plot(pred, "overall", xlab=variable, ylab="Relative Risk",
                           main=paste0("Cumulative ", variable, " Effect")))
  
  lags <- as.ggplot(~plot(pred, "slices", var=x, main=paste0("RR over lags at ", variable, " = ", x)))
  
  all <- as.ggplot(~plot(pred, "contour", xlab=variable, ylab="Lag (days)",
                         key.title=title("RR"),
                         main=paste0(variable, "-Lag-Response")))
  
  res <- level + lags + all
  
  return(list(res=res, level=level, lags=lags, all=all))
}

