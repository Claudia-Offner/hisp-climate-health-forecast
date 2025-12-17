### ------------------------------------------------------------------------ ### 
###  DNLM - TEST univariate model for distribution selection                 ###
### ------------------------------------------------------------------------ ### 

set.seed(123)

# LOAD THE PACKAGE
# install.packages(c('xgboost', 'SHAPforxgboost', 'data.table', 'dplyr', 'caret'))
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

data <- df
exposure <- 't2m'
outcome <- 'total_patients'

# Function to get reference meta-analysis
get_meta_analysis <- function(ymat, Slist, arg, ranges, cen){
  
  # Meta-analysis w/ prediction basis
  mv <- mvmeta(ymat, Slist, method='ml')
  # bound <- c(min(ranges[,1]), max(ranges[,2])) #c(colMeans(ranges))
  bound <- c(colMeans(ranges))
  mean <- seq(bound[1],bound[2],length=30)
  b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
  
  return(list(ymat=ymat, Slist=Slist, arg=arg, ranges=ranges, mv=mv, bound=bound, b_mean=b_mean, cen=cen))
  
}
# Function to calculate comparative metrics for meta-DLNM
comp_metrics <- function(model, cb=NULL, cen=NULL) {
  
  res <- list()
  
  # Get sum aic
  res[['aic']] <- summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion
  
  # Get sum of ACF
  res[['acf']] <- acf(residuals(model, type = 'pearson'), plot = FALSE, lag.max = 30)$acf
  
  # Get sum of PACF
  res[['pacf']] <- pacf(residuals(model, type = 'pearson'), plot = FALSE, lag.max = 30)$acf
  
  
  if(!is.null(cb)){
    
    # Get average of estimates (lag 4)
    pred <- crosspred(cb, model, cen=cen)
    res[['esti']] <- pred$coef[4]  
    
  }
  
  return(res) 
}

# 1. Data and dependency prep ----

# LIST OF DATAFRAMES FOR ALL SITES
datalist <- lapply(locations, function(region) data[data$location==region,])
names(datalist) <- locations

# VARIABLE
ranges <- t(sapply(datalist,function(x) range(x[[exposure]],na.rm=T)))
arg <- list(argvar=3, arglag=3)
df <- prod(unlist(arg))

# 2. Stage 1: Fit the DLNM ----

# BUILT OBJECTS WHERE outputs WILL BE STORED:

Slist_comp <- vector('list',m)
names(Slist_comp) <- locations

ymat <- matrix(NA,m,df,dimnames=list(locations,paste('spl',seq(df),sep='')))
Slist <- vector('list',m)
names(Slist) <- locations

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)


system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,'')
    
    # LOAD
    data <- datalist[[i]]
    
    # City specific parameters
    data <- data %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))
    cenrel <- mean(data[[exposure]], na.rm = TRUE)
    print(cenrel)
    
    # CREATE THE SPLINES
    cb <- crossbasis(data[[exposure]], lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag), group=data$year)

    # Set formula
    form <- as.formula(paste0(outcome, ' ~  cb + 
                                patients_lag + weekday + season + 
                                ns(time, df=8*length(unique(data$year)))'))
    
    # RUN THE MODEL
    model <- glm(formula=form, family=quasipoisson(), data=data)

    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    pred <- crosspred(cb, model, cen=cenrel)
    ymat[i,] <- pred$coef
    Slist[[i]] <- pred$vcov
    
    # EXTRACT COMPARATIVE METRICS
    Slist_comp[[i]] <- comp_metrics(model)
    
  })



# 3. Stage 2: Multivariate Meta-analysis ----

# Meta coefs
meta <- get_meta_analysis(ymat=ymat,  Slist=Slist, arg=arg, ranges=ranges, cen=cen)

# Meta metrics
meta_comp <- list()
meta_comp[['sum_aic']] <- sum(sapply(Slist_comp, function(x) x$aic))
meta_comp[['sum_pacf']] <- rowSums(sapply(Slist_comp, function(x) as.numeric(x$pacf)))
