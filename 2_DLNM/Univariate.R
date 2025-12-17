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

# FUNCTIONS ####

# Function to run model with static formula
run_uni_model <- function(data, locations, outcome, exposure, fam_type='binom'){
  
  # Function to get reference meta-analysis
  get_meta_analysis <- function(ymat, Slist, arg, ranges, cen){
    
    # Meta-analysis w/ prediction basis
    mv <- mvmeta(ymat, Slist, method='ml')
    # bound <- c(min(ranges[,1]), max(ranges[,2])) #c(colMeans(ranges))
    bound <- c(colMeans(ranges))
    mean <- seq(bound[1],bound[2],length=30)
    # b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
    b_mean <- onebasis(data[[exposure]], fun="ns", df=arg$argvar, intercept=FALSE)
    
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
  cen <- mean(data[[exposure]], na.rm = TRUE)
  arg <- list(argvar=5, arglag=1)
  df <- prod(unlist(arg))
  
  # # DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
  # knotperc <- c(10,25,75,90)
  # knots <- rowMeans(sapply(datalist,function(x) 
  #   quantile(x$t2m,knotperc/100,na.rm=T)))
  
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
      
      data <- data %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))
      
      # CREATE THE SPLINES
      # cb <- crossbasis(data[[exposure]], lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag), group=data$year)
      cb <- onebasis(data[[exposure]], fun = "ns", df = arg$argvar, intercept=FALSE)
      
      # Set formula
      form <- as.formula(paste0(outcome, ' ~  cb + 
                                patients_lag + weekday + season + 
                                ns(time, df=8*length(unique(data$year)))'))
      
      # RUN THE MODEL
      if(fam_type=='quasi'){
        model <- glm(formula=form, family=quasipoisson(), data=data)
      } else if (fam_type=='binom'){
        model <- glm.nb(formula=form,data=data)
      } else {
        stop('Unknown model type: ', fam_type)
      }
      
      # EXTRACT AND SAVE THE RELATED COEF AND VCOV
      pred <- crosspred(cb, model, cen=cen)
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
  
  # RESET WARNING
  options(warn=0)
  
  return(list(meta_comp=meta_comp, meta=meta))
  
}

# Function to extract and format relevant results and figures
format_uni_results <- function(output, outcome, exposure, color='#52C5B0FF'){
  
  # Function to get usable tables and figures
  get_results <- function(output, city, lags, probs=0.5, outcome){
    
    # Function to extract model comparison criteria
    ftab <- function(m, mref=NULL) {
      
      # HETEROGENEITY AND IC STATS
      q <- qtest(m)
      het <- c(q$Q[1],q$df[1],q$pvalue[1],(q$Q[1]-q$df[1])/q$Q[1]*100,AIC(m),BIC(m))
      # LR TEST (ONLY FOR META-REGRESSION)
      if(!is.null(mref)) {
        lrstat <- -2*(logLik(mref)-logLik(m))[1]
        df <- attr(logLik(m),'df')-attr(logLik(mref),'df')
        pvalue <- 1-pchisq(lrstat,df)
        lr <- c(lrstat,df,pvalue)
      }
      # WALD TEST (ONLY FOR META-REGRESSION)
      if(!is.null(mref)) {
        coef <- coef(m)[-grep('Int',names(coef(m)))]
        vcov <- vcov(m)[-grep('Int',names(coef(m))),-grep('Int',names(coef(m)))]
        waldstat <- coef%*%solve(vcov)%*%coef
        df <- length(coef)
        pvalue <- 1-pchisq(waldstat,df)
        wald <- c(waldstat,df,pvalue)
      }
      # outputS
      if(!is.null(mref)) {
        return(c(het,lr,wald))
      } else return(c(het,rep(NA,6))) 
    }
    
    meta <- output$meta
    
    # Model selection criteria
    tab <- matrix(NA,1,12)
    colnames(tab) <- c('Q','df','p','I-square','AIC','BIC','stat','df','p','stat','df','p')
    rownames(tab) <- c(exposure)
    tab[1,] <- ftab(m=meta$mv)

    # Exposure response (with city specs) - Overall
    p_er <- city_exp_res(city=city, lag=lags, meta=meta, col=color, xlab='Exposure', title=NULL)

    # Lag response (with city specs) - e.g. 95th percentile
    q <- quantile(df[[exposure]], probs=probs, na.rm = TRUE)
    p_lr <- city_lag_res(city=city, quant_v=q, meta=meta, col=color, xlab='Exposure', title=NULL)

    return(list(tab=tab, p_er=p_er, p_lr=p_lr))
    
  }
  
  overall <- get_results(output, city='ALL', lags=NULL, outcome=outcome)
  all <- overall[!names(overall) %in% 'tab']
  p_blup_all <- lapply(all, function(x) x[['p_blup']])
  p_blup_er <- p_blup_all[grepl('er', names(p_blup_all))]
  # p_blup_lr <- p_blup_all[grepl('lr', names(p_blup_all))]
  p_er_overall <- wrap_plots(p_blup_er, ncol=1)+plot_annotation(title=exposure, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')))
  
  # quant_test <- c(0.1, 0.25, 0.75, 0.9)
  # quant <- list()
  # for (q in quant_test){
  #   quant[[as.character(q)]] <- get_results(output, city='ALL', lags=NULL, probs=q, outcome=outcome)
  # }
  # p_lr <- (wrap_plots(lapply(quant, function(x) x$p_lr$p_blup), ncol = 4) +
  #            plot_annotation(title = exposure, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))) & 
  #   theme(legend.position = "none")
  # 
  # lag_test <- c(0, 5, 10, 14)
  # lags <- list()
  # for (l in lag_test){
  #   lags[[as.character(l)]] <- get_results(output, city='ALL', lags=l, outcome=outcome)
  # }
  # p_er <- (wrap_plots(lapply(lags, function(x) x$p_er$p_blup), ncol = 4) +
  #            plot_annotation(title = exposure, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))) &
  #   theme(legend.position = "none")
  
  city <- list()
  for (i in locations){
    city[[i]] <- get_results(output, city=i, lags=NULL, outcome=outcome)
  }
  p_er_city <- wrap_plots(lapply(city, function(x) x$p_er$p_city), ncol=4)+plot_annotation(title=exposure, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')))
  p_lr_city  <- wrap_plots(lapply(city, function(x) x$p_lr$p_city), ncol=4)+plot_annotation(title=exposure, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')))
 
  return(list(overall=overall, p_er_overall=p_er_overall, 
              city=city, 
              p_er_city=p_er_city, p_lr_city=p_lr_city))
}

# MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# 1. Check distributions ####

# Total patients
primary <- 'total_patients'
exposure <- list('t2m'=list(color='#D74B3C'),
                 'mn2t'=list(color='#A1C8F2'),
                 'mx2t'=list(color='#F1A6A6'),
                 'tp'=list(color='#1D98EB'), 
                 'fg10'=list(color='#A6B727'), 
                 'pm25'=list(color='#8968CD'), 
                 'rh'=list(color='#FEC306'), 
                 'dt'=list(color='#E569BF'), 
                 'wb'=list(color='#EC800A') )
total_overall <- list()
total_lr <- list()
total_er <- list()

for (e in names(exposure)){
  output <- run_uni_model(data=df, locations=locations, outcome=primary, exposure=e)
  output$meta_comp$sum_aic
  main_res <- format_uni_results(output=output, outcome=primary, exposure=e, color=exposure[[e]]$color)
  total_overall[[e]] <- main_res$p_er_overall + plot_annotation(title = replace_exact_matches(e, nm))
  # total_lr[[e]] <- main_res$p_lr + plot_annotation(title = replace_exact_matches(e, nm))
  # total_er[[e]] <- main_res$p_er + plot_annotation(title = replace_exact_matches(e, nm))
}

(p_er_overall <- wrap_plots(total_overall, ncol=3)+plot_annotation(title=primary, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))

# total_overall
# total_lr
# total_er
#######################################################


