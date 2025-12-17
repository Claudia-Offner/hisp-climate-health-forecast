### ------------------------------------------------------------------------ ### 
### HISP DLNM - Dependencies                                                 ###
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/

# Turn off scientific notation
options(scipen=999)

#### PATHS ####
user <- 'claer14'
### Project
path_data <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
path_ch2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 2 - Causal Inference/')
### Git
path_git <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path_project <- paste0(path_git, '2_DLNM/')
path_report <- paste0(path_project, '- Reporting/')

#### PACKAGES ####

# List packages
packages <- c('dplyr', 'tidyr', 'lubridate', 'zoo', 'dlnm', 'splines', 
              'ggplot2', 'patchwork')

# library(MASS); library(ggfortify); library(gridExtra)

# Install if needed
installed <- packages %in% rownames(installed.packages())
if(any(!installed)) {
  install.packages(packages[!installed])
}

# Load packages
lapply(packages, library, character.only=TRUE)

#### FUNCTIONS ####

# Function to replace exact matches from a list
replace_exact_matches <- function(df, dict) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      col <- ifelse(col %in% names(dict), dict[col], col)
    }
    col
  })
  df
}

# Function to calculate comparative metrics for meta-DLNM
comp_metrics <- function(model, cb=NULL, cen=NULL) {
  
  res <- list()
  
  # Get sum aic
  res[['aic']] <- summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion
  
  # Get sum of ACF
  res[['acf']] <- acf(residuals(model, type='pearson'), plot=FALSE, lag.max=30)$acf
  
  # Get sum of PACF
  res[['pacf']] <- pacf(residuals(model, type='pearson'), plot=FALSE, lag.max=30)$acf
  
  
  if(!is.null(cb)){
    
    # Get average of estimates (lag 4)
    pred <- crosspred(cb, model, cen=cen)
    res[['esti']] <- pred$coef[4]  
    
  }
  
  return(res) 
}

# Function to run model with static formula
run_model <- function(data, locations, outcome, fam_type='binom', reference='50%'){
  
  # Function to get reference meta-analysis
  get_meta_analysis <- function(exposure, ymat, Slist, arg, ranges, cen){
    
    # Meta-analysis w/ prediction basis
    mv <- mvmeta(ymat, Slist, method='ml')
    bound <- c(colMeans(ranges))
    mean <- seq(bound[1],bound[2],length=30)
    b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag)) #, bound=bound

    return(list(ymat=ymat, Slist=Slist, arg=arg, ranges=ranges, mv=mv, 
                bound=bound, b_mean=b_mean, cen=cen))
    
  }
  
  # 1. Data and dependency prep ----
  
  # LIST OF DATAFRAMES FOR ALL SITES
  datalist <- lapply(locations, function(region) data[data$location==region,])
  names(datalist) <- locations
  
  # # ADDITIONAL CITY-LEVEL META-PREDICTORS
  # lat <- unique(data$lat)
  # perclat <- round(quantile(lat,c(1,3)*0.25),1)
  
  # wb
  ranges_wb <- t(sapply(datalist,function(x) range(x$wb,na.rm=T)))
  cen_wb <- sapply(datalist,function(x) {
    quantile(x$wb,c(0:100)/100,na.rm=T)
  })
  arg_wb <- list(argvar=4, arglag=3)
  df_wb <- prod(unlist(arg_wb))
  
  # TP
  ranges_tp <- t(sapply(datalist,function(x) range(x$tp,na.rm=T)))
  cen_tp <- sapply(datalist,function(x) {
    quantile(x$tp,c(0:100)/100,na.rm=T)
  })
  arg_tp <- list(argvar=4, arglag=3)
  df_tp <- prod(unlist(arg_tp))
  
  # FG10
  ranges_fg10 <- t(sapply(datalist,function(x) range(x$fg10,na.rm=T)))
  cen_fg10 <- sapply(datalist,function(x) {
    quantile(x$fg10,c(0:100)/100,na.rm=T)
  })
  arg_fg10 <- list(argvar=4, arglag=3)
  df_fg10 <- prod(unlist(arg_fg10))
  
  # RH
  ranges_rh <- t(sapply(datalist,function(x) range(x$rh,na.rm=T)))
  cen_rh <- sapply(datalist,function(x) {
    quantile(x$rh,c(0:100)/100,na.rm=T)
  })
  arg_rh <- list(argvar=4, arglag=3)
  df_rh <- prod(unlist(arg_rh))

  # PM25
  ranges_pm25 <- t(sapply(datalist,function(x) range(x$pm25,na.rm=T)))
  cen_pm25 <- sapply(datalist,function(x) {
    quantile(x$pm25,c(0:100)/100,na.rm=T)
  })
  arg_pm25 <- list(argvar=4, arglag=3)
  df_pm25 <- prod(unlist(arg_pm25))
  
  # 2. Stage 1: Fit the DLNM ----
  
  # BUILT OBJECTS WHERE outputs WILL BE STORED:
  
  Slist_comp <- vector('list',m)
  names(Slist_comp) <- locations
  
  ymat_wb <- matrix(NA,m,df_wb,dimnames=list(locations,paste('spl',seq(df_wb),sep='')))
  Slist_wb <- vector('list',m)
  names(Slist_wb) <- locations
  
  ymat_tp <- matrix(NA,m,df_tp,dimnames=list(locations,paste('spl',seq(df_tp),sep='')))
  Slist_tp <- vector('list',m)
  names(Slist_tp) <- locations
  
  ymat_fg10 <- matrix(NA,m,df_fg10,dimnames=list(locations,paste('spl',seq(df_fg10),sep='')))
  Slist_fg10 <- vector('list',m)
  names(Slist_fg10) <- locations
  
  ymat_rh <- matrix(NA,m,df_rh,dimnames=list(locations,paste('spl',seq(df_rh),sep='')))
  Slist_rh <- vector('list',m)
  names(Slist_rh) <- locations
  
  ymat_pm25 <- matrix(NA,m,df_pm25,dimnames=list(locations,paste('spl',seq(df_pm25),sep='')))
  Slist_pm25 <- vector('list',m)
  names(Slist_pm25) <- locations
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  system.time(
    for(i in seq(m)) {
      
      # PRINT ITERATION
      cat(i,'')
      
      # LOAD
      data <- datalist[[i]]
      
      # City specific parameters
      data <- data %>% arrange(date) %>% mutate(patients_lag=lag(total_patients, 1))

      # CREATE THE SPLINES
      cb.wb <- crossbasis(data$wb, lag=c(0,14),argvar=list(fun='ns', df=arg_wb$argvar), arglag=list(fun='ns', df=arg_wb$arglag), group=data$year)
      cb.tp <- crossbasis(data$tp, lag=c(0,14),argvar=list(fun='ns', df=arg_tp$argvar), arglag=list(fun='ns', df=arg_tp$arglag), group=data$year)
      cb.fg10 <- crossbasis(data$fg10, lag=c(0,14),argvar=list(fun='ns', df=arg_fg10$argvar), arglag=list(fun='ns', df=arg_fg10$arglag), group=data$year)
      cb.rh <- crossbasis(data$rh, lag=c(0,14),argvar=list(fun='ns', df=arg_rh$argvar), arglag=list(fun='ns', df=arg_rh$arglag), group=data$year)
      cb.pm25 <- crossbasis(data$pm25, lag=c(0,14),argvar=list(fun='ns', df=arg_pm25$argvar), arglag=list(fun='ns', df=arg_pm25$arglag), group=data$year)
      
      # Set formula
      form <- as.formula(paste0(outcome, ' ~  cb.wb + cb.tp + cb.rh + cb.pm25 + cb.fg10 + 
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
      
      # wb
      predwb <- crosspred(cb.wb, model, cen=cen_wb[reference,i])
      ymat_wb[i,] <- predwb$coef
      Slist_wb[[i]] <- predwb$vcov
      
      # TP
      predtp <- crosspred(cb.tp, model, cen=cen_tp[reference,i])
      ymat_tp[i,] <- predtp$coef
      Slist_tp[[i]] <- predtp$vcov   
      
      # RH
      predrh <- crosspred(cb.rh, model, cen=cen_rh[reference,i])
      ymat_rh[i,] <- predrh$coef
      Slist_rh[[i]] <- predrh$vcov
      
      # FG10
      predfg10 <- crosspred(cb.fg10, model, cen=cen_fg10[reference,i])
      ymat_fg10[i,] <- predfg10$coef
      Slist_fg10[[i]] <- predfg10$vcov
      
      # PM25
      predpm25 <- crosspred(cb.pm25, model, cen=cen_pm25[reference,i])
      ymat_pm25[i,] <- predpm25$coef
      Slist_pm25[[i]] <- predpm25$vcov

      # EXTRACT COMPARATIVE METRICS
      Slist_comp[[i]] <- comp_metrics(model)
      
    })
  
  
  # 3. Stage 2: Multivariate Meta-analysis ----
  
  # Meta coefs
  meta_wb <- get_meta_analysis(exposure='wb', ymat=ymat_wb,  Slist=Slist_wb, arg=arg_wb, ranges=ranges_wb, cen=cen_wb)
  meta_tp <- get_meta_analysis(exposure='tp', ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp)
  meta_fg10 <- get_meta_analysis(exposure='fg10', ymat=ymat_fg10,  Slist=Slist_fg10, arg=arg_fg10, ranges=ranges_fg10, cen=cen_fg10)
  meta_rh <- get_meta_analysis(exposure='rh', ymat=ymat_rh,  Slist=Slist_rh, arg=arg_rh, ranges=ranges_rh, cen=cen_rh)
  meta_pm25 <- get_meta_analysis(exposure='pm25', ymat=ymat_pm25,  Slist=Slist_pm25, arg=arg_pm25, ranges=ranges_pm25, cen=cen_pm25)
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['avg_wb']] <- mean(sapply(Slist_comp, function(x) x$esti_wb))
  meta_comp[['avg_tp']] <- mean(sapply(Slist_comp, function(x) x$esti_tp))
  meta_comp[['avg_rh']] <- mean(sapply(Slist_comp, function(x) x$esti_rh))
  meta_comp[['avg_fg10']] <- mean(sapply(Slist_comp, function(x) x$esti_fg10))
  meta_comp[['avg_pm25']] <- mean(sapply(Slist_comp, function(x) x$esti_pm25))
  meta_comp[['sum_aic']] <- sum(sapply(Slist_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowSums(sapply(Slist_comp, function(x) as.numeric(x$pacf)))
  
  # RESET WARNING
  options(warn=0)
  
  return(list(meta_comp=meta_comp, meta_wb=meta_wb, meta_tp=meta_tp, 
              meta_fg10=meta_fg10, meta_rh=meta_rh, meta_pm25=meta_pm25))
  
}

# Function to extract and format relevant results and figures
format_results <- function(output, outcome, city=NULL){
  
  # Function to get usable tables and figures
  get_results <- function(output, city, lags=NULL, probs='90%', outcome){
    
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
    
    meta_wb <- output$meta_wb
    meta_tp <- output$meta_tp
    meta_rh <- output$meta_rh
    meta_fg10 <- output$meta_fg10
    meta_pm25 <- output$meta_pm25
    
    # Model selection criteria
    tab <- matrix(NA,5,12)
    colnames(tab) <- c('Q','df','p','I-square','AIC','BIC','stat','df','p','stat','df','p')
    rownames(tab) <- c('wb','tp', 'rh', 'fg10', 'pm25')
    tab[1,] <- ftab(m=meta_wb$mv)
    tab[2,] <- ftab(m=meta_tp$mv)
    tab[3,] <- ftab(m=meta_rh$mv)
    tab[4,] <- ftab(m=meta_fg10$mv)
    tab[5,] <- ftab(m=meta_pm25$mv)
    
    # Exposure response (with city specs) - Overall
    p_er_wb <- city_exp_res(city=city, lag=lags, meta=meta_wb, col='#D74B3C', xlab='Temperature (%)', title=NULL)
    p_er_tp <- city_exp_res(city=city, lag=lags, meta=meta_tp, col='#1D98EB', xlab='Rainfall (%)', title=NULL)
    p_er_rh <- city_exp_res(city=city, lag=lags, meta=meta_rh, col='#A6B727', xlab='Relative humidity (%)', title=NULL)
    p_er_fg10 <- city_exp_res(city=city, lag=lags, meta=meta_fg10, col='#8968CD', xlab='Wind speed (%)', title=NULL)
    p_er_pm25 <- city_exp_res(city=city, lag=lags, meta=meta_pm25, col='#E569BF', xlab='Air pollution (%)', title=NULL)
    
    # Lag response (with city specs) by percentile
    q_wb <- rowMeans(meta_wb$cen)[probs]
    q_tp <- rowMeans(meta_tp$cen)[probs]
    q_rh <- rowMeans(meta_rh$cen)[probs]
    q_fg10 <- rowMeans(meta_fg10$cen)[probs]
    q_pm25 <- rowMeans(meta_pm25$cen)[probs]
    p_lr_wb <- city_lag_res(city=city, quant_v=q_wb, meta=meta_wb, col='#D74B3C', xlab='Temperature (%)', title=NULL)
    p_lr_tp <- city_lag_res(city=city, quant_v=q_tp, meta=meta_tp, col='#1D98EB', xlab='Rainfall (%)', title=NULL)
    p_lr_rh <- city_lag_res(city=city, quant_v=q_rh, meta=meta_rh, col='#A6B727', xlab='Relative humidity (%)', title=NULL)
    p_lr_fg10 <- city_lag_res(city=city, quant_v=q_fg10, meta=meta_fg10, col='#8968CD', xlab='Wind speed (%)', title=NULL)
    p_lr_pm25 <- city_lag_res(city=city, quant_v=q_pm25, meta=meta_pm25, col='#E569BF', xlab='Air pollution (%)', title=NULL)
    
    return(list(tab=tab, 
                p_er_wb=p_er_wb, p_er_tp=p_er_tp, p_er_rh=p_er_rh, p_er_fg10=p_er_fg10, p_er_pm25=p_er_pm25,
                p_lr_wb=p_lr_wb, p_lr_tp=p_lr_tp, p_lr_rh=p_lr_rh, p_lr_fg10=p_lr_fg10, p_lr_pm25=p_lr_pm25))
    
  }
  
  overall <- get_results(output, city='ALL', lags=NULL, outcome=outcome)
  all <- overall[!names(overall) %in% 'tab']
  p_blup_all <- lapply(all, function(x) x[['p_blup']])
  p_blup_er <- p_blup_all[grepl('er', names(p_blup_all))] |> head(-2)
  p_blup_lr <- p_blup_all[grepl('lr', names(p_blup_all))] |> head(-2)
  p_er_overall <- wrap_plots(p_blup_er, ncol=3)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
  
  quant_test <- c('10%', '25%', '75%', '90%')
  quant <- list()
  for (q in quant_test){
    quant[[q]] <- get_results(output, city='ALL', lags=NULL, probs=q, outcome=outcome)
  }
  p_lr_wb <- wrap_plots(lapply(quant, function(x) x$p_lr_wb$p_blup), ncol=4) & theme(legend.position='none')
  p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_lr_rh <- wrap_plots(lapply(quant, function(x) x$p_lr_rh$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_lr_fg10 <- wrap_plots(lapply(quant, function(x) x$p_lr_fg10$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_lr_pm25 <- wrap_plots(lapply(quant, function(x) x$p_lr_pm25$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_lr <- wrap_plots(p_lr_wb, p_lr_tp, p_lr_rh, nrow=3)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
  
  lag_test <- c(0, 5, 10, 14)
  lags <- list()
  for (l in lag_test){
    lags[[as.character(l)]] <- get_results(output, city='ALL', lags=l, outcome=outcome)
  }
  p_er_wb <- wrap_plots(lapply(lags, function(x) x$p_er_wb$p_blup), ncol=4) & theme(legend.position='none')
  p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_er_rh <- wrap_plots(lapply(lags, function(x) x$p_er_rh$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_er_fg10 <- wrap_plots(lapply(lags, function(x) x$p_er_fg10$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_er_pm25 <- wrap_plots(lapply(lags, function(x) x$p_er_pm25$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
  p_er <- wrap_plots(p_er_wb, p_er_tp, p_er_rh, nrow=3)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
  
  if(!is.null(city)){
    
    locations <- city
    city <- list()
    city_er <- list()
    city_lr <- list()
    for (i in locations) {
      city_er[[i]] <- list()
      city_lr[[i]] <- list()
      for (l in lag_test) {
        city_er[[i]][[as.character(l)]] <- get_results(output, city=i, lags=l, outcome=outcome)
        print(paste0(i, ' ', l))
      }
      for (q in quant_test) {
        city_lr[[i]][[as.character(q)]] <- get_results(output, city=i, probs=q, outcome=outcome)
        print(paste0(i, ' ', q))
      }
      
      p_er_wb_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_wb$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      p_er_tp_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_tp$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      p_er_rh_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_rh$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      # p_er_fg10_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_fg10$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      # p_er_pm25_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_pm25$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      
      p_lr_wb_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_wb$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      p_lr_tp_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_tp$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      p_lr_rh_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_rh$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      # p_lr_fg10_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_fg10$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      # p_lr_pm25_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_pm25$p_city), ncol=4)+plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
      
    }
  
  } else {
    
    city <- list()
    city_er <- list()
    city_lr <- list()
    for (i in locations){
      city[[i]] <- get_results(output, city=i, probs='50%', outcome=outcome)
    }
    p_er_wb_city <- wrap_plots(lapply(city, function(x) x$p_er_wb$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    p_er_tp_city  <- wrap_plots(lapply(city, function(x) x$p_er_tp$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    p_er_rh_city <- wrap_plots(lapply(city, function(x) x$p_er_rh$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    # p_er_fg10_city  <- wrap_plots(lapply(city, function(x) x$p_er_fg10$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    # p_er_pm25_city  <- wrap_plots(lapply(city, function(x) x$p_er_pm25$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    p_lr_wb_city  <- wrap_plots(lapply(city, function(x) x$p_lr_wb$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    p_lr_tp_city <- wrap_plots(lapply(city, function(x) x$p_lr_tp$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    p_lr_rh_city  <- wrap_plots(lapply(city, function(x) x$p_lr_rh$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    # p_lr_fg10_city <- wrap_plots(lapply(city, function(x) x$p_lr_fg10$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    # p_lr_pm25_city <- wrap_plots(lapply(city, function(x) x$p_lr_pm25$p_city), ncol=4)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
  }

  return(list(overall=overall, p_er_overall=p_er_overall, 
              quant=quant, p_lr=p_lr, lags=lags, p_er=p_er,
              city=city, city_er=city_er, city_lr=city_lr,
              p_er_wb_city=p_er_wb_city, p_er_tp_city=p_er_tp_city, p_er_rh_city=p_er_rh_city, #p_er_fg10_city=p_er_fg10_city, p_er_pm25_city=p_er_pm25_city,
              p_lr_wb_city=p_lr_wb_city, p_lr_tp_city=p_lr_tp_city, p_lr_rh_city=p_lr_rh_city)) #, p_lr_fg10_city=p_lr_fg10_city, p_lr_pm25_city=p_lr_pm25_city))
}

# Function to extract and plot overall & city specific exposure-response relationships
city_exp_res <- function(city='ALL', meta, col, xlab=NULL, lag=NULL, title=NULL, reference='50%'){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col, title=NULL, legend='City-specific') {
    
    if(city=='ALL'){
      
      ggplot() +
        geom_hline(yintercept=1, color='gray20') +
        geom_line(data=df_specific1, aes(x=pred, y=RR, group=city, color=legend, linetype=legend)) +
        geom_ribbon(data=df_overall, aes(x=pred, ymin=low, ymax=high), fill=col, alpha=0.2) +
        geom_line(data=df_overall, aes(x=pred, y=RR, color='Pooled', linetype='Pooled'), size=1.2) +
        scale_color_manual(name=NULL, values=setNames(c('grey70', col), c(legend, 'Pooled'))) +
        scale_linetype_manual(name=NULL, values=setNames(c('dashed', 'solid'), c(legend, 'Pooled'))) +
        scale_x_continuous(breaks=pretty(df_overall$pred, n=5)) +
        scale_y_continuous(breaks=pretty(c(min(df_overall$low), max(df_overall$high)), n=5)) +
        labs(x=xlab, y='RR', title=title) +
        coord_cartesian(ylim=c(min(df_overall$low), max(df_overall$high)), xlim=range(df_overall$pred)) +
        theme_classic() +
        theme(
          axis.line=element_line(color='gray20'),
          axis.text=element_text(color='gray20'),
          axis.ticks=element_line(color='gray20'),
          plot.title=element_text(face='bold', hjust=0.5, size=12),
          legend.position='bottom'
        )
      
    } else{
      df_all <- rbind(
        data.frame(pred=df_overall$pred, RR=df_overall$RR, low=df_overall$low, high=df_overall$high, type='Pooled'),
        data.frame(pred=df_specific1$pred, RR=df_specific1$RR, low=df_specific1$low, high=df_specific1$high, type='Raw'),
        data.frame(pred=df_specific2$pred, RR=df_specific2$RR, low=df_specific2$low, high=df_specific2$high, type='BLUP')
      )
      
      ggplot() +
        geom_hline(yintercept=1, color='gray20') +
        geom_ribbon(data=df_specific1, aes(x=pred, ymin=low, ymax=high), fill=col, alpha=0.2) +
        geom_line(data=df_all, aes(x=pred, y=RR, color=type, linetype=type, size=type)) +
        scale_color_manual(name=NULL, values=c('Pooled'='black', 'Raw'=col, 'BLUP'='black')) +
        scale_linetype_manual(name=NULL, values=c('Pooled'='solid', 'Raw'='solid', 'BLUP'='dashed')) +
        scale_size_manual(name=NULL, values=c('Pooled'=0.7, 'Raw'=1, 'BLUP'=0.7)) +
        scale_x_continuous(breaks=pretty(df_all$pred, n=6)) +
        scale_y_continuous(breaks=pretty(c(min(df_all$low), max(df_all$high)), n=6)) +
        labs(x=xlab, y='RR', title=title) +
        coord_cartesian(ylim=c(min(df_all$low), max(df_all$high)), xlim=range(df_overall$pred)) +
        theme_classic() +
        theme(
          axis.line=element_line(color='gray20'),
          axis.text=element_text(color='gray20'),
          axis.ticks=element_line(color='gray20'),
          plot.title=element_text(face='bold', hjust=0.5, size=12),
          legend.position='bottom'
        )
    }
  }
  
  # Unlist meta-analysis w/ prediction basis
  list2env(meta, envir=environment())
  
  # Remove duplicates from cen list (especially for tp)
  rows_with_dup <- apply(cen, 1, function(row) any(duplicated(row)))
  cen <- cen[!rows_with_dup, ]
  perc <- as.numeric(sub("%", "", names(rowMeans(cen))))
  
  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link='log', at=rowMeans(cen), cen=rowMeans(cen)[reference])
  if(!is.null(lag)){
    df_overall <- data.frame(pred=perc, RR=cp$matRRfit[, lag + 1], low=cp$matRRlow[, lag + 1], high=cp$matRRhigh[, lag + 1], row.names=NULL)
  } else {
    df_overall <- data.frame(pred=perc, RR=cp$allRRfit, low=cp$allRRlow, high=cp$allRRhigh, row.names=NULL)
  }
  
  # if (grepl('Rainfall', xlab)) {
  #   df_overall$pred <- round(exp(df_overall$pred), 2)
  # }
  
  # Get city specific estimates (for specific lag, if applicable)
  df_city <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=ymat[i,], vcov=Slist[[i]], model.link='log', 
                         from=ranges[i, 1], to=ranges[i, 2], at=cen[,i], cen=cen[reference,i])
    if(!is.null(lag)){
      data.frame(city=rownames(ymat)[i], pred=perc, RR=cp_city$matRRfit[, lag + 1], low=cp_city$matRRlow[, lag + 1], high=cp_city$matRRhigh[, lag + 1], row.names=NULL)
    } else {
      data.frame(city=rownames(ymat)[i], pred=perc, RR=cp_city$allRRfit, low=cp_city$allRRlow, high=cp_city$allRRhigh, row.names=NULL)
    }
  })
  
  # Get BLUP estimates (for specific lag, if applicable)
  blup <- blup(mv, vcov=TRUE)
  df_blup <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=blup[[i]]$blup, vcov=blup[[i]]$vcov, model.link='log', 
                         from=ranges[i, 1], to=ranges[i, 2], at=cen[,i], cen=cen[reference,i])
    if(!is.null(lag)){
      data.frame(city=names(blup)[i], pred=perc, RR=cp_city$matRRfit[, lag + 1], low=cp_city$matRRlow[, lag + 1], high=cp_city$matRRhigh[, lag + 1], row.names=NULL)
    } else {
      data.frame(city=names(blup)[i], pred=perc, RR=cp_city$allRRfit, low=cp_city$allRRlow, high=cp_city$allRRhigh, row.names=NULL)
    }
  })
  
  # Get city specific Raw AND BLUP estimates
  if(city=='ALL'){
    
    # Set title based on lags
    if(!is.null(lag)){
      title <- paste0(title, ' (Lag=', lag, ')')
    } 
    
    # Plot average vs city estimates
    p_raw  <- plot_estimates(city=city, df_specific1=df_city, df_overall=df_overall, xlab=xlab, col=col, title=title, legend='City-specific (Raw)')
    # Plot average vs BLUP estimates
    p_blup <- plot_estimates(city=city, df_specific1=df_blup, df_overall=df_overall, xlab=xlab, col=col, title=title, legend='City-specific (BLUP)')
    p_city <- NULL
    
  } else{
    
    # Set title based on lags
    if(!is.null(lag)){
      title <- paste0('Lag=', lag)
    } else {
      title <- city
    }
    df_city <- df_city[df_city$city==city,]
    df_blup <- df_blup[df_blup$city==city,]
    
    # Plot average vs specific city estimates
    p_raw <- NULL
    p_blup <- NULL
    p_city  <- plot_estimates(city=city, df_specific1=df_city, df_specific2=df_blup, df_overall=df_overall, xlab=xlab, col=col, title=title)
    
  }
  
  return(list(df_overall=df_overall, df_city=df_city, df_blup=df_blup,
              p_raw=p_raw, p_blup=p_blup, p_city=p_city))
}

# Function to extract and plot overall & city specific lag-response relationships
city_lag_res <- function(city='ALL', meta, col, xlab=NULL, quant_v=NULL, title=NULL, reference='50%'){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col, title=NULL, legend='City-specific') {
    
    
    if(city=='ALL'){
      
      ylim <- c(min(df_overall$low)-0.01, max(df_overall$high)+0.01)
      
      ggplot() +
        geom_hline(yintercept=1, color='gray20') +
        geom_line(data=df_specific1, aes(x=lag, y=RR, group=city, color=legend, linetype=legend)) +
        geom_line(data=df_overall, aes(x=lag, y=RR, color='Pooled', linetype='Pooled'), size=1.2) +
        geom_errorbar(data=df_overall, aes(x=lag, ymin=low, ymax=high, color='Pooled'), width=0.2, size=0.2) +
        geom_point(data=df_overall, aes(x=lag, y=RR, color='Pooled'), size=2) +
        scale_color_manual(name=NULL, values=setNames(c('grey70', col), c(legend, 'Pooled'))) +
        scale_linetype_manual(name=NULL, values=setNames(c('dashed', 'solid'), c(legend, 'Pooled'))) +
        scale_x_continuous(breaks=seq(min(df_overall$lag), max(df_overall$lag), by=2)) +
        scale_y_continuous(breaks=pretty(ylim, n=6)) +
        coord_cartesian(ylim=ylim, xlim=range(df_overall$lag)) +
        labs(x='Lag (days)', y='RR', title=title) +
        theme_classic() +
        theme(
          axis.line=element_line(color='gray20'),
          axis.text=element_text(color='gray20'),
          axis.ticks=element_line(color='gray20'),
          plot.title=element_text(face='bold', hjust=0.5, size=12),
          legend.position='bottom'
        )
      
    } else{
      
      ylim <- c(min(df_specific1$low)-0.01, max(df_specific1$high)+0.01)
      
      ggplot() +
        geom_hline(yintercept=1, color='gray20') +
        geom_line(data=df_overall, aes(x=lag, y=RR, color='Pooled', linetype='Pooled'), size=0.5) +
        geom_line(data=df_specific1, aes(x=lag, y=RR, group=city, color='Raw', linetype='Raw'), size=0.7) +
        geom_line(data=df_specific2, aes(x=lag, y=RR, group=city, color='BLUP', linetype='BLUP'), size=0.5) +
        geom_errorbar(data=df_specific1, aes(x=lag, ymin=low, ymax=high, color='Raw'), width=0.2, size=0.7) +
        geom_point(data=df_specific1, aes(x=lag, y=RR, color='Raw'), size=2) +
        scale_color_manual(name=NULL, breaks= c('Pooled', 'BLUP', 'Raw'), values=c('Pooled'='black', 'Raw'=col, 'BLUP'='black')) +
        scale_linetype_manual(name=NULL, breaks= c('Pooled', 'BLUP', 'Raw'), values=c('Pooled'='solid', 'Raw'='solid', 'BLUP'='dashed')) +
        scale_x_continuous(breaks=seq(min(df_specific1$lag), max(df_specific1$lag), by=2)) +
        scale_y_continuous(breaks=pretty(ylim, n=6)) +
        coord_cartesian(ylim=ylim, xlim=range(df_specific1$lag)) +
        labs(x='Lag (days)', y='RR', title=title) +
        theme_classic() +
        theme(
          axis.line=element_line(color='gray20'),
          axis.text=element_text(color='gray20'),
          axis.ticks=element_line(color='gray20'),
          plot.title=element_text(face='bold', hjust=0.5, size=12),
          legend.position='bottom'
        )
      
      
    }
  }
  
  # Unlist meta-analysis w/ prediction basis
  list2env(meta, envir=environment())
  
  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link='log', cen=rowMeans(cen)[reference], at=quant_v)
  df_overall <- data.frame(lag=seq(cp$lag[1], cp$lag[2]), RR=t(cp$matRRfit), low=t(cp$matRRlow), high=t(cp$matRRhigh), row.names=NULL)
  names(df_overall) <- c('lag', 'RR', 'low','high')
  
  # Get city specific estimates (for specific lag, if applicable)
  df_city <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=ymat[i,], vcov=Slist[[i]], model.link='log', 
                         from=ranges[i, 1], to=ranges[i, 2], cen=cen[reference,i], at=quant_v)
    data.frame(city=rownames(ymat)[i], lag=seq(cp_city$lag[1], cp_city$lag[2]), RR=t(cp_city$matRRfit), low=t(cp_city$matRRlow), high=t(cp_city$matRRhigh), row.names=NULL)
  })
  names(df_city) <- c('city', 'lag', 'RR', 'low','high')
  
  # Get BLUP estimates (for specific lag, if applicable)
  blup <- blup(mv, vcov=TRUE)
  df_blup <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=blup[[i]]$blup, vcov=blup[[i]]$vcov,
                         model.link='log', from=ranges[i, 1], to=ranges[i, 2], cen=cen[reference,i], at=quant_v)
    data.frame(city=names(blup)[i], lag=seq(cp_city$lag[1], cp_city$lag[2]), RR=t(cp_city$matRRfit), low=t(cp_city$matRRlow), high=t(cp_city$matRRhigh), row.names=NULL)
  })
  names(df_blup) <- c('city', 'lag', 'RR', 'low','high')
  
  # Get city specific Raw AND BLUP estimates
  if(city=='ALL'){
    
    p_raw  <- plot_estimates(city=city, df_specific1=df_city, df_overall=df_overall, xlab=xlab, col=col, title=names(quant_v), legend='City-specific (Raw)')
    p_blup  <- plot_estimates(city=city, df_specific1=df_blup, df_overall=df_overall, xlab=xlab, col=col, title=names(quant_v), legend='City-specific (BLUP)')
    p_city <- NULL
    
  } else{
    
    # Set title based on lags
    if(names(quant_v) != '50%'){
      title <- names(quant_v)
    } else {
      title <- city
    }
    
    df_city <- df_city[df_city$city==city,]
    df_blup <- df_blup[df_blup$city==city,]
    
    # Plot average vs specific city estimates
    p_raw <- NULL
    p_blup <- NULL
    p_city  <- plot_estimates(city=city, df_specific1=df_city, df_specific2=df_blup, df_overall=df_overall, xlab=xlab, col=col, title=title)
    
  }
  
  return(list(df_overall=df_overall, df_city=df_city, df_blup=df_blup,
              p_raw=p_raw, p_blup=p_blup, p_city=p_city))
}

# Function to extract and plot overall & city specific 2D contour plots
city_con_res <- function(meta, xlab=NULL, title=NULL, reference='50%'){
  
  make_dlnm_plot <- function(df, xlab="", title="",
                             ncol=300,   # resolution of color ramp
                             nbreaks=15, # number of breaks for contour lines / ticks
                             contour_size=0.35) {
    

    # Colour palette (same as you used)
    pal <- colorRampPalette(c('#1D98EB', 'white', '#D74B3C'))(ncol)
    
    # Breaks to mimic filled.contour style
    brks <- pretty(range(df$z, na.rm=TRUE), n=nbreaks)
    
    # Base ggplot object
    p <- ggplot(df, aes(x=x, y=y, fill=z)) +
      
      # 1) raster background for smooth fill (interpolated)
      geom_raster(interpolate=TRUE) +
      
      # 2) contour lines coloured according to the same scale (map colour to ..level..)
      geom_contour(
        aes(z=z, colour=..level..),
        breaks=brks,
        size=contour_size,
        show.legend=FALSE    # legend comes from fill/colorbar only
      ) +
      
      # 3) consistent continuous fill scale (this drives the legend)
      scale_fill_gradientn(
        colours=pal,
        limits =range(brks),
        breaks =brks,
        name   ="RR",
        guide  =guide_colorbar(
          barwidth =unit(0.6, "lines"),
          barheight=unit(12,  "lines"),
          ticks=TRUE,
          frame.colour="black",
          ticks.colour="black"
        )
      ) +
      
      # 4) the contour lines must use the same palette so they match fill
      scale_colour_gradientn(
        colours=pal,
        limits =range(brks),
        breaks =brks,
        guide  ="none"   # hide separate colour legend
      ) +
      # optional band-outline trick (draw before the color contour)
      geom_contour(aes(z=z), breaks=brks, colour="grey90", size=0.7, alpha=0.6) +

      # 5) pretty axis breaks & labels (adjust n for desired granularity)
      scale_x_continuous(breaks=pretty(df$x, n=6)) +
      scale_y_continuous(breaks=pretty(df$y, n=6)) +
      
      labs(
        x=xlab,
        y="Lag (days)",
        title=title
      ) +
      
      coord_cartesian(expand=FALSE) +
      
      # 6) theme polish to resemble base filled.contour spacing
      theme_minimal(base_size=13) +
      theme(
        panel.grid=element_blank(),
        axis.title=element_text(size=12),
        axis.text =element_text(size=10),
        legend.title=element_text(size=12),
        legend.text =element_text(size=10),
        plot.title=element_text(size=14, face="plain")
      )
    
    return(p)
  }
  
  
  # Unlist meta-analysis w/ prediction basis
  list2env(meta, envir=environment())
  
  # Remove duplicates from cen list (especially for tp)
  rows_with_dup <- apply(cen, 1, function(row) any(duplicated(row)))
  cen <- cen[!rows_with_dup, ]
  perc <- as.numeric(sub("%", "", names(rowMeans(cen))))
  
  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link='log', at=rowMeans(cen), cen=rowMeans(cen)[reference])
  df_overall <- expand.grid(x=perc, y=seq(cp$lag[1], cp$lag[2]))
  df_overall$z <- as.vector(cp$matRRfit)

  p <- make_dlnm_plot(df_overall, xlab=xlab, title=paste0("Exposure-lag-response: ", title))
 
  return(p)
}


# OTHER ####

# Function to check significant increases by exp levels and get table
meta_checks <- function(meta, threshold){

  cp <- crosspred(meta$b_mean, coef=coef(meta$mv),vcov=vcov(meta$mv), model.link='log', cen=meta$cen)
  (checks <- apply(cp$matRRlow, 2, function(x) rownames(cp$matRRlow)[x > threshold]))
  (table <- as.data.frame(matrix(
    paste0(
      sprintf('%.3f', cp$matRRfit),   # RR with 3 decimals
      ' (',
      sprintf('%.3f', cp$matRRlow),   # lower CI with 3 decimals
      ', ',
      sprintf('%.3f', cp$matRRhigh),  # upper CI with 3 decimals
      ')'
    ),
    nrow=nrow(cp$matRRfit),
    ncol=ncol(cp$matRRfit),
    dimnames=dimnames(cp$matRRfit)  # keep original row/col names
  )))

  return(list(checks=checks, table=table))
}

# Function to extract model diagnostics
get_model_fit <- function(model,name) {
  modelfit <- data.frame(name,
                         summary(model)$df.null+1,
                         summary(model)$df[1],
                         summary(model)$dispersio,
                         summary(model)$deviance,
                         summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion)
  names(modelfit) <- c('model','obs','dfmodel','dispersion','deviance','AIC' )
  return(modelfit)
}

# Function to retrieve results from each predictor, with visuals and tables
get_results <- function(df, outcome, quant_v=0.50, lag_v=NULL){

  df <- df %>%
    # dplyr::group_by(location) %>%
    dplyr::arrange(date, .by_group=TRUE) %>%
    dplyr::mutate(patients_lag=lag(total_patients, 1))

  # Fit cross-basis
  cb.wb <- crossbasis(df$wb, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=3), group=df$year)
  cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=3), group=df$year)
  cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=3), group=df$year)
  cb.rh <- crossbasis(df$rh, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=3), group=df$year)
  cb.pm25 <- crossbasis(df$pm25, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=3), group=df$year)
  
  # Fit model
  model <- glm(formula=as.formula(paste0(outcome, ' ~  season + weekday + holiday + ns(time, df=10*9) + cb.wb + cb.tp + cb.fg10 + cb.rh + cb.pm25 + patients_lag')), data=df, family=quasipoisson(link='log'))

  # Select quantile
  q_wb <- quantile(df$wb, probs=quant_v, na.rm=TRUE)
  q_tp <- quantile(df$tp, probs=quant_v, na.rm=TRUE)
  q_fg10 <- quantile(df$fg10, probs=quant_v, na.rm=TRUE)
  q_rh <- quantile(df$rh, probs=quant_v, na.rm=TRUE)
  q_pm25<- quantile(df$pm25, probs=quant_v, na.rm=TRUE)
  
  # DISPLAY ESTIMATE OF RR FOR BENCHMARK TP INCREMENT OF 23.4 (99th perc)
  cp.wb <- crosspred(cb.wb, model, at=q_wb)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.tp <- crosspred(cb.tp, model, at=q_tp)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.fg10 <- crosspred(cb.fg10, model, at=q_fg10)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.rh <- crosspred(cb.rh, model, at=q_rh)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.pm25<- crosspred(cb.pm25, model, at=q_pm25)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp <- data.frame(rbind(cp.wb, cp.tp, cp.fg10, cp.rh, cp.pm25))
  df_cp <- data.frame(variables=rownames(cp), cp, row.names=NULL)

  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK wb value
  cp <- crosspred(cb.wb, model, at=q_wb)
  p_lr_wb <- get_visuals(cp, model, lag_rep=TRUE, col='#D74B3C', title='Lag-response of wb (at 95th percentile)')
  # GRAPH CUMULATIVE wb-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & FG10)
  cp <- crosspred(cb.wb, model, at=seq(-10, 35, by=1))
  p_er_wb <- get_visuals(cp, model, lag_rep=FALSE, lag=lag_v, col='#D74B3C', xlab='Temperature (C)', title=lag_v)

  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK TP value
  cp <- crosspred(cb.tp, model, at=q_tp)
  p_lr_tp <- get_visuals(cp, model, lag_rep=TRUE, col='#1D98EB', title='Lag-response of tp (at 95th percentile)')
  # GRAPH CUMULATIVE TP-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR wb & FG10)
  cp <- crosspred(cb.tp, model)
  p_er_tp <- get_visuals(cp, model, lag_rep=FALSE, lag=lag_v, col='#1D98EB', xlab='Rainfall (mm)', title=lag_v)

  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK FG10 value
  cp <- crosspred(cb.rh, model, at=q_rh)
  p_lr_rh <- get_visuals(cp, model, lag_rep=TRUE, col='#A6B727', title='Lag-response of fh (at 95th percentile)')
  # GRAPH CUMULATIVE FG10-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & wb)
  cp <- crosspred(cb.rh, model)
  p_er_rh <- get_visuals(cp, model, lag_rep=FALSE, lag=lag_v, col='#A6B727', xlab='Relative humidity (%)', title=lag_v)

  result <- list(model=model, overall=df_cp,
                 exp_rep=list(wb=p_er_wb, tp=p_er_tp, rh=p_er_rh),
                 lag_rep=list(wb=p_lr_wb, tp=p_lr_tp, rh=p_lr_rh))

  return(result)

}

# Function to get cumulative exposure-response and lag-response effects
get_visuals <- function(cp, model, xlab, title, col, lag_rep=FALSE, lag=NULL){


  if(lag_rep==FALSE){
    if(!is.null(lag)){
      cum <- data.frame(x=cp$predvar, RR=cp$matRRfit[, lag + 1], lower=cp$matRRlow[, lag + 1], upper=cp$matRRhigh[, lag + 1], row.names=NULL)
    } else {
      cum <- data.frame(x=cp$predvar, RR=cp$allRRfit, lower=cp$allRRlow, upper=cp$allRRhigh, row.names=NULL)
    }
    p <- ggplot(cum, aes(x=x, y=RR)) +
      geom_line(color=col) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill=col) +
      geom_hline(yintercept=1, color='gray20') +
      scale_x_continuous(breaks=pretty(cum$x, n=6)) +
      # scale_y_continuous(limits=c(-0.5, 2.5)) +
      labs(x=xlab, y='RR',
           title=title) +
      theme_classic() +
      theme(
        axis.line=element_line(color='gray20'),      # pale axis lines
        axis.text=element_text(color='gray20'),      # pale axis labels
        axis.ticks=element_line(color='gray20'),     # pale ticks
        plot.title=element_text(face='bold', hjust=0.5, size=12)  # center title
      )
  } else {

    slice <- data.frame(lag=seq(cp$lag[1], cp$lag[2]), RR=t(cp$matRRfit), lower=t(cp$matRRlow), upper=t(cp$matRRhigh), row.names=NULL)
    names(slice) <- c('lag', 'RR', 'lower','upper')
    p <- ggplot(slice, aes(x=lag, y=RR)) +
      geom_hline(yintercept=1, color='gray20') +
      geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2,  size=0.2) +
      geom_point(color=col,  size=2.5) +
      scale_x_continuous(breaks=seq(min(slice$lag), max(slice$lag), by=2)) +  # more breaks
      labs(x='Lag (days)', y='RR',
           title=title) +
      theme_classic() +
      theme(
        axis.line=element_line(color='gray20'),      # pale axis lines
        axis.text=element_text(color='gray20'),      # pale axis labels
        axis.ticks=element_line(color='gray20'),     # pale ticks
        plot.title=element_text(face='bold', hjust=0.5, size=12)  # center title
      )

  }

  return(p)
}

# Function to format and organizing visuals
get_all_plots <- function(plot_list, single_outcome=NULL, rows=3){

  # Function: Format plots for visualizing (2*4 grid)
  multi_plot_format <- function(plot_list){

    p <- plot_list

    # Loop to remove existing titles and replace with the list name
    p <- lapply(names(p), function(n) {
      p[[n]] + labs(title=n)
    })

    # Keep names after lapply
    names(p) <- names(p)

    # result <- (p[[1]] | p[[2]] | p[[3]] | p[[4]]) / (p[[5]] | p[[6]] | p[[7]] | p[[8]])

    result <- wrap_plots(p, nrow=rows)

    return(result)
  }


  # Visualise all outputs
  if(!is.null(single_outcome)){

    # Arrange your 2x3 grid
    exp_rep <- (plots$exp_rep_wb[[single_outcome]] | plots$exp_rep_tp[[single_outcome]] | plots$exp_rep_fg10[[single_outcome]]) + plot_annotation(title=paste0('Outcome: ', single_outcome), theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    lag_rep <- (plots$lag_rep_wb[[single_outcome]] | plots$lag_rep_tp[[single_outcome]] | plots$lag_rep_fg10[[single_outcome]]) + plot_annotation(title=paste0('Outcome: ', single_outcome), theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    all <- (exp_rep / lag_rep) + plot_annotation(title=paste0('Outcome: ', single_outcome), theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))

    # Add a main title
    grid_plots <- list(all=all, exp_rep=exp_rep, lag_rep=lag_rep)

  } else {

    exp_rep_wb <- multi_plot_format(plot_list$exp_rep_wb) + plot_annotation(title='RR of wb exposure-response (lag=5)', theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    lag_rep_wb <- multi_plot_format(plot_list$lag_rep_wb) + plot_annotation(title='Lag-response of wb (at 95th percentile)', theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    exp_rep_tp <- multi_plot_format(plot_list$exp_rep_tp) + plot_annotation(title='RR of tp exposure-response (lag=6)', theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    lag_rep_tp <- multi_plot_format(plot_list$lag_rep_tp) + plot_annotation(title='Lag-response of tp (at 95th percentile)', theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    exp_rep_fg10 <- multi_plot_format(plot_list$exp_rep_fg10) + plot_annotation(title='Overall RR of fg10 exposure-response', theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))
    lag_rep_fg10 <- multi_plot_format(plot_list$lag_rep_fg10) + plot_annotation(title='Lag-response of fg10 (at 95th percentile)', theme=theme(plot.title=element_text(size=20, hjust=0.5, face='bold')))

    grid_plots <- list(exp_rep_wb=exp_rep_wb, lag_rep_wb=lag_rep_wb, exp_rep_tp=exp_rep_tp,
                       lag_rep_tp=lag_rep_tp, exp_rep_fg10=exp_rep_fg10, lag_rep_fg10=lag_rep_fg10)

  }

  return(grid_plots)

}

#### DEPENDENCIES ####

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD','#52C5B0FF', '#910030FF', '#127804FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9','#52C5B080', '#91003080', '#12780480') 

w_cols <- c('wb', 'mn2t','mx2t','tp','fg10', 'rh_mean', 'pm25')
h_cols <- c('total_patients',
            'F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
            'gynaecology', 'orthopaedics', 'dermatology', 'surgery',
            'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
            'ncd', 'antiretroviral', 'substance_counselling', 'physiotherapy')

location_levels <- c('Bilaspur','Chamba', 'Dharamshala', 'Hamirpur', 'Kullu', 'Kyelang', 
                     'Mandi', 'Nahan', 'Reckong Peo', 'Shimla', 'Solan', 'Una')
season_levels   <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
month_levels    <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
weekday_levels  <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
age_levels      <- c('under_5', 'age5_14', 'age15_24', 'age25_64',  'over_65','Unknown')
gender_levels   <- c('F', 'M', 'O')
ward_levels     <- c('general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
                     'gynaecology', 'orthopaedics', 'dermatology', 'surgery',
                     'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
                     'ncd', 'mch', 'antiretroviral', 'substance_counselling', 'physiotherapy')

# Dictionary of names for figures and tables
nm <- c(season='Season', closed='Closed', holiday='Holiday', total_patients='Total patients', 
        elev='Elevation', popu='Population density',
        # POPULATION GROUPS
        M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        # WARDS
        general='General', emergency='Emergency', ophthalmology='Ophthalmology', 
        paediatrics='Paediatric', orthopaedics='Orthopaedics', ent='ENT', 
        gynaecology='Gynaecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        pulmonary='Pulmonary', disability='Disability', ncd='NCD', mch='MCH',
        antiretroviral='HIV Treatment', substance_counselling='Substance counselling',
        physiotherapy='Physiotherapy', 
        # CONTINUOUS WEATHER
        wb='Mean temp. (\u00B0C)', mx2t='Max temp. (\u00B0C)', mn2t='Min temp. (\u00B0C)',
        rh='Relative humidity (%)', tp='Rainfall (mm)', fg10='Windspeed (km/hr)', pm25='Air pollution (µg/m³)',
        dt="Diurnal temp. range (\u00B0C)", wb="Wet-bulb temp. (\u00B0C)",
        # IMD DEFINITIONS
        mx2t_IMD='IMD: Mx temp. >30\u00B0C', mn2t_IMD='IMD: Min temp. <0\u00B0C', 
        tp_IMD='IMD: Rainfall >15mm', fg10_IMD='IMD: Wind speed >40km/h', pm25_IMD='IMD: PM2.5 >60 µg/m³', 
        rh_low_IMD='IMD: Humidity <30%', rh_hig_IMD='IMD: Humidity >60%',
        # CLIM. DEFINITIONS
        mx2t_p90='Clim: Max temp. >90th', mn2t_p10='Clim: Min temp. <10th', 
        tp_p90='Clim: Rainfall >90th', fg10_p90='Clim: Wind speed >90th', pm25_p90='Clim: PM2.5 >90th', 
        rh_mean_p10 ='Clim: Humidity <10th', rh_mean_p90='Clim: Humidity >90th') 

