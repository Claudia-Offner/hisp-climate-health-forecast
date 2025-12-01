### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city meta-analysis                                         ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta)
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

#### FUNCTIONS ####

# Function to run model with static formula
run_model <- function(data, locations, outcome){
  
  # Function to get reference meta-analysis
  get_meta_analysis <- function(ymat, Slist, arg, ranges, cen){
    
    # Meta-analysis w/ prediction basis
    mv <- mvmeta(ymat, Slist, method="ml")
    bound <- c(colMeans(ranges))
    mean <- seq(bound[1],bound[2],length=30)
    b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
    
    return(list(ymat=ymat, Slist=Slist, arg=arg, ranges=ranges, mv=mv, bound=bound, b_mean=b_mean, cen=cen))
    
  }
  
  # 1. Data and dependency prep ----
  # LIST OF DATAFRAMES FOR ALL SITES
  datalist <- lapply(locations, function(region) data[data$location==region,])
  names(datalist) <- locations
  
  # ADDITIONAL CITY-LEVEL META-PREDICTORS
  lat <- unique(data$lat)
  perclat <- round(quantile(lat,c(1,3)*0.25),1)
  
  # T2M
  ranges_t2m <- t(sapply(datalist,function(x) range(x$t2m,na.rm=T)))
  cen_t2m <- mean(data$t2m)
  arg_t2m <- list(argvar=3, arglag=4)
  df_t2m <- prod(unlist(arg_t2m))
  
  # TP
  ranges_tp <- t(sapply(datalist,function(x) range(x$tp,na.rm=T)))
  cen_tp <- mean(data$tp)
  arg_tp <- list(argvar=3, arglag=4)
  df_tp <- prod(unlist(arg_tp))
  
  # RH
  ranges_rh <- t(sapply(datalist,function(x) range(x$rh,na.rm=T)))
  cen_rh <- mean(data$rh)
  arg_rh <- list(argvar=3, arglag=4)
  df_rh <- prod(unlist(arg_rh))
  
  # PM25
  ranges_pm25 <- t(sapply(datalist,function(x) range(x$pm25,na.rm=T)))
  cen_pm25 <- mean(data$pm25)
  arg_pm25 <- list(argvar=3, arglag=4)
  df_pm25 <- prod(unlist(arg_pm25))
  
  # # DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
  # knotperc <- c(10,30,70,90)
  # knots <- rowMeans(sapply(datalist,function(x) 
  #   quantile(x$t2m,knotperc/100,na.rm=T)))
  
  # 2. Stage 1: Fit the DLNM ----
  
  # BUILT OBJECTS WHERE outputs WILL BE STORED:
  ymat_t2m <- matrix(NA,m,df_t2m,dimnames=list(locations,paste("spl",seq(df_t2m),sep="")))
  Slist_t2m <- vector("list",m)
  names(Slist_t2m) <- locations
  
  ymat_tp <- matrix(NA,m,df_tp,dimnames=list(locations,paste("spl",seq(df_tp),sep="")))
  Slist_tp <- vector("list",m)
  names(Slist_tp) <- locations
  
  ymat_rh <- matrix(NA,m,df_rh,dimnames=list(locations,paste("spl",seq(df_rh),sep="")))
  Slist_rh <- vector("list",m)
  names(Slist_rh) <- locations
  
  ymat_pm25 <- matrix(NA,m,df_pm25,dimnames=list(locations,paste("spl",seq(df_pm25),sep="")))
  Slist_pm25 <- vector("list",m)
  names(Slist_pm25) <- locations
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  system.time(
    for(i in seq(m)) {
      
      # PRINT ITERATION
      cat(i,"")
      
      # LOAD
      data <- datalist[[i]]
      
      # CREATE THE SPLINES
      cb.t2m <- crossbasis(data$t2m, lag=c(0,14),argvar=list(fun='ns', df=arg_t2m$argvar), arglag=list(fun='ns', df=arg_t2m$arglag), group=data$year)
      cb.tp <- crossbasis(data$tp, lag=c(0,14),argvar=list(fun='ns', df=arg_tp$argvar), arglag=list(fun='ns', df=arg_tp$arglag), group=data$year)
      cb.rh <- crossbasis(data$rh, lag=c(0,14),argvar=list(fun='ns', df=arg_rh$argvar), arglag=list(fun='ns', df=arg_rh$arglag), group=data$year)
      cb.pm25 <- crossbasis(data$pm25, lag=c(0,14),argvar=list(fun='ns', df=arg_pm25$argvar), arglag=list(fun='ns', df=arg_pm25$arglag), group=data$year)
      
      # Set formula
      form <- as.formula(paste0(outcome, ' ~  season + weekday + closed + ns(time, df=7*length(unique(data$year))) + cb.t2m + cb.tp + cb.rh + cb.pm25'))
      
      # RUN THE MODEL
      model <- glm(formula=form,family=quasipoisson(),data)
      
      # EXTRACT AND SAVE THE RELATED COEF AND VCOV
      
      # T2M
      predt2m <- crosspred(cb.t2m, model, cen=cen_t2m)
      ymat_t2m[i,] <- predt2m$coef
      Slist_t2m[[i]] <- predt2m$vcov
      
      # TP
      predtp <- crosspred(cb.tp, model, cen=cen_tp)
      ymat_tp[i,] <- predtp$coef
      Slist_tp[[i]] <- predtp$vcov   
      
      # RH
      predrh <- crosspred(cb.rh, model, cen=cen_rh)
      ymat_rh[i,] <- predrh$coef
      Slist_rh[[i]] <- predrh$vcov  
      
      # PM25
      predpm25 <- crosspred(cb.pm25, model, cen=cen_pm25)
      ymat_pm25[i,] <- predpm25$coef
      Slist_pm25[[i]] <- predpm25$vcov  
      
    })

  
  # 3. Stage 2: Multivariate Meta-analysis ----
  
  meta_t2m <- get_meta_analysis(ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m)
  meta_tp <- get_meta_analysis(ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp)
  meta_rh <- get_meta_analysis(ymat=ymat_rh,  Slist=Slist_rh, arg=arg_rh, ranges=ranges_rh, cen=cen_rh)
  meta_pm25 <- get_meta_analysis(ymat=ymat_pm25,  Slist=Slist_pm25, arg=arg_pm25, ranges=ranges_pm25, cen=cen_pm25)
  
  # RESET WARNING
  options(warn=0)
  
  return(list(meta_t2m=meta_t2m, meta_tp=meta_tp, meta_rh=meta_rh, meta_pm25=meta_pm25))
  
}

# Function to extract and format relevant results and figures
format_results <- function(output, outcome){
  
  overall <- get_results(output, city='ALL', lags=NULL, outcome=outcome)
  all <- overall[!names(overall) %in% "tab"]
  p_blup_all <- lapply(all, function(x) x[["p_blup"]])
  p_blup_er <- p_blup_all[grepl("er", names(p_blup_all))]
  p_blup_lr <- p_blup_all[grepl("lr", names(p_blup_all))]
  p_er_overall <- wrap_plots(p_blup_er, ncol=4)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  quant_test <- c(0.1, 0.25, 0.75, 0.9)
  quant <- list()
  for (q in quant_test){
    quant[[as.character(q)]] <- get_results(output, city='ALL', lags=NULL, probs=q, outcome=outcome)
  }
  p_lr_t2m <- wrap_plots(lapply(quant, function(x) x$p_lr_t2m$p_blup), ncol=4) & theme(legend.position = "none")
  p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_blup), ncol=4) & theme(legend.position = "none", plot.title = element_blank())
  p_lr_rh <- wrap_plots(lapply(quant, function(x) x$p_lr_rh$p_blup), ncol=4) & theme(legend.position = "none", plot.title = element_blank())
  p_lr_pm25 <- wrap_plots(lapply(quant, function(x) x$p_lr_pm25$p_blup), ncol = 4) & theme(legend.position = "none", plot.title = element_blank())
  p_lr <- wrap_plots(p_lr_t2m, p_lr_tp, p_lr_rh, p_lr_pm25, nrow=4)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  lag_test <- c(0, 5, 10, 14)
  lags <- list()
  for (l in lag_test){
    lags[[as.character(l)]] <- get_results(output, city='ALL', lags=l, outcome=outcome)
  }
  p_er_t2m <- wrap_plots(lapply(lags, function(x) x$p_er_t2m$p_blup), ncol=4) & theme(legend.position = "none")
  p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_blup), ncol=4) & theme(legend.position = "none", plot.title = element_blank())
  p_er_rh <- wrap_plots(lapply(lags, function(x) x$p_er_rh$p_blup), ncol=4) & theme(legend.position = "none", plot.title = element_blank())
  p_er_pm25 <- wrap_plots(lapply(lags, function(x) x$p_er_pm25$p_blup), ncol=4) & theme(legend.position = "none", plot.title = element_blank())
  p_er <- wrap_plots(p_er_t2m, p_er_tp, p_er_rh, p_er_pm25, nrow=4)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  city <- list()
  for (i in locations){
    city[[i]] <- get_results(output, city=i, lags=NULL, outcome=outcome)
  }
  p_er_t2m_city <- wrap_plots(lapply(city, function(x) x$p_er_t2m$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_er_tp_city  <- wrap_plots(lapply(city, function(x) x$p_er_tp$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_er_rh_city <- wrap_plots(lapply(city, function(x) x$p_er_rh$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_er_pm25_city  <- wrap_plots(lapply(city, function(x) x$p_er_pm25$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  p_lr_t2m_city  <- wrap_plots(lapply(city, function(x) x$p_lr_t2m$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_lr_tp_city <- wrap_plots(lapply(city, function(x) x$p_lr_tp$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_lr_rh_city  <- wrap_plots(lapply(city, function(x) x$p_lr_rh$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_lr_pm25_city <- wrap_plots(lapply(city, function(x) x$p_lr_pm25$p_city), ncol=3)+plot_annotation(title=outcome, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  return(list(overall=overall, p_er_overall=p_er_overall, 
              quant=quant, p_lr=p_lr, lags=lags, p_er=p_er,
              city=city, 
              p_er_t2m_city=p_er_t2m_city, p_er_tp_city=p_er_tp_city, p_er_rh_city=p_er_rh_city, p_er_pm25_city=p_er_pm25_city,
              p_lr_t2m_city=p_lr_t2m_city, p_lr_tp_city=p_lr_tp_city, p_lr_rh_city=p_lr_rh_city, p_lr_pm25_city=p_lr_pm25_city))
}

# Function to get usable tables and figures
get_results <- function(output, city, lags, probs=0.5, outcome){
  
  meta_t2m <- output$meta_t2m
  meta_tp <- output$meta_tp
  meta_rh <- output$meta_rh
  meta_pm25 <- output$meta_pm25
  
  # Model selection criteria
  tab <- matrix(NA,4,12)
  colnames(tab) <- c("Q","df","p","I-square","AIC","BIC","stat","df","p","stat","df","p")
  rownames(tab) <- c("t2m","tp", 'rh', 'pm25')
  tab[1,] <- ftab(m=meta_t2m$mv)
  tab[2,] <- ftab(m=meta_tp$mv)
  tab[3,] <- ftab(m=meta_rh$mv)
  tab[4,] <- ftab(m=meta_pm25$mv)
  
  # Exposure response (with city specs) - Overall
  p_er_t2m <- city_exp_res(city=city, lag=lags, meta=meta_t2m, col="#D74B3C", xlab="Temperature (°C)", title=NULL)
  p_er_tp <- city_exp_res(city=city, lag=lags, meta=meta_tp, col="#1D98EB", xlab="Rainfall (mm)", title=NULL)
  p_er_rh <- city_exp_res(city=city, lag=lags, meta=meta_rh, col="#A6B727", xlab="Relative humidity (%)", title=NULL)
  p_er_pm25 <- city_exp_res(city=city, lag=lags, meta=meta_pm25, col="#8968CD", xlab="Air pollution (µg/m³)", title=NULL)
  
  # Lag response (with city specs) - e.g. 95th percentile
  q_t2m <- quantile(df$t2m, probs=probs, na.rm = TRUE)
  q_tp <- quantile(df$tp, probs=probs, na.rm = TRUE)
  q_rh <- quantile(df$rh, probs=probs, na.rm = TRUE)
  q_pm25 <- quantile(df$pm25, probs=probs, na.rm = TRUE)
  p_lr_t2m <- city_lag_res(city=city, quant_v=q_t2m, meta=meta_t2m, col="#D74B3C", xlab="Temperature (°C)", title=NULL)
  p_lr_tp <- city_lag_res(city=city, quant_v=q_tp, meta=meta_tp, col="#1D98EB", xlab="Rainfall (mm)", title=NULL)
  p_lr_rh <- city_lag_res(city=city, quant_v=q_rh, meta=meta_rh, col="#A6B727", xlab="Relative humidity (%)", title=NULL)
  p_lr_pm25 <- city_lag_res(city=city, quant_v=q_pm25, meta=meta_pm25, col="#8968CD", xlab="Air pollution (µg/m³)", title=NULL)
  
  return(list(tab=tab, 
              p_er_t2m=p_er_t2m, p_er_tp=p_er_tp, p_er_rh=p_er_rh, p_er_pm25=p_er_pm25, 
              p_lr_t2m=p_lr_t2m, p_lr_tp=p_lr_tp, p_lr_rh=p_lr_rh, p_lr_pm25=p_lr_pm25))
  
}

# Function to extract model comparison criteria
ftab <- function(m, mref=NULL) {
  
  # HETEROGENEITY AND IC STATS
  q <- qtest(m)
  het <- c(q$Q[1],q$df[1],q$pvalue[1],(q$Q[1]-q$df[1])/q$Q[1]*100,AIC(m),BIC(m))
  # LR TEST (ONLY FOR META-REGRESSION)
  if(!is.null(mref)) {
    lrstat <- -2*(logLik(mref)-logLik(m))[1]
    df <- attr(logLik(m),"df")-attr(logLik(mref),"df")
    pvalue <- 1-pchisq(lrstat,df)
    lr <- c(lrstat,df,pvalue)
  }
  # WALD TEST (ONLY FOR META-REGRESSION)
  if(!is.null(mref)) {
    coef <- coef(m)[-grep("Int",names(coef(m)))]
    vcov <- vcov(m)[-grep("Int",names(coef(m))),-grep("Int",names(coef(m)))]
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

# Function to extract and plot overall & city specific exposure-response relationships
city_exp_res <- function(city='ALL', meta, col, xlab=NULL, lag=NULL, title=NULL){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col, title=NULL, legend='City-specific') {
    
    if(city=='ALL'){
      
      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_line(data = df_specific1, aes(x = pred, y = RR, group = city, color = legend, linetype = legend)) +
        geom_ribbon(data = df_overall, aes(x = pred, ymin = low, ymax = high), fill = col, alpha = 0.2) +
        geom_line(data = df_overall, aes(x = pred, y = RR, color = "Pooled", linetype = "Pooled"), size = 1.2) +
        scale_color_manual(name = NULL, values = setNames(c("grey70", col), c(legend, "Pooled"))) +
        scale_linetype_manual(name = NULL, values = setNames(c("dashed", "solid"), c(legend, "Pooled"))) +
        scale_x_continuous(breaks = pretty(df_overall$pred, n = 6)) +
        scale_y_continuous(breaks = pretty(c(min(df_overall$low), max(df_overall$high)), n = 6)) +
        labs(x = xlab, y = "RR", title = title) +
        coord_cartesian(ylim=c(min(df_overall$low), max(df_overall$high)), xlim = range(df_overall$pred)) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "bottom"
        )
      
    } else{
      df_all <- rbind(
        data.frame(pred = df_overall$pred, RR = df_overall$RR, low = df_overall$low, high = df_overall$high, type = "Pooled"),
        data.frame(pred = df_specific1$pred, RR = df_specific1$RR, low = df_specific1$low, high = df_specific1$high, type = "Raw"),
        data.frame(pred = df_specific2$pred, RR = df_specific2$RR, low = df_specific2$low, high = df_specific2$high, type = "BLUP")
      )
      
      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_ribbon(data = df_specific2, aes(x = pred, ymin = low, ymax = high), fill = col, alpha = 0.2) +
        geom_line(data = df_all, aes(x = pred, y = RR, color = type, linetype = type, size = type)) +
        scale_color_manual(name = NULL, values = c("Pooled" = "black", "Raw" = "black", "BLUP" = col)) +
        scale_linetype_manual(name = NULL, values = c("Pooled" = "solid", "Raw" = "dashed", "BLUP" = "solid")) +
        scale_size_manual(name = NULL, values = c("Pooled" = 1.2, "Raw" = 0.7, "BLUP" = 0.7)) +
        scale_x_continuous(breaks = pretty(df_all$pred, n = 6)) +
        scale_y_continuous(breaks = pretty(c(min(df_all$low), max(df_all$high)), n = 6)) +
        labs(x = xlab, y = "RR", title = title) +
        coord_cartesian(ylim = c(min(df_all$low), max(df_all$high)), xlim = range(df_overall$pred)) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "bottom"
        )
    }
  }
  
  # Unlist meta-analysis w/ prediction basis
  list2env(meta, envir = environment())

  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link="log", cen=cen)
  if(!is.null(lag)){
    df_overall <- data.frame(pred=cp$predvar, RR=cp$matRRfit[, lag + 1], low=cp$matRRlow[, lag + 1], high=cp$matRRhigh[, lag + 1], row.names = NULL)
  } else {
    df_overall <- data.frame(pred=cp$predvar, RR=cp$allRRfit, low=cp$allRRlow, high=cp$allRRhigh, row.names = NULL)
  }
  
  # Get city specific estimates (for specific lag, if applicable)
  df_city <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef = ymat[i,], vcov = Slist[[i]], model.link = "log", 
                         from = ranges[i, 1], to = ranges[i, 2], cen = cen)
    if(!is.null(lag)){
    data.frame(city=rownames(ymat)[i], pred=cp_city$predvar, RR=cp_city$matRRfit[, lag + 1], low=cp_city$matRRlow[, lag + 1], high=cp_city$matRRhigh[, lag + 1], row.names = NULL)
    } else {
    data.frame(city=rownames(ymat)[i], pred=cp_city$predvar, RR=cp_city$allRRfit, low=cp_city$allRRlow, high=cp_city$allRRhigh, row.names = NULL)
    }
  })
  
  # Get BLUP estimates (for specific lag, if applicable)
  blup <- blup(mv, vcov=TRUE)
  df_blup <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef = blup[[i]]$blup, vcov = blup[[i]]$vcov,
                         model.link = "log", from = ranges[i, 1], to = ranges[i, 2], cen = cen)
    if(!is.null(lag)){
      data.frame(city=names(blup)[i], pred=cp_city$predvar, RR=cp_city$matRRfit[, lag + 1], low=cp_city$matRRlow[, lag + 1], high=cp_city$matRRhigh[, lag + 1], row.names = NULL)
    } else {
      data.frame(city=names(blup)[i], pred=cp_city$predvar, RR=cp_city$allRRfit, low=cp_city$allRRlow, high=cp_city$allRRhigh, row.names = NULL)
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
      title <- paste0(city, ' (Lag=', lag, ')')
    } else {
      title <- paste0(city)
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
city_lag_res <- function(city='ALL', meta, col, xlab=NULL, quant_v=NULL, title=NULL){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col, legend='City-specific') {

    
    if(city=='ALL'){

    ylim <- c(min(df_overall$low)-0.01, max(df_overall$high)+0.01)

      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_line(data = df_specific1, aes(x = lag, y = RR, group = city, color = legend, linetype = legend)) +
        geom_line(data = df_overall, aes(x = lag, y = RR, color = "Pooled", linetype = "Pooled"), size = 1.2) +
        geom_errorbar(data = df_overall, aes(x = lag, ymin = low, ymax = high, color = "Pooled"), width = 0.2, size = 0.2) +
        geom_point(data = df_overall, aes(x = lag, y = RR, color = "Pooled"), size=2) +
        scale_color_manual(name = NULL, values = setNames(c("grey70", col), c(legend, "Pooled"))) +
        scale_linetype_manual(name = NULL, values = setNames(c("dashed", "solid"), c(legend, "Pooled"))) +
        scale_x_continuous(breaks = seq(min(df_overall$lag), max(df_overall$lag), by = 2)) +
        scale_y_continuous(breaks = pretty(ylim, n = 6)) +
        coord_cartesian(ylim=ylim, xlim = range(df_overall$lag)) +
        labs(x = "Lag (days)", y = "RR", title = names(quant_v)) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "bottom"
        )
      
    } else{

      ylim <- c(min(df_specific2$low)-0.01, max(df_specific2$high)+0.01)
      
      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_line(data = df_overall, aes(x = lag, y = RR, color = "Pooled", linetype = "Pooled"), size = 1.2) +
        geom_line(data = df_specific1, aes(x = lag, y = RR, group = city, color = "Raw", linetype = "Raw")) +
        geom_line(data = df_specific2, aes(x = lag, y = RR, group = city, color = "BLUP", linetype = "BLUP"), size = 1) +
        geom_errorbar(data = df_specific2, aes(x = lag, ymin = low, ymax = high, color = "BLUP"), width = 0.2, size = 0.2) +
        geom_point(data = df_specific2, aes(x = lag, y = RR, color = "BLUP"), size=2) +
        scale_color_manual(name = NULL, values = c("Pooled" = "black", "Raw" = col, "BLUP" = col)) +
        scale_linetype_manual(name = NULL, values = c("Pooled" = "solid", "Raw" = "dashed", "BLUP" = "solid")) +
        scale_x_continuous(breaks = seq(min(df_overall$lag), max(df_overall$lag), by = 2)) +
        scale_y_continuous(breaks = pretty(ylim, n = 6)) +
        coord_cartesian(ylim=ylim, xlim = range(df_overall$lag)) +
        labs(x = "Lag (days)", y = "RR", title = paste0(city, " - ", names(quant_v))) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "bottom"
        )
      
      
    }
  }
  
  # Unlist meta-analysis w/ prediction basis
  list2env(meta, envir = environment())
  
  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link="log", cen=cen, at=quant_v)
  df_overall <- data.frame(lag=seq(cp$lag[1], cp$lag[2]), RR=t(cp$matRRfit), low=t(cp$matRRlow), high=t(cp$matRRhigh), row.names = NULL)
  names(df_overall) <- c('lag', 'RR', 'low','high')

  # Get city specific estimates (for specific lag, if applicable)
  df_city <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=ymat[i,], vcov=Slist[[i]], model.link="log", 
                         from=ranges[i, 1], to=ranges[i, 2], cen=cen, at=quant_v)
    data.frame(city=rownames(ymat)[i], lag=seq(cp_city$lag[1], cp_city$lag[2]), RR=t(cp_city$matRRfit), low=t(cp_city$matRRlow), high=t(cp_city$matRRhigh), row.names = NULL)
  })
  names(df_city) <- c('city', 'lag', 'RR', 'low','high')
  
  # Get BLUP estimates (for specific lag, if applicable)
  blup <- blup(mv, vcov=TRUE)
  df_blup <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=blup[[i]]$blup, vcov=blup[[i]]$vcov,
                         model.link="log", from=ranges[i, 1], to=ranges[i, 2], cen=cen, at=quant_v)
    data.frame(city=names(blup)[i], lag=seq(cp_city$lag[1], cp_city$lag[2]), RR=t(cp_city$matRRfit), low=t(cp_city$matRRlow), high=t(cp_city$matRRhigh), row.names = NULL)
  })
  names(df_blup) <- c('city', 'lag', 'RR', 'low','high')
  
  # Get city specific Raw AND BLUP estimates
  if(city=='ALL'){
    
    p_raw  <- plot_estimates(city=city, df_specific1=df_city, df_overall=df_overall, xlab=xlab, col=col, legend='City-specific (Raw)')
    p_blup  <- plot_estimates(city=city, df_specific1=df_blup, df_overall=df_overall, xlab=xlab, col=col, legend='City-specific (BLUP)')
    p_city <- NULL
    
  } else{
    
    df_city <- df_city[df_city$city==city,]
    df_blup <- df_blup[df_blup$city==city,]
    
    # Plot average vs specific city estimates
    p_raw <- NULL
    p_blup <- NULL
    p_city  <- plot_estimates(city=city, df_specific1=df_city, df_specific2=df_blup, df_overall=df_overall, xlab=xlab, col=col)
    
  }
  
  return(list(df_overall=df_overall, df_city=df_city, df_blup=df_blup,
              p_raw=p_raw, p_blup=p_blup, p_city=p_city))
}


#### MAIN CODE ####

load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# Run model ----

# Total patients
primary <- 'total_patients'
output <- run_model(df, locations, primary)
main_res <- format_results(output, primary)
total_er_overall <- main_res$p_er_overall + plot_annotation(title = replace_exact_matches(primary, nm))
total_lr <- main_res$p_lr + plot_annotation(title = replace_exact_matches(primary, nm))
total_er <- main_res$p_er + plot_annotation(title = replace_exact_matches(primary, nm))

# Total patients, by city
city_t2m_er <- main_res$p_er_t2m_city + plot_annotation(title = replace_exact_matches(primary, nm))
city_tp_er <- main_res$p_er_tp_city + plot_annotation(title = replace_exact_matches(primary, nm))
city_rh_er <- main_res$p_er_rh_city + plot_annotation(title = replace_exact_matches(primary, nm))
city_pm25_er <- main_res$p_er_pm25_city + plot_annotation(title = replace_exact_matches(primary, nm))

all_cities <- list(city_t2m_er=city_t2m_er, city_tp_er=city_tp_er, 
                  city_rh_er=city_rh_er, city_pm25_er=city_pm25_er)

# Multiple outcomes
outcomes <- h_cols[1:8]
er5_t2m <- list()
er10_t2m <- list()
er0_tp <- list()
er14_rh <- list()
er14_pm25 <- list()

for (i in outcomes){
  output <- run_model(df, locations, i)
  res <- format_results(output, i)
  er5_t2m[[i]] <- res$lags$`5`$p_er_t2m$p_blup + labs(title=replace_exact_matches(i, nm))
  er10_t2m[[i]] <- res$lags$`10`$p_er_t2m$p_blup + labs(title=replace_exact_matches(i, nm))
  er0_tp[[i]] <- res$lags$`0`$p_er_tp$p_blup + labs(title=replace_exact_matches(i, nm))
  er14_rh[[i]] <- res$lags$`14`$p_er_rh$p_blup + labs(title=replace_exact_matches(i, nm))
  er14_pm25[[i]] <- res$lags$`14`$p_er_pm25$p_blup + labs(title=replace_exact_matches(i, nm))
}

er5_t2m <- (wrap_plots(er5_t2m, ncol = 4)& theme(legend.position = "none"))+plot_annotation(title='Expsoure-response (Lag=5)', theme = theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
er10_t2m <- (wrap_plots(er10_t2m, ncol = 4)& theme(legend.position = "none"))+plot_annotation(title='Expsoure-response (Lag=10)', theme = theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
er0_tp <- (wrap_plots(er0_tp, ncol = 4)& theme(legend.position = "none"))+plot_annotation(title='Expsoure-response (Lag=0)', theme = theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
er14_rh <- (wrap_plots(er14_rh, ncol = 4)& theme(legend.position = "none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme = theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
er14_pm25 <- (wrap_plots(er14_pm25, ncol = 4)& theme(legend.position = "none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme = theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))

all_outcomes <- list(er5_t2m=er5_t2m, er10_t2m=er10_t2m, 
                     er0_tp=er0_tp, er14_rh=er14_rh, er14_pm25=er14_pm25)

# Specific location

#### EXPORT #####
width <- 10
height <- 6

# Save cumulative exposure-response (total_patient)
ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_total_er_overall.png'), plot = total_er_overall, width=10, height=5, dpi = 300)
ggsave(filename = paste0(path_report, '- figures/MF_total_er_overall.png'), plot = total_er_overall, width=10, height=5, dpi = 300)

# Save lag-response by percentile (total_patient)
ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_total_lr.png'), plot = total_lr, width=width, height=height, dpi = 300)
ggsave(filename = paste0(path_report, '- figures/MF_total_lr.png'), plot = total_lr, width=width, height=height, dpi = 300)

# Save exposure-response by lag (total_patient)
ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_total_er.png'), plot = total_er, width=width, height=height, dpi = 300)
ggsave(filename = paste0(path_report, '- figures/MF_total_er.png'), plot = total_er, width=width, height=height, dpi = 300)

# Save cumulative exposure-response by city (total_patient)
for (n in names(all_cities)) {
  ggsave(filename = paste0(path_ch2, '/IV. Figures/SF_', n, '.png'), plot = all_cities[[n]], width=10, height=8, dpi = 300)
  ggsave(filename = paste0(path_report, '- figures/SF_', n, '.png'), plot = all_cities[[n]], width=10, height=6, dpi = 300)
}

# Population plots
names(all_outcomes) <- paste0("disag_", names(all_outcomes))
for (n in names(all_outcomes)) {
  ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_', n, '.png'), plot = all_outcomes[[n]], width=10, height=8, dpi = 300)
  ggsave(filename = paste0(path_report, '- figures/MF_', n, '.png'), plot = all_outcomes[[n]], width=10, height=6, dpi = 300)
}

################################################################################
# FIG: Meta-regression outputs (e.g. elevation) ----







