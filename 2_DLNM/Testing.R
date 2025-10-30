### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city meta-analysis                                         ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta)
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
# source('2_DLNM/HISP_DLNM - Dependencies.R')

#### FUNCTIONS ####

# Function to run model with static formula
run_model <- function(data, locations, outcome){
  
  # 1. Data and dependency prep ----
  # LIST OF DATAFRAMES FOR ALL SITES
  datalist <- lapply(locations, function(region) data[data$location==region,])
  names(datalist) <- locations
  
  # ADDITIONAL CITY-LEVEL META-PREDICTORS
  lat <- unique(data$lat)
  perclat <- round(quantile(lat,c(1,3)*0.25),1)
  
  # TEMPERATURE
  ranges_t2m <- t(sapply(datalist,function(x) range(x$t2m,na.rm=T)))
  cen_t2m <- mean(data$t2m)
  bound_t2m <- c(colMeans(ranges_t2m))
  arg_t2m <- list(argvar=3, arglag=4)
  df_t2m <- prod(unlist(arg_t2m))
  
  # RAINFALL
  ranges_tp <- t(sapply(datalist,function(x) range(x$tp,na.rm=T)))
  cen_tp <- mean(data$tp)
  bound_tp <- c(colMeans(ranges_tp))
  arg_tp <- list(argvar=3, arglag=4)
  df_tp <- prod(unlist(arg_tp))
  
  # # DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
  # knotperc <- c(10,30,70,90)
  # knots <- rowMeans(sapply(datalist,function(x) 
  #   quantile(x$t2m,knotperc/100,na.rm=T)))
  
  # 2. Stage 1: Fit the DLNM ----
  
  # BUILT OBJECTS WHERE outputS WILL BE STORED:
  ymat_t2m <- matrix(NA,m,df_t2m,dimnames=list(locations,paste("spl",seq(df_t2m),sep="")))
  Slist_t2m <- vector("list",m)
  names(Slist_t2m) <- locations
  
  ymat_tp <- matrix(NA,m,df_tp,dimnames=list(locations,paste("spl",seq(df_tp),sep="")))
  Slist_tp <- vector("list",m)
  names(Slist_tp) <- locations
  
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
      
      # Set formula
      form <- as.formula(paste0(outcome, ' ~  season + weekday + closed + ns(time, df=9*9) + cb.t2m + cb.tp'))
      
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
      
    })

  
  # 3. Stage 2: Multivariate Meta-analysis ----
  
  meta_t2m <- get_meta_analysis(ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m)
  meta_tp <- get_meta_analysis(ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_t2m)
  
  # RESET WARNING
  options(warn=0)
  
  return(list(meta_t2m=meta_t2m, meta_tp=meta_tp))
  
}

# Function to get reference meta-analysis 
get_meta_analysis <- function(ymat, Slist, arg, ranges, cen){
  
  # Meta-analysis w/ prediction basis
  mv <- mvmeta(ymat, Slist, method="ml")
  bound <- c(colMeans(ranges))
  mean <- seq(bound[1],bound[2],length=30)
  b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
  
  return(list(ymat=ymat, Slist=Slist, arg=arg, ranges=ranges, mv=mv, bound=bound, b_mean=b_mean, cen=cen))
  
}

# Function to get usable tables and figures
get_outputs <- function(output, city, lags, probs=0.5, outcome){
  
  meta_t2m <- output$meta_t2m
  meta_tp <- output$meta_tp
  
  # Model selection criteria
  tab <- matrix(NA,2,12)
  colnames(tab) <- c("Q","df","p","I-square","AIC","BIC","stat","df","p","stat","df","p")
  rownames(tab) <- c("t2m","tp")
  tab[1,] <- ftab(m=meta_t2m$mv)
  tab[2,] <- ftab(m=meta_tp$mv)
  
  # Exposure response (with city specs) - Overall
  p_er_t2m <- city_exp_res(city=city, lag=lags, meta=meta_t2m, col="#D74B3C", xlab="Temperature (°C)", title=NULL)
  p_er_tp <- city_exp_res(city=city, lag=lags, meta=meta_tp, col="#1D98EB", xlab="Rainfall (mm)", title=NULL)
  
  # Lag response (with city specs) - e.g. 95th percentile
  q_t2m <- quantile(df$t2m, probs=probs, na.rm = TRUE)
  q_tp <- quantile(df$tp, probs=probs, na.rm = TRUE)
  p_lr_t2m <- city_lag_res(city=city, quant_v=q_t2m, meta=meta_t2m, col="#D74B3C", xlab="Temperature (°C)", title=NULL)
  p_lr_tp <- city_lag_res(city=city, quant_v=q_tp, meta=meta_tp, col="#1D98EB", xlab="Rainfall (mm)", title=NULL)
  
  return(list(tab=tab, 
              p_er_t2m=p_er_t2m, p_lr_t2m=p_lr_t2m, 
              p_er_tp=p_er_tp, p_lr_tp=p_lr_tp))
  
}

# 
get_results <- function(output, outcome){
  
  overall <- get_outputs(output, city='ALL', lags=NULL, outcome=outcome)
  overall <- overall[!names(overall) %in% "tab"]
  p_blup_all <- lapply(overall, function(x) x[["p_blup"]])
  p_blup_er <- p_blup_all[grepl("er", names(p_blup_all))]
  p_blup_lr <- p_blup_all[grepl("lr", names(p_blup_all))]
  p_er_overall <- wrap_plots(p_blup_er, ncol=2)+plot_annotation(title=outcome,
                                                                theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  quant_test <- c(0.1, 0.25, 0.75, 0.9)
  quant <- list()
  for (q in quant_test){
    quant[[as.character(q)]] <- get_outputs(output, city='ALL', lags=NULL, probs=q, outcome=outcome)
  }
  p_lr_t2m <- wrap_plots(lapply(quant, function(x) x$p_lr_t2m$p_blup), ncol=4)
  p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_blup), ncol=4)
  p_lr <- wrap_plots(p_lr_t2m, p_lr_tp, nrow=2)+plot_annotation(title=outcome,
                                                                theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  lag_test <- c(0, 5, 10, 14)
  lags <- list()
  for (l in lag_test){
    lags[[as.character(l)]] <- get_outputs(output, city='ALL', lags=l, outcome=outcome)
  }
  p_er_t2m <- wrap_plots(lapply(lags, function(x) x$p_er_t2m$p_blup), ncol=4)
  p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_blup), ncol=4)
  p_er <- wrap_plots(p_er_t2m, p_er_tp, nrow=2)+plot_annotation(title=outcome,
                                                                theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  city <- list()
  for (i in locations){
    city[[i]] <- get_outputs(output, city=i, lags=NULL, outcome=outcome)
  }
  p_er_t2m_city <- wrap_plots(lapply(city, function(x) x$p_er_t2m$p_city), ncol=3)+plot_annotation(title=outcome,
                                                                                                   theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_er_tp_city  <- wrap_plots(lapply(city, function(x) x$p_er_tp$p_city), ncol=3)+plot_annotation(title=outcome,
                                                                                                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_lr_t2m_city  <- wrap_plots(lapply(city, function(x) x$p_lr_t2m$p_city), ncol=3)+plot_annotation(title=outcome,
                                                                                                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  p_lr_tp_city <- wrap_plots(lapply(city, function(x) x$p_lr_tp$p_city), ncol=3)+plot_annotation(title=outcome,
                                                                                                 theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  return(list(overall=overall, p_er_overall=p_er_overall, 
              quant=quant, p_lr=p_lr,
              lags=lags, p_er=p_er,
              city=city, p_er_t2m_city=p_er_t2m_city, p_er_tp_city=p_er_tp_city,
              p_lr_t2m_city=p_lr_t2m_city, p_lr_tp_city=p_lr_tp_city))
}

# Function to get model comparison criteria
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

# Function to get overall & city specific exposure-response relationships
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

# Function to get overall & city specific lag-response relationships
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
        geom_point(data = df_overall, aes(x = lag, y = RR, color = "Pooled"), size = 3.5) +
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
        geom_point(data = df_specific2, aes(x = lag, y = RR, color = "BLUP"), size = 3.5) +
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

# Run model ---

output <- run_model(df, locations, 'total_patients')
main_res <- get_results(output, 'total_patients')
print(main_res$p_er_overall)
print(main_res$p_lr)
print(main_res$p_er)

outcomes <- c('total_patients', 'F', 'M')
out_er5_t2m <- list()
out_er10_t2m <- list()
out_er5_tp <- list()

for (i in outcomes){
  output <- run_model(df, locations, i)
  out_res <- get_results(output, i)
  print(out_res$p_er)
  out_er5_t2m[[i]] <- out_res$lags$`5`$p_er_t2m$p_blup + labs(title=i)
  out_er10_t2m[[i]] <- out_res$lags$`10`$p_er_t2m$p_blup + labs(title=i)
  out_er5_tp[[i]] <- out_res$lags$`14`$p_er_tp$p_blup + labs(title=i)
}

wrap_plots(out_er5_t2m, ncol = 3) 
wrap_plots(out_er10_t2m, ncol = 3) 
wrap_plots(out_er5_tp, ncol = 3) 





# FIG: Meta-regression outputs (e.g. elevation) ----







# TAB: Cochran Q test, I2 and information criteria for model selection ----
tab2 <- matrix(NA,4,12)
colnames(tab2) <- c("Q","df","p","I-square","AIC","BIC","stat","df","p","stat","df","p")
rownames(tab2) <- c("intercept-only","with latitude","rel intercept-only","rel with latitude")

# FUNCTION TO COMPUTE THE STATISTICS FOR EACH MODEL
# THESE COMPUTATIONS ARE LIKELY TO BE REPLACED BY PROPER anova METHODS IN
#   FUTURE RELEASES OF THE PACKAGE
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

mv_t2m <- mvmeta(ymat_t2m, Slist_t2m, method="ml")
tab2[1,] <- ftab(ymat_t2m, Slist_t2m)

tab2[1,] <- ftab(mv_t2m)
tab2[2,] <- ftab(mvlat_t2m,mv_t2m)
# THE TABLE WILL BE COMPLETED LATER


# # GET TABLE
# cp_t2m <- crosspred(b_t2m_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m),model.link="log")[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
# cp_tp <- crosspred(b_tp_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m),model.link="log")[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
# cp <- data.frame(rbind(cp_t2m, cp_tp))
# df_cp <- data.frame(variables = rownames(cp), cp, row.names = NULL)


# COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS
#   1ST STEP: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF META-PREDICTOR
#   2ND STEP: PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS
# (mvlat_t2m <- mvmeta(ymat_t2m~lat,Slist_t2m,method="ml")) # meta-regression
# predlat <- predict(mvlat_t2m,data.frame(lat=perclat),vcov=T)
# cplat25 <- crosspred(b_t2m_mean,coef=predlat[[1]]$fit,vcov=predlat[[1]]$vcov, model.link="log",by=0.1,cen=cen_t2m)
# cplat75 <- crosspred(b_t2m_mean,coef=predlat[[2]]$fit,vcov=predlat[[2]]$vcov, model.link="log",by=0.1,cen=cen_t2m)


# mv_t2m <- mvmeta(ymat_t2m, Slist_t2m, method="ml")
# t2m_mean <- seq(bound_t2m[1],bound_t2m[2],length=30)
# b_t2m_mean <- crossbasis(t2m_mean, lag=c(0,14),argvar=list(fun='ns', df=arg_t2m$argvar), arglag=list(fun='ns', df=arg_t2m$arglag))
# 
# # SET par OPTIONS AND MULTIPANEL
# par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
# layout(1)
# 
# plot(cplat25,"overall",type="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
#      xlab="Temperature (C)",ci.arg=list(density=20,col=grey(0.7)))
# lines(cplat25,"overall",col=1,lty=4,lwd=2)
# lines(cplat75,"overall",col=1,lty=5,lwd=2,ci="area",
#       ci.arg=list(density=20,angle=-45,col=grey(0.7)))
# abline(h=1)
# legend("top",paste(perclat),lty=c(4,5),lwd=2,cex=1,inset=0.1,bty="n",
#        title="Latitude (degree North)")
# 

