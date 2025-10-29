### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city meta-analysis                                         ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta)
library(ggplot2); library(dplyr); library(purrr)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
# source('2_DLNM/HISP_DLNM - Dependencies.R')

city_exp_res <- function(city='ALL', ymat, Slist, arg, ranges, cen, col, xlab=NULL, lag=NULL){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col, title=NULL) {
    
    if(city=='ALL'){
      
      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_line(data = df_specific1, aes(x = pred, y = RR, group = city, color = "City-specific", linetype = "City-specific")) +
        geom_ribbon(data = df_overall, aes(x = pred, ymin = low, ymax = high), fill = col, alpha = 0.2) +
        geom_line(data = df_overall, aes(x = pred, y = RR, color = "Pooled", linetype = "Pooled"), size = 1.2) +
        scale_color_manual(name = NULL, values = c("City-specific" = "grey70", "Pooled" = col)) +
        scale_linetype_manual(name = NULL, values = c("City-specific" = "dashed", "Pooled" = "solid")) +
        scale_x_continuous(breaks = pretty(df_overall$pred, n = 6)) +
        labs(x = xlab, y = "RR", title = title) +
        coord_cartesian(ylim = range(df_specific1$RR), xlim = range(df_overall$pred)) +
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
        data.frame(pred = df_specific1$pred, RR = df_specific1$RR, low = df_specific1$low, high = df_specific1$high, type = "First-stage"),
        data.frame(pred = df_specific2$pred, RR = df_specific2$RR, low = df_specific2$low, high = df_specific2$high, type = "BLUP")
      )
      
      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_ribbon(data = df_specific2, aes(x = pred, ymin = low, ymax = high), fill = col, alpha = 0.2) +
        geom_line(data = df_all, aes(x = pred, y = RR, color = type, linetype = type, size = type)) +
        scale_color_manual(name = NULL, values = c("Pooled" = "black", "First-stage" = "black", "BLUP" = col)) +
        scale_linetype_manual(name = NULL, values = c("Pooled" = "solid", "First-stage" = "dashed", "BLUP" = "solid")) +
        scale_size_manual(name = NULL, values = c("Pooled" = 1.2, "First-stage" = 0.7, "BLUP" = 0.7)) +
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
      
      # df_all <- rbind(
      #   data.frame(pred=df_overall$pred, RR=df_overall$RR, low=df_overall$low, high=df_overall$high, type = "Overall"),
      #   data.frame(pred=df_specific1$pred, RR=df_specific1$RR, low=df_specific1$low, high=df_specific1$high, type = "First-stage"),
      #   data.frame(pred=df_specific2$pred, RR=df_specific2$RR, low=df_specific2$low, high=df_specific2$high, type = "BLUP")
      # )
      # 
      # ggplot() +
      #   geom_ribbon(data = df_specific2, aes(x = pred, ymin = low, ymax = high), fill = col, alpha = 0.2) +
      #   geom_line(data = df_all, aes(x = pred, y = RR, color = type, linetype = type, size = type)) +
      #   geom_hline(yintercept = 1, color = "gray20") +
      #   scale_color_manual(values = c("Overall" = 'black', "First-stage" = 'black', "BLUP" = col)) +
      #   scale_linetype_manual(values = c("Overall" = "solid", "First-stage" = "dashed", "BLUP" = "solid")) +
      #   scale_size_manual(values = c("Overall" = 1.2, "First-stage" = 0.7, "BLUP" = 0.7)) +
      #   scale_x_continuous(breaks = pretty(df_all$pred, n = 6)) +
      #   scale_y_continuous(breaks = pretty(c(min(df_all$low), max(df_all$high)), n = 6)) +
      #   labs(x = xlab, y = "RR", title = title, color = "Estimate type", linetype = "Estimate type", size = "Estimate type") +
      #   coord_cartesian(ylim=c(min(df_all$low), max(df_all$high)), xlim = range(df_overall$pred)) +
      #   theme_classic() +
      #   theme(
      #     axis.line = element_line(color = "gray20"),
      #     axis.text = element_text(color = "gray20"),
      #     axis.ticks = element_line(color = "gray20"),
      #     plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      #     legend.position = "bottom"
      #   )
      
    }
  }
  
  # Meta-analysis w/ prediction basis
  mv <- mvmeta(ymat, Slist, method="ml")
  bound <- c(colMeans(ranges))
  mean <- seq(bound[1],bound[2],length=30)
  b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
  
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
  
  # Get city specific first stage AND BLUP estimates
  if(city=='ALL'){
    
    # Set title based on lags
    if(!is.null(lag)){
      title <- paste('Overall (Lag=', lag, ')')
    } else {
      title <- 'Overall'
    }
    # Plot average vs city estimates
    p_raw  <- plot_estimates(city=city, df_specific1=df_city, df_overall=df_overall, xlab=xlab, col=col, title=paste0(title, ' (first stage estimates)'))
    # Plot average vs BLUP estimates
    p_blup <- plot_estimates(city=city, df_specific1=df_blup, df_overall=df_overall, xlab=xlab, col=col, title=paste0(title, ' (BLUP estimates)'))
    p_city <- NULL

  } else{
    
    # Set title based on lags
    if(!is.null(lag)){
      title <- paste(city, ' (Lag=', lag, ')')
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

city_lag_res <- function(city='ALL', ymat, Slist, arg, ranges, cen, col, xlab=NULL, quant_v=NULL){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col) {
    
    if(city=='ALL'){
      
      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_line(data = df_specific1, aes(x = lag, y = RR, group = city, color = "City-specific", linetype = "City-specific")) +
        geom_line(data = df_overall, aes(x = lag, y = RR, color = "Pooled", linetype = "Pooled"), size = 1.2) +
        geom_errorbar(data = df_overall, aes(x = lag, ymin = lower, ymax = upper, color = "Pooled"), width = 0.2, size = 0.2) +
        geom_point(data = df_overall, aes(x = lag, y = RR, color = "Pooled"), size = 3.5) +
        scale_color_manual(name = NULL, values = c("City-specific" = "grey70", "Pooled" = col)) +
        scale_linetype_manual(name = NULL, values = c("City-specific" = "dashed", "Pooled" = "solid")) +
        scale_x_continuous(breaks = seq(min(df_overall$lag), max(df_overall$lag), by = 2)) +
        labs(x = "Lag (days)", y = "RR", title = paste0("Overall (", names(quant_v), ")")) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "bottom"
        )
      
    } else{

      ggplot() +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_line(data = df_overall, aes(x = lag, y = RR, color = "Pooled", linetype = "Pooled"), size = 1.2) +
        geom_line(data = df_specific1, aes(x = lag, y = RR, group = city, color = "First stage", linetype = "First stage")) +
        geom_line(data = df_specific2, aes(x = lag, y = RR, group = city, color = "BLUP", linetype = "BLUP"), size = 1) +
        geom_errorbar(data = df_specific2, aes(x = lag, ymin = lower, ymax = upper, color = "BLUP"), width = 0.2, size = 0.2) +
        geom_point(data = df_specific2, aes(x = lag, y = RR, color = "BLUP"), size = 3.5) +
        scale_color_manual(name = NULL, values = c("Pooled" = "black", "First stage" = col, "BLUP" = col)) +
        scale_linetype_manual(name = NULL, values = c("Pooled" = "solid", "First stage" = "dashed", "BLUP" = "solid")) +
        scale_x_continuous(breaks = seq(min(df_overall$lag), max(df_overall$lag), by = 2)) +
        labs(x = "Lag (days)", y = "RR", title = paste0(city, " (", names(quant_v), ")")) +
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
  
  # Meta-analysis w/ prediction basis
  mv <- mvmeta(ymat, Slist, method="ml")
  bound <- c(colMeans(ranges))
  mean <- seq(bound[1],bound[2],length=30)
  b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
  
  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link="log", cen=cen, at=quant_v)
  df_overall <- data.frame(lag=seq(cp$lag[1], cp$lag[2]), RR=t(cp$matRRfit), lower=t(cp$matRRlow), upper=t(cp$matRRhigh), row.names = NULL)
  names(df_overall) <- c('lag', 'RR', 'lower','upper')

  # Get city specific estimates (for specific lag, if applicable)
  df_city <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=ymat[i,], vcov=Slist[[i]], model.link="log", 
                         from=ranges[i, 1], to=ranges[i, 2], cen=cen, at=quant_v)
    data.frame(city=rownames(ymat)[i], lag=seq(cp_city$lag[1], cp_city$lag[2]), RR=t(cp_city$matRRfit), lower=t(cp_city$matRRlow), upper=t(cp_city$matRRhigh), row.names = NULL)
  })
  names(df_city) <- c('city', 'lag', 'RR', 'lower','upper')
  
  # Get BLUP estimates (for specific lag, if applicable)
  blup <- blup(mv, vcov=TRUE)
  df_blup <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef=blup[[i]]$blup, vcov=blup[[i]]$vcov,
                         model.link="log", from=ranges[i, 1], to=ranges[i, 2], cen=cen, at=quant_v)
    data.frame(city=names(blup)[i], lag=seq(cp_city$lag[1], cp_city$lag[2]), RR=t(cp_city$matRRfit), lower=t(cp_city$matRRlow), upper=t(cp_city$matRRhigh), row.names = NULL)
  })
  names(df_blup) <- c('city', 'lag', 'RR', 'lower','upper')
  
  # Get city specific first stage AND BLUP estimates
  if(city=='ALL'){
    
    p_raw  <- plot_estimates(city=city, df_specific1=df_city, df_overall=df_overall, xlab=xlab, col=col)
    p_blup  <- plot_estimates(city=city, df_specific1=df_blup, df_overall=df_overall, xlab=xlab, col=col)
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
data <- df

# 1. Data and dependency prep #####

# REGIONS
regions <- as.character(unique(data$location))

# LIST OF DATAFRAMES FOR ALL SITES
datalist <- lapply(regions, function(region) data[data$location==region,])
names(datalist) <- regions
m <- length(datalist)

# ADDITIONAL CITY-LEVEL META-PREDICTORS
lat <- unique(data$lat)
perclat <- round(quantile(lat,c(1,3)*0.25),1)

# # MOVING AVERAGE OF t2m_mean OVER LAG 0-6
# for(i in seq(datalist)) datalist[[i]]$t2m <- 
#   filter(datalist[[i]]$t2m,rep(1,6)/6,side=1)


# For all predictors...
## - Calculate the site specific range
## - Overall bounds, centering point and degrees of freedom
## - Compute the 25-75th percentiles

# MEAN TEMPERATURE
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

# COMPUTE 25TH-75TH PERCENTILES OF META-VARIABLES

# DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
# (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)

# degree <- 3
# type <- "bs"
# df <- 12

# DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
knotperc <- c(10,30,70,90)
knots <- rowMeans(sapply(datalist,function(x) 
  quantile(x$t2m,knotperc/100,na.rm=T)))

# 2. Stage 1: Fit the DLNM #####

# BUILT OBJECTS WHERE RESULTS WILL BE STORED:
#   ymat_t2m IS THE MATRIX FOR THE T2M PARAMETERS
#   Slist_t2m IS THE LIST WITH (CO)VARIANCE MATRICES FOR T2M

ymat_t2m <- matrix(NA,m,df_t2m,dimnames=list(regions,paste("spl",seq(df_t2m),sep="")))
Slist_t2m <- vector("list",m)
names(Slist_t2m) <- regions

ymat_tp <- matrix(NA,m,df_tp,dimnames=list(regions,paste("spl",seq(df_tp),sep="")))
Slist_tp <- vector("list",m)
names(Slist_tp) <- regions

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
    
    # RUN THE MODEL
    model <- glm(total_patients ~ cb.t2m + cb.tp + weekday + ns(time,7*14),family=quasipoisson(),data)
    
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


# RESET WARNING
options(warn=0)

# 3. Stage 2: Multivariate Meta-analysis results ####

# Exposure response (with city specs) - Overall
(t <- city_exp_res(city='ALL', ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- city_exp_res(city='Chamba', ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- city_exp_res(city='ALL', ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))
(t <- city_exp_res(city='Chamba', ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))

# Exposure response (with city specs) - At lag 5
(t <- city_exp_res(city='ALL', lag=5, ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- city_exp_res(city='Chamba', lag=5, ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- city_exp_res(city='ALL', lag=5, ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))
(t <- city_exp_res(city='Chamba', lag=5, ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))

# Specific city exposure response - Overall
for (i in regions){
  t <- city_exp_res(city=i, ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)")
  print(t$p_city)
}

# Lag response (with city specs) - 95th percentile
q95_t2m <- quantile(df$t2m, probs = 0.95, na.rm = TRUE)
q95_tp <- quantile(df$tp, probs = 0.95, na.rm = TRUE)
(t <- city_lag_res(city='ALL', quant_v=q95_t2m, ymat=ymat_t2m, Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- city_lag_res(city='Chamba', quant_v=q95_t2m, ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- city_lag_res(city='ALL', quant_v=q95_tp, ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))
(t <- city_lag_res(city='Chamba', quant_v=q95_tp, ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))


# FIG: Meta-regression outputs (e.g. elevation) ----

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


mv_t2m <- mvmeta(ymat_t2m, Slist_t2m, method="ml")
t2m_mean <- seq(bound_t2m[1],bound_t2m[2],length=30)
b_t2m_mean <- crossbasis(t2m_mean, lag=c(0,14),argvar=list(fun='ns', df=arg_t2m$argvar), arglag=list(fun='ns', df=arg_t2m$arglag))

# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(1)

plot(cplat25,"overall",type="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
     xlab="Temperature (C)",ci.arg=list(density=20,col=grey(0.7)))
lines(cplat25,"overall",col=1,lty=4,lwd=2)
lines(cplat75,"overall",col=1,lty=5,lwd=2,ci="area",
      ci.arg=list(density=20,angle=-45,col=grey(0.7)))
abline(h=1)
legend("top",paste(perclat),lty=c(4,5),lwd=2,cex=1,inset=0.1,bty="n",
       title="Latitude (degree North)")


# TAB: Cochran Q test, I2 and information criteria for model selection ----
tab2 <- matrix(NA,4,12)
colnames(tab2) <- c("Q","df","p","I-square","AIC","BIC","stat","df","p",
                    "stat","df","p")
rownames(tab2) <- c("intercept-only","with latitude","rel intercept-only",
                    "rel with latitude")

# FUNCTION TO COMPUTE THE STATISTICS FOR EACH MODEL
# THESE COMPUTATIONS ARE LIKELY TO BE REPLACED BY PROPER anova METHODS IN
#   FUTURE RELEASES OF THE PACKAGE
ftab <- function(m,mref=NULL) {
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
  # RESULTS
  if(!is.null(mref)) {
    return(c(het,lr,wald))
  } else return(c(het,rep(NA,6))) 
}
tab2[1,] <- ftab(mv_t2m)
tab2[2,] <- ftab(mvlat_t2m,mv_t2m)
# THE TABLE WILL BE COMPLETED LATER