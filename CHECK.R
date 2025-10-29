
# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta)
library(ggplot2); library(dplyr); library(purrr)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
# source('2_DLNM/HISP_DLNM - Dependencies.R')

get_visuals <- function(cp, model, xlab, title, col, lag_rep=FALSE, lag=NULL){
  
  
  if(lag_rep==FALSE){
    if(!is.null(lag)){
      cum <- data.frame(x = cp$predvar, RR = cp$matRRfit[, lag + 1], lower = cp$matRRlow[, lag + 1], upper = cp$matRRhigh[, lag + 1], row.names = NULL)
    } else {
      cum <- data.frame(x = cp$predvar, RR = cp$allRRfit, lower = cp$allRRlow, upper = cp$allRRhigh, row.names = NULL)
    }
    p <- ggplot(cum, aes(x = x, y = RR)) +
      geom_line(color = col) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = col) +
      geom_hline(yintercept = 1, color = "gray20") +
      scale_x_continuous(breaks = pretty(cum$x, n = 6)) +
      # scale_y_continuous(limits = c(-0.5, 2.5)) +
      labs(x = xlab, y = "RR", title = title) +
      theme_classic() +
      theme(
        axis.line = element_line(color = "gray20"),      # pale axis lines
        axis.text = element_text(color = "gray20"),      # pale axis labels
        axis.ticks = element_line(color = "gray20"),     # pale ticks
        plot.title = element_text(face = "bold", hjust = 0.5, size=12)  # center title
      )  
  } else {
    
    slice <- data.frame(lag=seq(cp$lag[1], cp$lag[2]), RR=t(cp$matRRfit), lower=t(cp$matRRlow), upper=t(cp$matRRhigh), row.names = NULL)
    names(slice) <- c('lag', 'RR', 'lower','upper')
    p <- ggplot(slice, aes(x = lag, y = RR)) +
      geom_hline(yintercept = 1, color = "gray20") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,  size = 0.2) +
      geom_point(color = col,  size = 2.5) +
      scale_x_continuous(breaks = seq(min(slice$lag), max(slice$lag), by = 2)) +  # more breaks
      labs(x = "Lag (days)", y = "RR", title = title) +
      theme_classic() +
      theme(
        axis.line = element_line(color = "gray20"),      # pale axis lines
        axis.text = element_text(color = "gray20"),      # pale axis labels
        axis.ticks = element_line(color = "gray20"),     # pale ticks
        plot.title = element_text(face = "bold", hjust = 0.5, size=12)  # center title
      )
    
  }
  
  return(p)
}

get_city_visuals <- function(city='ALL', ymat, Slist, arg, ranges, cen, col, xlab){
  
  # FUNCTION: Formatting for multi-city plots
  plot_estimates <- function(city=city, df_specific1, df_specific2=NULL, df_overall, xlab, col, title=NULL) {
    
    if(city=='ALL'){

      ggplot() +
        geom_line(data = df_specific1, aes(x = pred, y = RR, group = city),
                  color = "grey70", linetype = "dashed") +
        geom_ribbon(data = df_overall, aes(x = pred, ymin = low, ymax = high), fill = col, alpha = 0.2) +
        geom_line(data = df_overall, aes(x = pred, y = RR), color = col, size = 1.2) +
        geom_hline(yintercept = 1, linetype = "solid") +
        # Dummy lines for legend only
        geom_line(aes(x = 0, y = 0, color = "City-specific", linetype = "City-specific"), size = 0.7) +
        geom_line(aes(x = 0, y = 0, color = "Overall", linetype = "Overall"), size = 1.2) +
        scale_color_manual(name = "Estimate type", values = c("City-specific" = "grey70", "Overall" = col)) +
        scale_linetype_manual(name = "Estimate type", values = c("City-specific" = "dashed", "Overall" = "solid")) +
        scale_x_continuous(breaks = pretty(df_overall$pred, n = 6)) +
        labs(x = xlab, y = "RR", title = title) +
        coord_cartesian(ylim = c(0.5, 1.5), xlim = range(df_overall$pred)) +
        theme_classic() +
        theme(axis.line = element_line(color = "gray20"),
              axis.text = element_text(color = "gray20"),
              axis.ticks = element_line(color = "gray20"),
              plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
      
      
    } else{
      df_all <- rbind(
        data.frame(pred=df_overall$pred, RR=df_overall$RR, low=df_overall$low, high=df_overall$high, type = "Overall"),
        data.frame(pred=df_specific1$pred, RR=df_specific1$RR, low=df_specific1$low, high=df_specific1$high, type = "First-stage"),
        data.frame(pred=df_specific2$pred, RR=df_specific2$RR, low=df_specific2$low, high=df_specific2$high, type = "BLUP")
      )
      # Create a numeric mapping for size
      df_all$size <- ifelse(df_all$type == "Overall", 1.2,
                            ifelse(df_all$type == "First-stage", 0.7, 0.7))  # BLUP thinner
      
      ggplot() +
        geom_ribbon(data = df_specific2, aes(x = pred, ymin = low, ymax = high), fill = 'grey50', alpha = 0.2) +
        geom_line(data = df_all, aes(x = pred, y = RR, color = type, linetype = type, size = type)) +
        geom_hline(yintercept = 1, color = "gray20") +
        scale_color_manual(values = c("Overall" = col, "First-stage" = col, "BLUP" = 'black')) +
        scale_linetype_manual(values = c("Overall" = "solid", "First-stage" = "dashed", "BLUP" = "solid")) +
        scale_size_manual(values = c("Overall" = 1.2, "First-stage" = 1, "BLUP" = 0.7)) +
        scale_x_continuous(breaks = pretty(df_overall$pred, n = 6)) +
        labs(x = xlab, y = "RR", title = city, color = "Estimate type", linetype = "Estimate type", size = "Estimate type") +
        coord_cartesian(ylim = c(0.5, 1.5), xlim = range(df_overall$pred)) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
        )
      
    }
  }
  
  # Meta-analysis w/ prediction basis
  mv <- mvmeta(ymat, Slist, method="ml")
  bound <- c(colMeans(ranges))
  mean <- seq(bound[1],bound[2],length=30)
  b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
  
  # Get overall estimates
  cp <- crosspred(b_mean, coef=coef(mv),vcov=vcov(mv), model.link="log", cen=cen)
  df_overall <- data.frame(pred = cp$predvar, RR = cp$allRRfit,
                           low = cp$allRRlow, high = cp$allRRhigh)
  
  # Get city specific estimates
  df_city <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef = ymat[i,], vcov = Slist[[i]], model.link = "log", 
                         from = ranges[i, 1], to = ranges[i, 2], cen = cen)
    data.frame(city = rownames(ymat)[i], pred = cp_city$predvar, RR = cp_city$allRRfit,
               low = cp_city$allRRlow, high = cp_city$allRRhigh)
  })
  
  # Get BLUP estimates
  blup <- blup(mv, vcov=TRUE)
  df_blup <- map_dfr(seq_len(m), function(i) {
    cp_city <- crosspred(b_mean, coef = blup[[i]]$blup, vcov = blup[[i]]$vcov,
                         model.link = "log", from = ranges[i, 1], to = ranges[i, 2], cen = cen)
    data.frame(city = names(blup)[i], pred = cp_city$predvar, RR = cp_city$allRRfit,
               low = cp_city$allRRlow, high = cp_city$allRRhigh)
  })
  
  
  if(city=='ALL'){
    # Plot average vs city estimates
    p_raw  <- plot_estimates(city=city, df_specific1=df_city, df_overall=df_overall, xlab=xlab, col=col, title="Study-specific & Overall Estimates")
    # Plot average vs BLUP estimates
    p_blup <- plot_estimates(city=city, df_specific1=df_blup, df_overall=df_overall, xlab=xlab, col=col, title="BLUP & Overall Estimates")
    p_city <- NULL
    
  } else{
    
    df_city <- df_city[df_city$city==city,]
    df_blup <- df_blup[df_blup$city==city,]
    
    # Plot average vs city specific estimates
    p_raw <- NULL
    p_blup <- NULL
    p_city  <- plot_estimates(city=city, df_specific1=df_city, df_specific2=df_blup, df_overall=df_overall, xlab=xlab, col=col)

  }

  return(list(df_overall=df_overall, df_city=df_city, p_raw=p_raw, p_blup=p_blup, p_city=p_city))
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

# 3. Stage 2: Multivariate Meta-analysis/regression ####

# COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS
#   1ST STEP: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF META-PREDICTOR
#   2ND STEP: PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS
# (mvlat_t2m <- mvmeta(ymat_t2m~lat,Slist_t2m,method="ml")) # meta-regression
# predlat <- predict(mvlat_t2m,data.frame(lat=perclat),vcov=T)
# cplat25 <- crosspred(b_t2m_mean,coef=predlat[[1]]$fit,vcov=predlat[[1]]$vcov, model.link="log",by=0.1,cen=cen_t2m)
# cplat75 <- crosspred(b_t2m_mean,coef=predlat[[2]]$fit,vcov=predlat[[2]]$vcov, model.link="log",by=0.1,cen=cen_t2m)

# 4. Results ####



# CITY SPECIFIC
(t <- get_city_visuals(city='ALL', ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- get_city_visuals(city='Chamba', ymat=ymat_t2m,  Slist=Slist_t2m, arg=arg_t2m, ranges=ranges_t2m, cen=cen_t2m, col="#D74B3C", xlab="Temperature (°C)"))
(t <- get_city_visuals(city='ALL', ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)"))

for (i in regions){
  t <- get_city_visuals(city=i, ymat=ymat_tp,  Slist=Slist_tp, arg=arg_tp, ranges=ranges_tp, cen=cen_tp, col="#1D98EB", xlab="Rainfall (mm)")
  print(t$p_city)
}




# Set overall quantile
q95_t2m <- quantile(df$t2m, probs = 0.95, na.rm = TRUE)
q05_t2m <- quantile(df$t2m, probs = 0.05, na.rm = TRUE)
q95_tp <- quantile(df$tp, probs = 0.95, na.rm = TRUE)

# GET TABLE
cp_t2m <- crosspred(b_t2m_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m),model.link="log")[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
cp_tp <- crosspred(b_tp_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m),model.link="log")[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
cp <- data.frame(rbind(cp_t2m, cp_tp))
df_cp <- data.frame(variables = rownames(cp), cp, row.names = NULL)

# GET FIGURES
# T2M
cp <- crosspred(b_t2m_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m), model.link="log", cen=cen_t2m, at=q95_t2m)
(p_lr_t2m95 <- get_visuals(cp, mv_t2m, lag_rep=TRUE, col="#D74B3C", title='Lag-response of t2m (at 95th percentile)'))
cp <- crosspred(b_t2m_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m), model.link="log", cen=cen_t2m, at=q05_t2m)
(p_lr_t2m05 <- get_visuals(cp, mv_t2m, lag_rep=TRUE, col="#D74B3C", title='Lag-response of t2m (at 5th percentile)'))
cp <- crosspred(b_t2m_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m), model.link="log", cen=cen_t2m) # at=seq(-10, 40, by=1)
(p_er_t2m <- get_visuals(cp, mv_t2m, lag_rep=FALSE, lag=NULL, col="#D74B3C", xlab="Temperature (C)", title="RR of t2m exposure-response (lag=5)"))
# TP
cp <- crosspred(b_tp_mean, coef=coef(mv_tp),vcov=vcov(mv_tp), model.link="log", at=q95_tp)
(p_lr_tp <- get_visuals(cp, mv_tp, lag_rep=TRUE, col="#1D98EB", title='Lag-response of tp (at 95th percentile)'))
cp <- crosspred(b_tp_mean, coef=coef(mv_tp),vcov=vcov(mv_tp), model.link="log", cen=cen_tp) #at=seq(0, 100, by=1)
(p_er_tp <- get_visuals(cp, mv_tp, lag_rep=FALSE, lag=NULL, col="#1D98EB", xlab="Rainfall (mm)", title="RR of tp exposure-response (lag=0)"))



# FIG: Pooled relationship ----




# ----------------------------------


# TEMPERATURE
# Multivariate meta-analysis/regression
# Set prediction basis
t2m_mean <- seq(bound_t2m[1],bound_t2m[2],length=30)
b_t2m_mean <- crossbasis(t2m_mean, lag=c(0,14),argvar=list(fun='ns', df=arg_t2m$argvar), arglag=list(fun='ns', df=arg_t2m$arglag))
cp_t2m <- crosspred(b_t2m_mean, coef=coef(mv_t2m),vcov=vcov(mv_t2m), model.link="log", at=seq(-10, 40, by=1), cen=cen_t2m)

plot(cp_t2m,"overall",col=1,lwd=2,ylab="RR",ylim=c(.5,2.5),xlim=c(-10,40), xlab="Temperature (C)")

# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(matrix(1:2,1,2))

plot(cp_t2m,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),ci="n",
     xlab="Temperature (C)")
lines(crosspred(b_t2m_mean,coef=ymat_t2m["Chamba",],vcov=Slist_t2m[["Chamba"]],
                model.link="log",by=0.1,cen=cen_t2m),lty=4)
lines(crosspred(b_t2m_mean,coef=blup[["Chamba"]]$blup,vcov=blup[["Chamba"]]$vcov,
                model.link="log",by=0.1,cen=cen_t2m),lty=5)
mtext("Chamba",cex=1.3)
legend("top",c("Overall","First-stage","BLUP"),lty=c(1,4,5),
       lwd=c(2,1,1),cex=0.8,inset=0.1,bty="n")

plot(cp_t2m,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),ci="n",
     xlab="Temperature (C)")
lines(crosspred(b_t2m_mean,coef=ymat_t2m["Kullu",],vcov=Slist_t2m[["Kullu"]],
                model.link="log",by=0.1,cen=cen_t2m),lty=4)
lines(crosspred(b_t2m_mean,coef=blup[["Kullu"]]$blup,vcov=blup[["Kullu"]]$vcov,
                model.link="log",by=0.1,cen=cen_t2m),lty=5)
mtext("Kullu",cex=1.3)
legend("top",c("Overall","First-stage","BLUP"),lty=c(1,4,5),
       lwd=c(2,1,1),cex=0.8,inset=0.1,bty="n")

# FIG: Meta-regression outputs (e.g. latitude) ----

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