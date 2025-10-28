
# LOAD THE PACKAGE
library(dlnm) ; library(splines) ; library(xtable) ; library(mvmeta)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
# source('2_DLNM/HISP_DLNM - Dependencies.R')

#### MAIN CODE ####

load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
data <- df

# 1. Data and dependency prep #####

# REGIONS
regions <- as.character(unique(data$location))

# LIST OF DATAFRAMES FOR 10 REGIONS
datalist <- lapply(regions, function(region) data[data$location==region,])
names(datalist) <- regions

# CITY-LEVEL META-PREDICTORS
lat <- unique(data$lat)
perclat <- round(quantile(lat,c(1,3)*0.25),1)

# ADDITIONAL INFO
m <- length(datalist)

# # MOVING AVERAGE OF TMEAN OVER LAG 0-6
# for(i in seq(datalist)) datalist[[i]]$t2m <- 
#   filter(datalist[[i]]$t2m,rep(1,6)/6,side=1)

# TEMPERATURE RANGES (FOR LAG 0-5)
ranges <- t(sapply(datalist,function(x) range(x$t2m,na.rm=T)))

# COMPUTE 25TH-75TH PERCENTILES OF META-VARIABLES

# DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
# (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
cen <- mean(data$t2m)
bound <- c(colMeans(ranges))
degree <- 2
type <- "bs"
df <- 6

# DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
knotperc <- c(10,30,70,90)
knots <- rowMeans(sapply(datalist,function(x) 
  quantile(x$t2m,knotperc/100,na.rm=T)))

# 2. Stage 1: Fit the DLNM #####

# BUILT OBJECTS WHERE RESULTS WILL BE STORED:
#   ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   Slist IS THE LIST WITH (CO)VARIANCE MATRICES

ymat <- matrix(NA,m,df,dimnames=list(regions,paste("spl",seq(df),sep="")))
Slist <- vector("list",m)
names(Slist) <- regions

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,"")
    
    # LOAD
    data <- datalist[[i]]
    data <- data[complete.cases(data$t2m, data$total_patients, data$weekday, data$time), ]
    
    # CREATE THE SPLINE
    # NB: KNOTS AND BOUNDARIES FIXED AT SAME VALUES
    bt2m <- onebasis(data$t2m,fun=type,degree=degree,knots=knots,Bound=bound)
    # bt2m <- crossbasis(data$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=data$year)
    
    # DEFINE THE KNOTS AT STUDY-SPECIFIC PERCENTILES
    # tempknots <- quantile(data$t2m,knotperc/100,na.rm=TRUE)
    # cenrel <- quantile(data$t2m,0.75,na.rm=TRUE)
    # 
    # CREATE THE QUADRATIC SPLINE: STUDY-SPECIFIC KNOTS
    # bt2m <- onebasis(data$t2m,fun=type,degree=degree,knots=tempknots)
    
    # RUN THE MODEL
    model <- glm(total_patients ~ bt2m + weekday + ns(time,7*14),family=quasipoisson(),data)

    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    predt2m <- crosspred(bt2m,model,cen=cen)
    ymat[i,] <- predt2m$coef
    Slist[[i]] <- predt2m$vcov
  })


# RESET WARNING
options(warn=0)

# 3. Stage 2: Multivariate Meta-analysis/regression ####

# MULTIVARIATE META-ANALYSIS
(mv <- mvmeta(ymat,Slist,method="ml"))
summary(mv)
# MULTIVARIATE META-REGRESSION MODELS
(mvlat <- mvmeta(ymat~lat,Slist,method="ml"))
summary(mvlat)

# BASIS USED TO PREDICT TEMPERATURE, EQUAL TO THAT USED FOR ESTIMATION
#   NOTE: INTERNAL AND BOUNDARY KNOTS PLACED AT SAME VALUES AS IN ESTIMATION
tmean <- seq(bound[1],bound[2],length=30)
btmean <- onebasis(tmean,fun=type,degree=degree,knots=knots,Bound=bound)


# USE OF crosspred TO PREDICT THE EFFECTS FOR THE CHOSEN VALUES

# PREDICTION FROM SIMPLE META-ANALYSES WITH NO PREDICTORS
# CENTERED TO SPECIFIC VALUE
cp <- crosspred(btmean,coef=coef(mv),vcov=vcov(mv),model.link="log",
                by=0.1,cen=cen)

# COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS
#   1ST STEP: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF META-PREDICTOR
#   2ND STEP: PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS

predlat <- predict(mvlat,data.frame(lat=perclat),vcov=T)
cplat25 <- crosspred(btmean,coef=predlat[[1]]$fit,vcov=predlat[[1]]$vcov,
                     model.link="log",by=0.1,cen=cen)
cplat75 <- crosspred(btmean,coef=predlat[[2]]$fit,vcov=predlat[[2]]$vcov,
                     model.link="log",by=0.1,cen=cen)

# 4. Results ####

# FIG:Temperature Ranges ----
par(mar=c(5,6,4,1)+0.1)
layout(1)
plot(seq(-10,30,length=m),seq(m),type="n",yaxt="n",
     ylab="",xlab="Temperature range",main="Temperature ranges (C) - Lag 0-5")
axis(2,at=seq(m),labels=regions,las=1,cex.axis=0.7)
arrows(ranges[,1],seq(m),ranges[,2],seq(m),angle=90,length=0.05,code=3)
title(ylab="Cities",mgp=c(5,1,0))
abline(v=knots)
abline(v=bound,lty=2)
mtext("Vertical lines identify internal (continuous) and boundary (dashed) knots",
      cex=0.8)

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
tab2[1,] <- ftab(mv)
tab2[2,] <- ftab(mvlat,mv)
# THE TABLE WILL BE COMPLETED LATER

# FIG: Pooled relationship ----

# SET par OPTIONS
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(1)

plot(cp,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
     xlab="Temperature (C)")
points(round(knots,1),cp$allRRfit[as.character(round(knots,1))],pch=19,cex=0.6)

# FIG: BLUP ----

# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(matrix(1:2,1,2))

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

# PLOT OF AVERAGE AND CITY-SPECIFIC ESTIMATES
plot(cp,"overall",type="n",ci="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
     xlim=range(ranges),xlab="Temperature (C)")
for(i in seq(m)) {
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  lines(crosspred(btmean,coef=ymat[i,],vcov=Slist[[i]],model.link="log",
                  from=ranges[i,1],to=ranges[i,2],cen=cen),col=grey(0.8),lty=5)
}
lines(cp,"overall",col=1,lwd=2)
abline(h=1)
mtext("Study-specific",cex=1.3)

# PLOT OF AVERAGE AND BLUP ESTIMATES
blup <- blup(mv,vcov=TRUE)
plot(cp,"overall",type="n",ci="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
     xlim=range(ranges),xlab="Temperature (C)")
for(i in seq(m)) {
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  lines(crosspred(btmean,coef=blup[[i]]$blup,vcov=blup[[i]]$vcov,
                  model.link="log",from=ranges[i,1],to=ranges[i,2],cen=cen),col=grey(0.8),
        lty=5)
}
lines(cp,"overall",col=1,lwd=2)
abline(h=1)
mtext("BLUP",cex=1.3)
options(warn=0)

# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(matrix(1:2,1,2))

plot(cp,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),ci="n",
     xlab="Temperature (C)")
lines(crosspred(btmean,coef=ymat["Chamba",],vcov=Slist[["Chamba"]],
                model.link="log",by=0.1,cen=cen),lty=4)
lines(crosspred(btmean,coef=blup[["Chamba"]]$blup,vcov=blup[["Chamba"]]$vcov,
                model.link="log",by=0.1,cen=cen),lty=5)
mtext("Chamba",cex=1.3)
legend("top",c("Population-average","First-stage","BLUP"),lty=c(1,4,5),
       lwd=c(2,1,1),cex=0.8,inset=0.1,bty="n")

plot(cp,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),ci="n",
     xlab="Temperature (C)")
lines(crosspred(btmean,coef=ymat["Kullu",],vcov=Slist[["Kullu"]],
                model.link="log",by=0.1,cen=cen),lty=4)
lines(crosspred(btmean,coef=blup[["Kullu"]]$blup,vcov=blup[["Kullu"]]$vcov,
                model.link="log",by=0.1,cen=cen),lty=5)
mtext("Kullu",cex=1.3)
legend("top",c("Population-average","First-stage","BLUP"),lty=c(1,4,5),
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



