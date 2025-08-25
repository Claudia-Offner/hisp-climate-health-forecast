### ------------------------------------------------------------------------ ### 
###  DNLM Multivariate Model - Results                                       ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/HISP_TSGC/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 2 - Causal Inference/')
setwd(path2)

####  PACKAGES & FUNCTIONS  ####
library(dlnm)
library(splines)
library(dplyr)
library(lubridate)

#### DEPENDENCIES ####
weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')

# ------------------------------------------------------------------------------ 
# ---------------------------- DNLM Development --------------------------------
#  Multivariate model of mx2t, mn2t, tp and fg10 on hospital admissions
# ------------------------------------------------------------------------------

# Load data
df <- read.csv('203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')
# df <- read.csv('203_comb/ind_10km_HISP-HP_2015-2017_hosp-bigthaw_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')

# Set date type and factors 
df$date <- as.Date(df$date, format='%Y-%m-%d')  # Adjust format as needed
df <- df %>%
  mutate(weekday=factor(weekday, levels=weekday_levels),
         month=factor(month, levels=month_levels), 
         season=factor(season, levels=season_levels))

# Select variables of interest
df <- df %>%
  mutate(doy = yday(date),
         time = row_number()) %>%
  select(
    location, time, date, year, season, month, week, weekday, doy, closed, holiday, total_patients, 
    F, M, under_5, age5_14, age15_24, age25_64, over_65,
    general, emergency, ophthalmology, pediatrics, orthopedics, ent, gynecology, dermatology, surgery, dentistry, psychiatry, oncology, 
    # total_patients_c, F_c, M_c, under_5_c, age5_14_c, age15_24_c, age25_64_c, over_65_c,
    # general_c, emergency_c, ophthalmology_c, pediatrics_c, orthopedics_c, ent_c, gynecology_c, dermatology_c, surgery_c, dentistry_c, psychiatry_c, oncology_c, 
    t2m, mn2t, mx2t, tp, fg10,
    mx2t_IMD, mx2t_IMD_b, mn2t_IMD, mn2t_IMD_b, tp_IMD, tp_IMD_b, fg10_IMD, fg10_IMD_b, 
    t2m_p10, t2m_p90, mn2t_p10, mn2t_p90, mx2t_p10, mx2t_p90, tp_p10, tp_p90, fg10_p10, fg10_p90
  )

df$tp <- exp(df$tp)

####################################################
install.packages('TRR')
library(dlnm)
library(glmmTMB)
library(TRR)
# Create crossbasis objects
cb.t2m <- crossbasis(df$t2m, lag=14, argvar=list(fun="ns", df=4), arglag=list(fun="ns", df=5))
# cb.tp <- EMA(df$tp, n = 14)
# cb.fg10 <- EMA(df$fg10, n = 14)


# cb.tp <- crossbasis(log(df$tp + 0.1), lag=14, argvar=list(fun="lin"), arglag=list(fun="ns", df=5))
# cb.fg10 <- crossbasis(df$fg10, lag=14, argvar=list(fun="lin"), arglag=list(fun="ns", df=5))

model <- glm(total_patients ~ cb.t2m + cb.tp + cb.fg10 + ns(time, 7*14) + weekday, 
              family=quasipoisson(), df, na.action=na.omit)

cp.t2m <- crosspred(cb.t2m, model, cen=median(df$t2m, na.rm=TRUE), cumul=TRUE)
cp.tp <- crosspred(cb.tp, model, cen=median(log(df$tp + 0.1), na.rm=TRUE), cumul=TRUE)

plot(cp.t2m, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
     main="Association with a 10-unit increase in temperature")
plot(cp.t2m, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
     main="Cumulative association with a 10-unit increase in temperature")

plot(cp.tp, "slices", var=2, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
     main="Association with a 10-unit increase in rainfall")
plot(cp.tp, "slices", var=2, col=2, cumul=TRUE, ylab="Cumulative RR",
     main="Cumulative association with a 10-unit increase in rainfall")


plot(cp.t2m, "slices", var=-5, ci="n", col=1, ylim=c(0.95,1.25), lwd=1.5,
     main="Lag-response curves for different temperatures, ref. 10C")
for(i in 1:3) lines(cp.t2m, "slices", var=c(0, 15, 20)[i], col=i+1, lwd=1.5)
legend("topright",paste("Temperature =",c(-5, 0, 15, 20)), col=1:4, lwd=1.5)

# Figure 4b: exposure-response relationships specific to lag 0 and 5 (left column), 
# and lag-response relationships specific to temperatures -20◦C and 33◦C (right column)
plot(cp.t2m, "slices", var=c(-5,20), lag=c(0,10), col=4,
     ci.arg=list(density=40,col=grey(0.7)))



#### LINEAR: Simple DLM (t2m VS total_patients - control tp and fg10) ####
# # Create cross basis functions for each predictor (exposures pm10; control temp)
# cb1.temp <- crossbasis(df$t2m, lag=7, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
# cb1.tp <- crossbasis(df$tp, lag=7, argvar=list(df=5), arglag=list(fun="strata",breaks=1))
# cb1.fg10 <- crossbasis(df$fg10, lag=7, argvar=list(df=5), arglag=list(fun="strata",breaks=1))
# 
# # Fit model 
# model1 <- glm(total_patients ~ cb1.temp + cb1.tp + cb1.fg10 + ns(time, 7*14) + weekday, 
#               family=quasipoisson(), df, na.action=na.omit)
# # over dispersed Poisson distribution, including a smooth function of time with 7 df/year
# # (in order to correct for seasonality and long time trend) and day of the week as factor
# 
# # Given this model, estimate the association of specific levels of PM10 on mortality
# pred1.pm <- crosspred(cb1.temp, model1, at=0:30, bylag=0.2, cumul=TRUE)
# # 'at=0:20' states that the prediction must be computed for each integer value from 0 to 20 µgr/m3. 
# # By setting bylag=0.2, the prediction is computed along the lag space with an increment of 0.2. 
# 
# # Lag−response curves
# plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
#      main="Association with a 10-unit increase in temperature")
# plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
#      main="Cumulative association with a 10-unit increase in temperature")
# 
# # Get the overall cumulative effect of a 10-unit increase in PM10 over 15 days of lag & CI
# pred1.pm$allRRfit["10"]
# cbind(pred1.pm$allRRlow, pred1.pm$allRRhigh)["10",]

#### NON-LINEAR: Bi-dimensional (t2m vs total_patients) ####

# cross-basis for temperature is specied through a natural and non-natural splines
# cb3.pm <- crossbasis(chicagoNMMAPS$pm10, lag=1, argvar=list(fun="lin"), arglag=list(fun="strata"))

cb3.temp <- crossbasis(df$t2m, lag=30, 
                       argvar=list(fun="bs", knots=equalknots(df$t2m,fun="bs",df=5,degree=2)), 
                       arglag=list(knots=logknots(30, 3)))
cb3.tp <- crossbasis(df$tp, lag=30, 
                       argvar=list(fun="bs", knots=equalknots(df$tp,fun="bs",df=5,degree=2)), 
                       arglag=list(knots=logknots(30, 3)))
cb3.fg10 <- crossbasis(df$fg10, lag=30, 
                     argvar=list(fun="bs", knots=equalknots(df$fg10,fun="bs",df=5,degree=2)), 
                     arglag=list(knots=logknots(30, 3)))

# Fit model
model3 <- glm(total_patients ~ cb3.temp*cb3.tp + cb3.fg10 + ns(time, 7*14) + weekday + holiday,
              family=quasipoisson(), df)

# Predict
pred3.temp <- crosspred(cb3.temp, model3, cen=10, by=1)
pred3.tp <- crosspred(cb3.tp, model3, cen=0, by=1)
pred3.fg10 <- crosspred(cb3.fg10, model3, cen=30, by=1)

# Contour plots
plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30,
     main="3D graph of temperature effect")
plot(pred3.temp, "contour", xlab="Temperature", key.title=title("RR"),
     plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

# Figure 4a: illustrates lag-response curves specific to mild and extreme cold 
# and hot temperatures of -20◦C, 0◦C, 27◦C, and 33◦C (with reference at 21◦C)
plot(pred3.temp, "slices", var=-5, ci="n", col=1, ylim=c(0.95,1.25), lwd=1.5,
     main="Lag-response curves for different temperatures, ref. 10C")
for(i in 1:3) lines(pred3.temp, "slices", var=c(0, 15, 20)[i], col=i+1, lwd=1.5)
legend("topright",paste("Temperature =",c(-5, 0, 15, 20)), col=1:4, lwd=1.5)

# Figure 4b: exposure-response relationships specific to lag 0 and 5 (left column), 
# and lag-response relationships specific to temperatures -20◦C and 33◦C (right column)
plot(pred3.temp, "slices", var=c(-5,20), lag=c(0,10), col=4,
     ci.arg=list(density=40,col=grey(0.7)))






# ------------------------------------------------------------------------------ 
# ---------------------------- DNLM Development --------------------------------
#  Multivariate model of temperature, pm10 and o3 on CVD and respiratory deaths
#  tutorial: https://cran.r-project.org/web/packages/dlnm/vignettes/dlnmTS.pdf
# ------------------------------------------------------------------------------

# Load data
data(chicagoNMMAPS)

#### LINEAR: Simple DLM (PM10 & TEMP VS DEATH) ####
# Create cross basis functions for each predictor (exposures pm10; control temp)
cb1.pm <- crossbasis(chicagoNMMAPS$pm10, lag=15, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
cb1.temp <- crossbasis(chicagoNMMAPS$temp, lag=3, argvar=list(df=5), arglag=list(fun="strata",breaks=1))


# Fit model 
model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7*14) + dow, 
              family=quasipoisson(), chicagoNMMAPS, na.action=na.omit)
# over dispersed Poisson distribution, including a smooth function of time with 7 df/year
# (in order to correct for seasonality and long time trend) and day of the week as factor

# Given this model, estimate the association of specific levels of PM10 on mortality
pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.2, cumul=TRUE)
# 'at=0:20' states that the prediction must be computed for each integer value from 0 to 20 µgr/m3. 
# By setting bylag=0.2, the prediction is computed along the lag space with an increment of 0.2. 

# Lag−response curves
plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
       main="Association with a 10-unit increase in PM10")
plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
       main="Cumulative association with a 10-unit increase in PM10")

# Get the overall cumulative effect of a 10-unit increase in PM10 over 15 days of lag & CI
pred1.pm$allRRfit["10"]
cbind(pred1.pm$allRRlow, pred1.pm$allRRhigh)["10",]

#### LINEAR: Seasonal analysis (OZONE & TEMP VS DEATH) ####

chicagoNMMAPSseas <- subset(chicagoNMMAPS, month %in% 6:9)

# Create cross basis functions for each predictor (exposures o3 & temp)
cb2.o3 <- crossbasis(chicagoNMMAPSseas$o3, lag=5, argvar=list(fun="thr",thr=40.3), 
                     arglag=list(fun="integer"),group=chicagoNMMAPSseas$year)
cb2.temp <- crossbasis(chicagoNMMAPSseas$temp, lag=10, argvar=list(fun="thr",thr=c(15,25)), 
                       arglag=list(fun="strata",breaks=c(2,6)), group=chicagoNMMAPSseas$year)

# Fit model 
model2 <- glm(death ~ cb2.o3 + cb2.temp + ns(doy, 4) + ns(time,3) + dow,
              family=quasipoisson(), chicagoNMMAPSseas)
pred2.o3 <- crosspred(cb2.o3, model2, at=c(0:65,40.3,50.3))

# Lag response and overal associations
plot(pred2.o3, "slices", var=50.3, ci="bars", type="p", col=2, pch=19,
     ci.level=0.80, main="Lag-response a 10-unit increase above threshold (80CI)")
plot(pred2.o3,"overall",xlab="Ozone", ci="l", col=3, ylim=c(0.9,1.3), lwd=2,
     ci.arg=list(col=1,lty=3), main="Overall cumulative association for 5 lags")

# Get the overall cumulative effect of a 10-unit increase in ozone above 50.3 with CI
pred2.o3$allRRfit["50.3"]
cbind(pred2.o3$allRRlow, pred2.o3$allRRhigh)["50.3",]

#### NON-LINEAR: Bi-dimensional (PM10 & TEMP VS DEATH) ####

# cross-basis for temperature is specied through a natural and non-natural splines
cb3.pm <- crossbasis(chicagoNMMAPS$pm10, lag=1, argvar=list(fun="lin"), arglag=list(fun="strata"))
varknots <- equalknots(chicagoNMMAPS$temp,fun="bs",df=5,degree=2)
lagknots <- logknots(30, 3)
cb3.temp <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(fun="bs", knots=varknots), arglag=list(knots=lagknots))

# Fit model
model3 <- glm(death ~ cb3.pm + cb3.temp + ns(time, 7*14) + dow,
              family=quasipoisson(), chicagoNMMAPS)
pred3.temp <- crosspred(cb3.temp, model3, cen=21, by=1)

# Contour plots
plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30,
       main="3D graph of temperature effect")
plot(pred3.temp, "contour", xlab="Temperature", key.title=title("RR"),
       plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

# Figure 4a: illustrates lag-response curves specific to mild and extreme cold 
# and hot temperatures of -20◦C, 0◦C, 27◦C, and 33◦C (with reference at 21◦C)
plot(pred3.temp, "slices", var=-20, ci="n", col=1, ylim=c(0.95,1.25), lwd=1.5,
       main="Lag-response curves for different temperatures, ref. 21C")
for(i in 1:3) lines(pred3.temp, "slices", var=c(0,27,33)[i], col=i+1, lwd=1.5)
legend("topright",paste("Temperature =",c(-20,0,27,33)), col=1:4, lwd=1.5)

# Figure 4b: exposure-response relationships specific to lag 0 and 5 (left column), 
# and lag-response relationships specific to temperatures -20◦C and 33◦C (right column)
plot(pred3.temp, "slices", var=c(-20,33), lag=c(0,5), col=4,
       ci.arg=list(density=40,col=grey(0.7)))

#### NON-LINEAR: Reducing a DLNM ####

# WHAT IS THE BEENFIT?
# cb4 <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(fun="thr",thr=c(10,25)), arglag=list(knots=lagknots))
# model4 <- glm(death ~ cb4 + ns(time, 7*14) + dow, family=quasipoisson(), chicagoNMMAPS)
# pred4 <- crosspred(cb4, model4, by=1)
# 
# redall <- crossreduce(cb4, model4) # overall cumulative
# redlag <- crossreduce(cb4, model4, type="lag", value=5) # lag specific
# redvar <- crossreduce(cb4, model4, type="var", value=33) # predictor specific








