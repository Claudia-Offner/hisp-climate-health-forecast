
install.packages('splines')

library(dlnm)
library(splines)


# Load dataset
data(chicagoNMMAPS)


# ------------------------------------------------------------------------------ 
# ---------------------------- DNLM Development --------------------------------
#  Multivariate model of temperature, pm10 and o3 on CVD and respiratory deaths
# ------------------------------------------------------------------------------

# Create cross basis functions
cb1.pm <- crossbasis(chicagoNMMAPS$pm10, lag=15, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
cb1.temp <- crossbasis(chicagoNMMAPS$temp, lag=3, argvar=list(df=5), arglag=list(fun="strata",breaks=1))

model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7*14) + dow, 
              family=quasipoisson(), chicagoNMMAPS, na.action=na.omit)


# Given this model, estimate the association of specific levels of PM10 on mortality
pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.2, cumul=TRUE)
# 'at=0:20' states that the prediction must be computed for each integer value from 0 to 20 µgr/m3. 
# By setting bylag=0.2, the prediction is computed along the lag space with an increment of 0.2. 

plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
       main="Association with a 10-unit increase in PM10")
plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
       main="Cumulative association with a 10-unit increase in PM10")

# Get the overall cumulative effect of a 10-unit increase in PM10 over 15 days of lag & CI
pred1.pm$allRRfit["10"]
cbind(pred1.pm$allRRlow, pred1.pm$allRRhigh)["10",]




