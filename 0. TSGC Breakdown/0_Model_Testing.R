### ---------------------------------### 
### Predicting Daily National Deaths ###
### ---------------------------------### 

# SET WORKING DIRECTORY
setwd('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/201_heal')

# 0. LOAD DEPENDENCIES ####
library(tsgc)
library(here)
library(zoo)
library(dplyr)
library(tibble)
library(kableExtra)
library(ggplot2)
library(xts)
library(readxl)

# 1. IMPORT & CLEAN DATA ####

# Load data
# df <- read.csv('eng_adm0_ONS_2019-2010_daily-national-mortality.csv')
# df <- read.csv('agy_adm4_NHS-Scotland_2015-2025_daily-national-ae.csv')
df <- read.csv('agy_adm4_NHS-Scotland_2015-2025_daily-hosptial-ae.csv')

# Select hospital
df <- df %>% 
  filter(sitecode=='A111H') %>% # relevant sites
  select(-sitecode)

# Set date type, create XTS object & select outcome
df <- df[ , -2] # Drop irrelevant data
df$date <- as.Date(df$date, format = "%Y-%m-%d")  # Adjust format as needed
Y <- xts(df[ , -1], order.by = df$date)

# Select outcome and subset desired dates
Y <- Y["2017-02-01/2020-02-01"]
Y <- Y$c_NumberOfAttendancesEpisode

# 2. DESCRIPTIVES ####
summary(df$w_NumberOfAttendancesEpisode)

# Calculate the 7 day centered moving average
ma.cent <- zoo::rollmean(diff(Y), 7, align = "center") 

# Create df for plot (actual data vs. MA data)
d <- setNames(cbind(diff(Y), ma.cent), c('New Cases', 'Centered 7-day MA'))
df.f <- data.frame(Date = index(d), 
                   New.Cases = coredata(d[,1]),
                   Centered.7.day.MA = coredata(d[,2]))

# Plot Figure 1
ggplot(data = df.f, aes(x = Date)) +
  geom_line(aes(y = New.Cases, color = "Real data"), linewidth=0.1) +
  geom_line(aes(y = Centered.7.day.MA, color = "Centered 7 day MA"), linewidth=1) +
  scale_y_continuous(n.breaks = 10) +
  xlab("Day") +
  ylab("New cases") +
  scale_x_date(date_breaks = "60 days") +
  scale_color_manual(name='', values=c('Real data'='black', 'Centered 7 day MA'='red')) +
  theme_light(base_size = 12) +
  theme(legend.position = c(.2, 0.85)) +
  theme(legend.title=element_text(size=2),legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(face = "bold"))

# NB: There will be a challenge in modelling the twin dips/peaks of the Christmas holidays
# https://www.ahajournals.org/doi/full/10.1161/01.CIR.0000151424.02045.F7

# 3. PARAMETERS & ESTIMATION ####

# Set sample period (dates)
estimation.date.start <- as.Date("2018-08-01")
estimation.date.end <- as.Date("2018-12-10")

# Set forecast period (days)
n.forecasts <- 14

# Set smoothness, confidence levels & plot visibility
q <- 0.005
confidence.level <- 0.68
plt.length <- 30

# Subset sample period 
y <- Y[estimation.date.start <= index(Y) & index(Y) <= estimation.date.end]

# Fit model with subset & smoothing ratio
model <- SSModelDynamicGompertz$new(Y = y, q = q)
res <- model$estimate()
res$output$model$y

# From the original data, subset dates greater than estimation end date and convert to growth rates
y.eval <- Y %>% 
  subset(index(.) > estimation.date.end) %>% # subset
  tsgc::df2ldl() # log growth rate conversion

res$get_gy_ci(smoothed = TRUE, confidence_level = 0.68)
x <- res$predict_all(n.ahead = 7, return.all = TRUE)
plot_gy_ci(res, plt.start.date = as.Date("2018-01-09"))
plot_gy_components(res, plt.start.date = as.Date("2018-01-09"))
res$predict_level(Y, n.ahead = 7, confidence_level = 0.68)

# 4. RESULTS ####


tsgc::plot_forecast(
  res = res, # projected data (model estimation sample & projection)
  y.eval = y.eval, # actual data (after estimation period)
  n.ahead = n.forecasts,
  plt.start.date = estimation.date.end - plt.length
)

tsgc::plot_new_cases(
  res = res, # projected data (model estimation sample & projection)
  Y = y, # actual data (from estimation period)
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = "%Y-%m-%d",
  plt.start.date = estimation.date.end - plt.length
)

tsgc::plot_holdout(res = res, 
                   Y = y,
                   Y.eval = Y[(estimation.date.end+0:n.forecasts)],
                   confidence.level = 0.68,
                   date_format = "%Y-%m-%d")

