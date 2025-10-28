setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# Set data path
setwd(path_git)

# install.packages(c("dlnm", "mvmeta", "splines", "dplyr"))  # if needed
library(dlnm)
library(mvmeta)
library(splines)
library(dplyr)

## ---- 1. Data ----

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

## ---- 2. First-stage models (per location) ----

outcome <- 'total_patients'
w_cols <- w_cols[w_cols != "sh_mean"]


fit_one <- function(df, outcome) {
  
  df <- df %>%
    group_by(location) %>%            # ensure lags are per city
    arrange(date) %>%                  # make sure data is sorted by time
    mutate(
      patient_lag = lag(.data[[outcome]], n = 1),
    ) %>%
    ungroup()
  
  # Fit cross-basis
  cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.mn2t <- crossbasis(df$mn2t, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.mx2t <- crossbasis(df$mx2t, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df$year)
  cb.rh_mean <- crossbasis(df$rh_mean, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df$year)
  cb.pm25 <- crossbasis(df$pm25, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df$year)
  
  # Fit model with crossbases in environment
  fmla <- as.formula(paste0(outcome, ' ~  season + weekday + holiday + ns(time, df=9*9) + 
                            cb.t2m + cb.mn2t + cb.mx2t +cb.tp + cb.fg10 + cb.rh_mean + cb.pm25 + patient_lag'))
  model <- glm(formula=fmla, data = df, family = quasipoisson(link="log"))
  
  # Reduce
  rd.t2m <- crossreduce(cb.t2m, model, type = "overall")
  rd.mn2t  <- crossreduce(cb.mn2t, model, type = "overall")
  rd.mx2t <- crossreduce(cb.mx2t, model, type = "overall")
  rd.tp <- crossreduce(cb.tp, model, type = "overall")
  rd.fg10  <- crossreduce(cb.fg10, model, type = "overall")
  rd.rh_mean <- crossreduce(cb.rh_mean, model, type = "overall")  
  rd.pm25 <- crossreduce(cb.pm25, model, type = "overall")  
  
  # Coefficients
  coef_t2m <- coef(model)[colnames(cb.t2m)]
  
  # Covariance
  cb_cols <- grep("^cb\\.t2m", colnames(model.matrix(model)))
  vcov_t2m <- vcov(model)[cb_cols, cb_cols]
  
  
  return(list(model=model, 
              cb=list(t2m=cb.t2m, mn2t=cb.mn2t, mx2t=cb.mx2t, tp=cb.tp, 
                      fg10=cb.fg10, rh_mean=cb.rh_mean, pm25=cb.pm25), 
              reduce=list(t2m=rd.t2m, mn2t=rd.mn2t, mx2t=rd.mx2t, tp=rd.tp, 
                          fg10=rd.fg10, rh_mean=rd.rh_mean, pm25=rd.pm25),
              coef=list(t2m=coef_t2m),
              vcov=list(t2m=vcov_t2m)))
  
}

fits <- df %>%
  group_split(location) %>%
  setNames(unique(df$location)) %>%
  lapply(fit_one, 'total_patients')

## ---- 3. Second-stage meta-analysis ----

pred_name <- "t2m"

# Y: m x k matrix of coefficients
Y <- do.call(rbind, lapply(fits, function(x) coef(x$reduce[[pred_name]])))

# S: k x k x m array of vcov matrices
S_list <- lapply(fits, function(x) vcov(x$reduce[[pred_name]]))
S <- simplify2array(S_list)

# Run meta-analysis
meta_t2m <- mvmeta(Y, S = S, method = "reml")
summary(meta_t2m)

# ------------------------------------------------------


locations <- unique(df$location)
coef_list <- list()
vcov_list <- list()

for (loc in locations) {
  df_loc <- subset(df, location == loc)
  
  # Fit cross-basis
  cb.t2m <- crossbasis(df_loc$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df_loc$year)
  cb.mn2t <- crossbasis(df_loc$mn2t, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df_loc$year)
  cb.mx2t <- crossbasis(df_loc$mx2t, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df_loc$year)
  cb.tp <- crossbasis(df_loc$tp, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df_loc$year)
  cb.fg10 <- crossbasis(df_loc$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df_loc$year)
  cb.rh_mean <- crossbasis(df_loc$rh_mean, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df_loc$year)
  cb.pm25 <- crossbasis(df_loc$pm25, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df_loc$year)
  
  # Fit model with crossbasis in environment
  fmla <- as.formula(paste0(outcome, ' ~  season + weekday + holiday + ns(time, df=9*9) + 
                            cb.t2m + cb.mn2t + cb.mx2t +cb.tp + cb.fg10 + cb.rh_mean + cb.pm25'))
  model <- glm(formula=fmla, data = df_loc, family = quasipoisson(link="log"))
  
  # Extract coefficients and covariance matrices for all predictors
  cb_coefs <- coef(model)[grepl("cb", names(coef(model)))]
  cb_vcov <- vcov(model)[grepl("cb", names(coef(model))), grepl("cb", names(coef(model)))]
  
  coef_list[[loc]] <- cb_coefs
  vcov_list[[loc]] <- cb_vcov
}


