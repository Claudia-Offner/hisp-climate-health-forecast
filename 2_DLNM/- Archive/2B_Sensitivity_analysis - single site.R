### ------------------------------------------------------------------------ ### 
###  DNLM - Sensitivity Analysis                                             ###
### ------------------------------------------------------------------------ ### 

library(gridExtra)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# Set data path
setwd(path_git)

#### FUNCTIONS ####

# Function to test model parameters for optimal fit 
get_parameter_test <- function(variables, v_dist=NULL, v_lag=14, df_lag=4, time_df=9, baseline=T){
  
  # Set baseline time parameters
  years <- length(unique(df$year))
  baseformula <- paste0('total_patients ~  location + season + weekday + holiday + ns(time, df=', time_df,'*years)')
  
  # BASELINE TESTS
  if(baseline==TRUE){
    # Only fit the baseline model
    model <- glm(formula=as.formula(baseformula), data=df,  family = quasipoisson(link="log"), na.action='na.exclude') 
    test_result <- get_model_fit(model,paste0('baseline (',time_df,')'))
  } else{
    
    # LAG/SHAPE TESTS
    if(length(variables)==1){
      # Fit univariate model
      cb <- crossbasis(df[[variables]], lag=c(0,v_lag),argvar=v_dist, arglag=list(fun='ns', df=df_lag))
      model <- glm(formula=as.formula(paste(baseformula,'+ cb')),data=df, family = quasipoisson(link="log"), na.action='na.exclude') 
      test_result <- get_model_fit(model,paste0(variables,'(',v_dist$fun,'[',v_dist$df , '], ', df_lag,')'))
      
      # AUTOCORR TESTS
    } else{
      # Fit multivariate model
      cb.t2m <- crossbasis(df$t2m, lag=c(0,v_lag),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=df_lag))
      cb.tp <- crossbasis(df$tp, lag=c(0,v_lag),argvar=list(fun='lin'), arglag=list(fun='ns', df=df_lag))
      cb.fg10 <- crossbasis(df$fg10, lag=c(0,v_lag),argvar=list(fun='lin'), arglag=list(fun='ns', df=df_lag))
      model <- glm(formula=as.formula(paste(baseformula,'+ cb.t2m + cb.tp + cb.fg10')), data=df, family = quasipoisson(link="log"), na.action='na.exclude') 
      test_result <- get_model_fit(model,paste0('autocorr(', 1 ,')'))
    }
    
  }
  
  return(list(test=test_result, model=model))
  
}

# Function to test outcome correlation lags
get_autocorr_test <- function(lag, df){
  
  if(lag==0){
    # Fit base model (no lag)
    model <- glm(formula=as.formula('total_patients ~  location + season + weekday + holiday + ns(time, df=9*9)'),family=quasipoisson,data=df, na.action='na.exclude') 
    test_result <- get_model_fit(model,paste0('autocorr(', 0 ,')'))
    
  } else{
    # Add lags to df
    df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, lag))
    
    # Fit multivariate model
    model <- glm(formula=as.formula('total_patients ~  season + weekday + holiday + ns(time, df=10*9) + patients_lag'),family=quasipoisson,data=df, na.action='na.exclude') 
    test_result <- get_model_fit(model,paste0('autocorr(', lag ,')'))
    
  }
  
  return(list(test=test_result, model=model))
}

# Function to plot AIC results (comparison)
plot_AIC_res <- function(modelfitstats, best_fit, main_title){
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, "yes", "no")
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels = modelfitstats$model)
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x = model, y = AIC, group = 1)) +
    geom_line(color = "grey20") +
    geom_point(aes(color = highlight), size = 3) +
    scale_color_manual(values = c("no" = "black", "yes" = "red")) +
    labs(title = main_title,
         x = "Model",
         y = "AIC") +
    theme_classic() +
    theme(
      axis.line = element_line(color = "gray20"),
      axis.text = element_text(color = "gray20"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_line(color = "gray20"),
      plot.title = element_text(face = "bold", hjust = 0.5, size=12),
      legend.position = "none"  # hide legend if not needed
    )
  return(p)
  
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# 1. Test time df for capturing long-term/seasonal (time_df=9) ####
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 1:10){
  test_result <- get_parameter_test(time_df=i, baseline=T)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
modelfitstats
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_season <- plot_AIC_res(modelfitstats, best_fit, "Model AIC (seasonality)"))

# 2. Test autocorrelation of hospital admissions (AR=1) #####

# Check autocorrelation
test_result <- get_autocorr_test(lag=0, df)
resid <- residuals(test_result$model, type = "pearson")
resid <- resid[!is.na(resid)]

acf_plot <- autoplot(acf(resid, plot = FALSE)) + 
  ggtitle("ACF of DLNM Model Residuals") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
  scale_y_continuous(limits = c(-0.2, 1)) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "gray20"),
    axis.text = element_text(color = "gray20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "gray20"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"  # hide legend if not needed
  )
pacf_plot <- autoplot(pacf(resid, plot = FALSE)) + 
  ggtitle("PACF of DLNM Model Residuals")+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
  scale_y_continuous(limits = c(-0.2, 1)) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "gray20"),
    axis.text = element_text(color = "gray20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "gray20"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"  # hide legend if not needed
  )

AIC_ACF <- grid.arrange(acf_plot, pacf_plot, ncol = 2)


# Run lag tests
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 0:7){
  test_result <- get_autocorr_test(lag=i, df)
  modelfitstats <- rbind(modelfitstats, test_result$test)
  resid <- residuals(test_result$model, type = "pearson")
  resid <- resid[!is.na(resid)]
  acf(resid, plot = TRUE)
}
best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
(AIC_2_autocorr <- plot_AIC_res(modelfitstats, best_fit, "AIC by Model: Autocorrleation lags"))

# 3. Test single weather variables for lag structures ####

# STATIC CHOICE: We are interested in two week periods (i.e. 14 days)
# 
# # TEMPERATURE (lag=0-2)
# setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
# for(i in 3:8){ # test max lags
#   test_result <- get_parameter_test(variables=c('t2m'), v_dist=list(fun='lin'), df_lag=i, baseline=F)
#   modelfitstats <- rbind(modelfitstats, test_result$test)
# }
# best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
# AIC_t2m_lag <- plot_AIC_res(modelfitstats, best_fit, "AIC by Model: Temperature lags")
# 
# # RAINFALL (lag=0-1)
# setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
# for(i in 3:8){ # test max lags
#   test_result <- get_parameter_test(variables=c('tp'), v_dist=list(fun='lin'), df_lag=i, baseline=F)
#   modelfitstats <- rbind(modelfitstats, test_result$test)
# }
# best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
# AIC_tp_lag <- plot_AIC_res(modelfitstats, best_fit, "AIC by Model: Rainfall lags")
# 
# # WIND SPEED (lag=0-11)
# setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
# for(i in 3:8){ # test max lags
#   test_result <- get_parameter_test(variables=c('fg10'), v_dist=list(fun='lin'), df_lag=i, baseline=F)
#   modelfitstats <- rbind(modelfitstats, test_result$test)
# }
# best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
# AIC_fg10_lag <- plot_AIC_res(modelfitstats, best_fit, "AIC by Model: Wind speed lags")
# 
# # Compile plots
# AIC_3_lags <- grid.arrange(AIC_t2m_lag, AIC_tp_lag, AIC_fg10_lag, ncol = 3)


# 4. Test single weather variables for shape structures ####

# TEMPERATURE (ns, df=3)
i <- 4
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
test_result <- get_parameter_test(variables=c('t2m'), v_dist=list(fun='lin'), df_lag=i, baseline=F)
modelfitstats <- rbind(modelfitstats, test_result$test)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(variables=c('t2m'), v_dist=list(fun='ns', df=j), df_lag=i, baseline=F)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
AIC_t2m_shp <- plot_AIC_res(modelfitstats, best_fit, "Model AIC (t2m shape)")

# RAINFALL (ns, df=4)
# i <- 4
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
test_result <- get_parameter_test(variables=c('tp'), v_dist=list(fun='lin'), df_lag=i, baseline=F)
modelfitstats <- rbind(modelfitstats, test_result$test)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(variables=c('tp'), v_dist=list(fun='ns', df=j), df_lag=i, baseline=F)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
AIC_tp_shp <- plot_AIC_res(modelfitstats, best_fit, "Model AIC (tp shape)")

# WIND SPEED (ns, df=5)
# i <- 5
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
test_result <- get_parameter_test(variables=c('fg10'), v_dist=list(fun='lin'), df_lag=i, baseline=F)
modelfitstats <- rbind(modelfitstats, test_result$test)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(variables=c('fg10'), v_dist=list(fun='ns', df=j), df_lag=i, baseline=F)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
AIC_fg10_shp <- plot_AIC_res(modelfitstats, best_fit, "Model AIC (fg10 shape)")

# Compile plots
AIC_shape <- grid.arrange(AIC_t2m_shp, AIC_tp_shp, AIC_fg10_shp, ncol = 3)

# 5. Test shape of response (negative binom) #####

# Fit appropriate parameters
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))
cb.t2m <- crossbasis(df$mx2t, lag=c(0,14), argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=8))
cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=4))
cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=5))

# Fit multivariate model (quasipoisson)
model <- glm(formula=as.formula('total_patients ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag'),family=quasipoisson,data=df, na.action='na.exclude') 
test_result <- get_model_fit(model,'final(quasipoisson)')
modelfitstats <- rbind(modelfitstats, test_result)
# Get exposure plot for temp only
q_t2m <- quantile(df$t2m, probs = 0.95, na.rm = TRUE)
cp.t2m <- crosspred(cb.t2m, model, at=q_t2m)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
cp <- crosspred(cb.t2m, model, at=seq(-10, 35, by=1))
qp_plot <- get_visuals(cp, model, lag_rep=FALSE, col="#D74B3C", xlab="Temperature (C)", title="Quassi-poisson: Cummulative RR")

# Fit multivariate model (negative binomial)
model <- glm.nb(formula=as.formula('total_patients ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag'),
                data = df)
test_result <- get_model_fit(model,'final(neg binom)')
modelfitstats <- rbind(modelfitstats, test_result)
# Get exposure plot for temp only
q_t2m <- quantile(df$t2m, probs = 0.95, na.rm = TRUE)
cp.t2m <- crosspred(cb.t2m, model, at=q_t2m)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
cp <- crosspred(cb.t2m, model, at=seq(-10, 35, by=1))
nb_plot <- get_visuals(cp, model, lag_rep=FALSE, col="#D74B3C", xlab="Temperature (C)", title="Negative Binomial: Cummulative RR")

AIC_model <- grid.arrange(qp_plot, nb_plot, ncol = 2)


#### EXPORT AIC plots ####
sens_analysis <- list(AIC_season=AIC_season, AIC_ACF=AIC_ACF, 
                      # AIC_3_lags=AIC_3_lags, 
                      AIC_shape=AIC_shape, AIC_model=AIC_model)

for (nm in names(sens_analysis)) {
  ggsave(filename = paste0(path_ch2, '/IV. Figures/SM_', nm, '.png'), plot = sens_analysis[[nm]],
         width = 8, height = 5, dpi = 300)
  ggsave(filename = paste0(path_report, '- figures/SM_', nm, '.png'), plot = sens_analysis[[nm]],
         width = 8, height = 5, dpi = 300)
}






