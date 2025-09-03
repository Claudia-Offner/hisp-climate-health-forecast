### ------------------------------------------------------------------------ ### 
###  DNLM - Sensitivity Analysis                                             ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Path
user <- 'claer14'
file_import <- 'ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_wide'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 2 - Causal Inference/')
setwd(path1)

# Packages
# install.packages('')
library(dlnm); library(splines); library(dplyr); library(lubridate); 
library(MASS); library(patchwork); library(ggplot2)

# Function to extract model diagnostics
get_model_fit <- function(model,name) {
  modelfit <- data.frame(name,
                         summary(model)$df.null+1,                               
                         summary(model)$df[1],
                         summary(model)$dispersio,       
                         summary(model)$deviance, 
                         summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion) 
  names(modelfit) <- c('model','obs','dfmodel','dispersion','deviance','AIC' )
  return(modelfit)
}

# Function to test model parameters for optimal fit 
get_parameter_test <- function(variables, v_dist=NULL, v_lag=NULL, time_df=9, baseline=T){
  
  # Set baseline time parameters
  years <- length(unique(df$year))
  baseformula <- paste0('total_patients ~  season + weekday + ns(time, df=', time_df,'*years)')
  
  if(baseline==TRUE){
    # Only fit the baseline model
    model <- glm.nb(formula=as.formula(baseformula), data=df, na.action='na.exclude') 
    test_result <- get_model_fit(model,paste0('baseline(',time_df,')'))
  } else{
    
    if(length(variables)==1){
      # Fit univariate model
      cb <- crossbasis(df[[variables]], lag=c(0,v_lag),argvar=v_dist, arglag=list(fun='ns', df=4), group=df$year)
      model <- glm.nb(formula=as.formula(paste(baseformula,'+ cb')),data=df, na.action='na.exclude') 
      test_result <- get_model_fit(model,paste0(variables,'(',v_dist$fun,'[',v_dist$df , '], ', v_lag,')'))
      
    } else{
      # Fit multivariate model
      cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
      cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
      cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
      model <- glm.nb(formula=as.formula(paste(baseformula,'+ cb.t2m + cb.tp + cb.fg10')),data=df, na.action='na.exclude') 
      test_result <- get_model_fit(model,paste0('autocorr(', 1 ,')'))
    }
    
  }
  
  return(list(test=test_result, model=model))
  
}

# Function to test outcome correlation lags
get_autocorr_test <- function(lag, df){
  
  if(lag==0){
    # Fit base model (no lag)
    cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
    cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
    cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
    model <- glm(formula=as.formula('total_patients ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10'),family=quasipoisson,data=df, na.action='na.exclude') 
    test_result <- get_model_fit(model,paste0('autocorr(', 0 ,')'))
    
  } else{
    # Add lags to df
    df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, lag))
    
    # Fit multivariate model
    cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
    cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
    cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
    model <- glm(formula=as.formula('total_patients ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag'),family=quasipoisson,data=df, na.action='na.exclude') 
    test_result <- get_model_fit(model,paste0('autocorr(', lag ,')'))
    
  }
  
  return(list(test=test_result, model=model))
}

# Function to retrieve results from each predictor, with visuals and tables
get_results <- function(df, outcome){
  
  # Function: Get cumulative exposure-response and lag-response effects
  get_visuals <- function(cp, model, xlab, title, col, lag_rep=FALSE, lag=NULL){
    
    if(lag_rep==FALSE){
      if(!is.null(lag)){
        cum <- data.frame(x = cp$predvar, RR = cp$matRRfit[, lag + 1], lower = cp$matRRlow[, lag + 1], upper = cp$matRRhigh[, lag + 1], row.names = NULL)
      } else {
        cum <- data.frame(x = cp$predvar, RR = cp$allRRfit, lower = cp$allRRlow, upper = cp$allRRhigh, row.names = NULL)
      }
      p <- ggplot(cum, aes(x = x, y = RR)) +
        geom_line(color = col) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, fill = col) +
        geom_hline(yintercept = 1, color = "gray20") +
        scale_x_continuous(breaks = pretty(cum$x, n = 6)) +
        # scale_y_continuous(limits = c(-0.5, 2.5)) +
        labs(x = xlab, y = "RR",
             title = title) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),      # pale axis lines
          axis.text = element_text(color = "gray20"),      # pale axis labels
          axis.ticks = element_line(color = "gray20"),     # pale ticks
          plot.title = element_text(face = "bold", hjust = 0.5)  # center title
        )  
    } else {
      
      slice <- data.frame(lag=seq(cp$lag[1], cp$lag[2]), RR=t(cp$matRRfit), lower=t(cp$matRRlow), upper=t(cp$matRRhigh), row.names = NULL)
      names(slice) <- c('lag', 'RR', 'lower','upper')
      p <- ggplot(slice, aes(x = lag, y = RR)) +
        geom_hline(yintercept = 1, color = "gray20") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,  size = 0.2) +
        geom_point(color = col,  size = 2.5) +
        scale_x_continuous(breaks = seq(min(slice$lag), max(slice$lag), by = 2)) +  # more breaks
        labs(x = "Lag (days)", y = "RR",
             title = title) +
        theme_classic() +
        theme(
          axis.line = element_line(color = "gray20"),      # pale axis lines
          axis.text = element_text(color = "gray20"),      # pale axis labels
          axis.ticks = element_line(color = "gray20"),     # pale ticks
          plot.title = element_text(face = "bold", hjust = 0.5)  # center title
        )
      
    }
    
    return(p)
  }
  
  # Fit cross-basis
  cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=1), arglag=list(fun='ns', df=4), group=df$year)
  cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='lin'), arglag=list(fun='ns', df=4), group=df$year)
  vars <- list(t2m = cb.t2m, tp = cb.tp, fg10 = cb.fg10)
  
  # Fit model
  model <- glm.nb(formula=as.formula(paste0(outcome, ' ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag')), data = df)

  # Select quantile
  q_t2m <- quantile(df$t2m, probs = 0.95, na.rm = TRUE)
  q_tp <- quantile(df$tp, probs = 0.95, na.rm = TRUE)
  q_fg10 <- quantile(df$fg10, probs = 0.95, na.rm = TRUE)
  
  # DISPLAY ESTIMATE OF RR FOR BENCHMARK TP INCREMENT OF 23.4 (99th perc) 
  cp.t2m <- crosspred(cb.t2m, model, at=q_t2m)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.tp <- crosspred(cb.tp, model, at=q_tp)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.fg10 <- crosspred(cb.fg10, model, at=q_fg10)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp <- data.frame(rbind(cp.t2m, cp.tp, cp.fg10))
  df_cp <- data.frame(variables = rownames(cp), cp, row.names = NULL)
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK T2M value
  cp <- crosspred(cb.t2m, model, at=q_t2m)
  p_lr_t2m <- get_visuals(cp, model, lag_rep=TRUE, col="#DF536B", title='Lag-response of t2m (at 95th percentile)')
  # GRAPH CUMULATIVE T2M-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & FG10) 
  cp <- crosspred(cb.t2m, model, at=seq(-10, 35, by=1))
  p_er_t2m <- get_visuals(cp, model, lag_rep=FALSE, lag=5, col="#DF536B", xlab="Temperature (C)", title="RR of t2m exposure-response (lag=5)")

  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK TP value
  cp <- crosspred(cb.tp, model, at=q_tp)
  p_lr_tp <- get_visuals(cp, model, lag_rep=TRUE, col="#1D98EB", title='Lag-response of tp (at 95th percentile)')
  # GRAPH CUMULATIVE TP-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR T2M & FG10) 
  cp <- crosspred(cb.tp, model)
  p_er_tp <- get_visuals(cp, model, lag_rep=FALSE, lag=6, col="#1D98EB", xlab="Rainfall (mm)", title="RR of tp exposure-response (lag=6)")

  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK FG10 value
  cp <- crosspred(cb.fg10, model, at=q_fg10)
  p_lr_fg10 <- get_visuals(cp, model, lag_rep=TRUE, col="#A6B727", title='Lag-response of fg10 (at 95th percentile)')
  # GRAPH CUMULATIVE FG10-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & T2M) 
  cp <- crosspred(cb.fg10, model)
  p_er_fg10 <- get_visuals(cp, model, lag_rep=FALSE, lag=NULL, col="#A6B727", xlab="Wind speed (km/hr)", title="Overall RR of fg10 exposure-response")

  result <- list(model=model, overall=df_cp, 
                 exp_rep = list(t2m=p_er_t2m, tp=p_er_tp, fg10=p_er_fg10),
                 lag_rep = list(t2m=p_lr_t2m, tp=p_lr_tp, fg10=p_lr_fg10))
  
  return(result)
  
}

# Function to format and organizing visuals 
get_all_plots <- function(plot_list, single_outcome=NULL, rows=3){
  
  # Function: Format plots for visualizing (2*4 grid)
  multi_plot_format <- function(plot_list){
    
    p <- plot_list
    
    # Loop to remove existing titles and replace with the list name
    p <- lapply(names(p), function(n) {
      p[[n]] + labs(title = n)
    })
    
    # Keep names after lapply
    names(p) <- names(p)
    
    # result <- (p[[1]] | p[[2]] | p[[3]] | p[[4]]) / (p[[5]] | p[[6]] | p[[7]] | p[[8]])
    
    result <- wrap_plots(p, nrow = rows)

    return(result)
  }
  
  
  # Visualise all outputs
  if(!is.null(single_outcome)){
    
    # Arrange your 2x3 grid
    exp_rep <- (plots$exp_rep_t2m[[single_outcome]] | plots$exp_rep_tp[[single_outcome]] | plots$exp_rep_fg10[[single_outcome]]) + plot_annotation(title = paste0('Outcome: ', single_outcome), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep <- (plots$lag_rep_t2m[[single_outcome]] | plots$lag_rep_tp[[single_outcome]] | plots$lag_rep_fg10[[single_outcome]]) + plot_annotation(title = paste0('Outcome: ', single_outcome), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    all <- (exp_rep / lag_rep) + plot_annotation(title = paste0('Outcome: ', single_outcome), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
      
    # Add a main title
    grid_plots <- list(all=all, exp_rep=exp_rep, lag_rep=lag_rep)
    
  } else {
    
    exp_rep_t2m <- multi_plot_format(plot_list$exp_rep_t2m) + plot_annotation(title = 'RR of t2m exposure-response (lag=5)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep_t2m <- multi_plot_format(plot_list$lag_rep_t2m) + plot_annotation(title = 'Lag-response of t2m (at 95th percentile)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    exp_rep_tp <- multi_plot_format(plot_list$exp_rep_tp) + plot_annotation(title = 'RR of tp exposure-response (lag=6)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep_tp <- multi_plot_format(plot_list$lag_rep_tp) + plot_annotation(title = 'Lag-response of tp (at 95th percentile)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    exp_rep_fg10 <- multi_plot_format(plot_list$exp_rep_fg10) + plot_annotation(title = 'Overall RR of fg10 exposure-response', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep_fg10 <- multi_plot_format(plot_list$lag_rep_fg10) + plot_annotation(title = 'Lag-response of fg10 (at 95th percentile)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    
    grid_plots <- list(exp_rep_t2m=exp_rep_t2m, lag_rep_t2m=lag_rep_t2m, exp_rep_tp=exp_rep_tp, 
                       lag_rep_tp=lag_rep_tp, exp_rep_fg10=exp_rep_fg10, lag_rep_fg10=lag_rep_fg10)
    
  }
  
  return(grid_plots)
  
}


#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/', file_import, '.RData'))

#### Test time df for capturing long-term/seasonal (time_df=9) ####
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 1:10){
  test_result <- get_parameter_test(time_df=i, baseline=T)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)


#### Test single weather variables for distribution & lag structures ####

# Temperature (lag=1:14, argvar=ns[df=3]) NB: updated data reflect this and so does lit
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 1:14){ # test linear
  test_result <- get_parameter_test(variables=c('t2m'), v_dist=list(fun='lin'), v_lag=i, baseline=F)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
for(i in 1:14){ # test non-linear
  for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
    test_result <- get_parameter_test(variables=c('t2m'), v_dist=list(fun='ns', df=j), v_lag=i, baseline=F)
    modelfitstats <- rbind(modelfitstats, test_result$test)
  }
}
(modelfitstats)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)

# Rainfall (lag=1:14, argvar=lin) - very similar to ns but opt for simpler model
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 1:14){ # test linear
  test_result <- get_parameter_test(variables=c('tp'), v_dist=list(fun='lin'), v_lag=i, baseline=F)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
for(i in 1:14){ # test non-linear
  for (j in 3:8){ # test exp-resp dfs
    test_result <- get_parameter_test(variables=c('tp'), v_dist=list(fun='ns', df=j), v_lag=i, baseline=F)
    modelfitstats <- rbind(modelfitstats, test_result$test)
  }
}
(modelfitstats)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)

# Wind speed (lag=1:14, argvar=lin) - very similar to ns but opt for simpler model
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 1:14){ # test linear
  test_result <- get_parameter_test(variables=c('fg10'), v_dist=list(fun='lin'), v_lag=i, baseline=F)
  modelfitstats <- rbind(modelfitstats, test_result$test)
}
for(i in 1:14){ # test non-linear
  for (j in 3:8){ # test exp-resp dfs
    test_result <- get_parameter_test(variables=c('fg10'), v_dist=list(fun='ns', df=j), v_lag=i, baseline=F)
    modelfitstats <- rbind(modelfitstats, test_result$test)
  }
}
(modelfitstats)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)

#### Test autocorrelation of hospital admissions (AR=1) #####

setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
for(i in 0:7){
  test_result <- get_autocorr_test(lag=i, df)
  modelfitstats <- rbind(modelfitstats, test_result$test)
  # Check autocorrelation
  resid <- residuals(test_result$model, type = "pearson")
  resid <- resid[!is.na(resid)]
  acf(resid, main="ACF of DLNM Model Residuals")
  pacf(resid, main="PACF of DLNM Model Residuals")
}
modelfitstats
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)

#### Test shape of response (negative binom) #####
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))

# Add lags to df
df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))

# Fit multivariate model (quasipoisson)
cb.t2m <- crossbasis(df$mx2t, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=4), group=df$year)
cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=7), arglag=list(fun='ns', df=4), group=df$year)
model <- glm(formula=as.formula('total_patients ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag'),family=quasipoisson,data=df, na.action='na.exclude') 
test_result <- get_model_fit(model,'final(quasipoisson)')
modelfitstats <- rbind(modelfitstats, test_result)

# Fit multivariate model (negative binomial)
model <- glm.nb(formula=as.formula('total_patients ~  season + weekday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag'),
                data = df)
test_result <- get_model_fit(model,'final(neg binom)')
modelfitstats <- rbind(modelfitstats, test_result)

#### Check differences across responses ####

# Add lags to df
df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))

# Reinitialize dfs and lists
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
exp_rep_t2m <- list()
lag_rep_t2m <- list()
exp_rep_tp <- list()
lag_rep_tp <- list()
exp_rep_fg10 <- list()
lag_rep_fg10 <- list()

# Extract data for all outcomes
outcomes <- c('total_patients', 'F', 'M', 'under_5', 'age5_14', 'age15_24', 'age25_64', 'over_65', 
              'general', 'ophthalmology', 'pediatrics', 'ent')

for (i in outcomes) {

  results <- get_results(df, i)
  # Get model fit
  test_result <- get_model_fit(results$model, i)
  modelfitstats <- rbind(modelfitstats, test_result)
  
  # Get model plots
  exp_rep_t2m[[i]] <- results$exp_rep$t2m
  lag_rep_t2m[[i]] <- results$lag_rep$t2m
  exp_rep_tp[[i]] <- results$exp_rep$tp
  lag_rep_tp[[i]] <- results$lag_rep$tp
  exp_rep_fg10[[i]] <- results$exp_rep$fg10
  lag_rep_fg10[[i]] <- results$lag_rep$fg10

}

plots <- list(exp_rep_t2m=exp_rep_t2m, lag_rep_t2m=lag_rep_t2m, exp_rep_tp=exp_rep_tp, 
              lag_rep_tp=lag_rep_tp, exp_rep_fg10=exp_rep_fg10, lag_rep_fg10=lag_rep_fg10)


(total_patients <- get_all_plots(total_patients, single_outcome = 'total_patients'))
(all_results <- get_all_plots(plots))

#### Export outcome plots #####
for (nm in names(total_patients)) {
  ggsave(filename = paste0(path3, '/IV. Figures/total_patients_', nm, '.png'), plot = total_patients[[nm]],
         width = 12, height = 6, dpi = 300)
}
for (nm in names(all_results)) {
  ggsave(filename = paste0(path3, '/IV. Figures/', nm, '.png'), plot = all_results[[nm]],
         width = 9, height = 6, dpi = 300)
}


