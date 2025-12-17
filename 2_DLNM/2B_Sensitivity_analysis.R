### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city Sensitivity  Analysis                                 ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

#### FUNCTIONS ####

# Function for parameter testing
get_parameter_test <- function(data, locations, reference='50%',
                               test='Variable Lag', variables=NULL, 
                               v_lag=14, time_df=8, auto_lag=1,
                               arg_var=list(fun='ns', df=3), 
                               arg_lag=list(fun='ns', df=3)){
  
  # SET PARAMETERS
  datalist <- lapply(locations, function(region) data[data$location==region,])
  names(datalist) <- locations
  years <- length(unique(data$year))
  
  # EMPTY OBJECTS TO STORE OUTPUTS
  list_comp <- vector('list',m)
  list_model <- vector('list',m)
  names(list_comp) <- locations
  names(list_model) <- locations
  
  # IF variable = test parameters
  if(!is.null(variables)){
    
    ranges <- t(sapply(datalist,function(x) range(x[[variables]],na.rm=T)))
    cen <- sapply(datalist,function(x) {
      quantile(x[[variables]],c(0:100)/100,na.rm=T)
    })
    arg <- list(argvar=arg_var$df, arglag=arg_lag$df)
    df <- prod(unlist(arg))
    
    # Set storage matrix
    ymat <- matrix(NA,m,df,dimnames=list(locations,paste("spl",seq(df),sep="")))
    Slist <- vector("list",m)
    names(Slist) <- locations
    
  }
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  system.time(
    for(i in seq(m)) {
      
      # PRINT ITERATION
      cat(i,'')
      
      # LOAD
      data <- datalist[[i]]

      # RUN THE MODEL & EXTRACT COMPARATIVE METRICS
      if(test=='Longterm'){
        # Fit model
        baseformula <- as.formula(paste0('total_patients ~  season + weekday + ns(time, df=', time_df,'*', years, ')'))
        model <- glm(formula=baseformula,family=quasipoisson(),data)
        list_comp[[i]] <- comp_metrics(model)
        list_model[[i]] <- model
      } else if (test == 'Autocorrelation') {
        
        if (auto_lag==0) {
          # Fit model
          baseformula <- as.formula(paste0('total_patients ~  season + weekday +  ns(time, df=', time_df,'*', years, ')'))
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model)
          list_model[[i]] <- model
        } else{
          # Add lags to df
          data <- data %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, auto_lag))
          # Fit model
          baseformula <- as.formula(paste0('total_patients ~  season + weekday + ns(time, df=', time_df,'*', years, ') + patients_lag '))
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model)
          list_model[[i]] <- model
        }
      } else if (grepl('Variable', test)){
        # Add lags to df
        data <- data %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, auto_lag))
        # Fit cross basis
        cb <- crossbasis(data[[variables]], lag=c(0,v_lag),argvar=arg_var, arglag=arg_lag, group=data$year)
        # Fit model
        baseformula <- as.formula(paste0('total_patients ~  season + weekday + ns(time, df=', time_df,'*', years, ') + cb + patients_lag'))
        model <- glm(formula=as.formula(baseformula),family=quasipoisson(),data)
        list_comp[[i]] <- comp_metrics(model, cb=cb, cen=cen[reference,i])
        list_model[[i]] <- model
      }
    })
  
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['sum_aic']] <- sum(sapply(list_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(list_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']] <- rowMeans(sapply(list_comp, function(x) as.numeric(x$acf)))
  
  if (grepl('Variable', test)) {
    
    meta_comp[['avg']] <- mean(sapply(list_comp, function(x) x$esti))
    
  }
  
  # Format AIC
  if(test=='Longterm'){
    aic_result <- data.frame(paste0(test, ' (',time_df,')'), meta_comp[['sum_aic']])
  } else if (test == 'Autocorrelation') {
    aic_result <- data.frame(paste0(test, ' (',auto_lag,')'), meta_comp[['sum_aic']])
  } else if (test == 'Variable Lag'){
    aic_result <- data.frame(paste0(variables, 
                                    ' (arglag=', arg_lag$fun, '|', arg_lag$df,  ')'), meta_comp[['sum_aic']])
  } else if (test == 'Variable Shape'){
    aic_result <- data.frame(paste0(variables, 
                                    ' (argvar=', arg_var$fun, '|', arg_var$df,')'), meta_comp[['sum_aic']])
  }
  names(aic_result) <- c('model','AIC' )
  
  return(list(aic_result=aic_result, models=list_model,
              pacf_result=meta_comp[['sum_pacf']], 
              acf_result=meta_comp[['sum_acf']]))
  
}

# Function to plot AIC results (comparison)
plot_AIC_res <- function(modelfitstats, best_fit, main_title){
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, 'yes', 'no')
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels = modelfitstats$model)
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x = model, y = AIC, group = 1)) +
    geom_line(color = 'grey20') +
    geom_point(aes(color = highlight), size = 3) +
    scale_color_manual(values = c('no' = 'black', 'yes' = 'red')) +
    labs(title = main_title,
         x = 'Model',
         y = 'Pooled AIC') +
    theme_classic() +
    theme(
      axis.line = element_line(color = 'gray20'),
      axis.text = element_text(color = 'gray20'),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_line(color = 'gray20'),
      plot.title = element_text(face = 'bold', hjust = 0.5, size=12),
      legend.position = 'none'  # hide legend if not needed
    )
  return(p)
  
}

# Function to plot pooled ACF/PACF results
plot_pooled_acf <- function(test_result, plot='PACF'){
  
  # Conditionally set plot
  value_col <- if(plot=='ACF') test_result$acf_result else test_result$pacf_result
  if(plot=='ACF') value_col <- value_col[-1]   # Skip lag 0 for ACF
  
  # Set data frame
  df <- data.frame(Lag = 1:length(value_col), Value = value_col)
  n_obs <- sum(sapply(test_result$models, function(m) length(m$residuals)))

  # Plot
  ggplot(df, aes(x = Lag, y = Value)) +
    geom_segment(aes(x = Lag, xend = Lag, y = 0, yend = Value), linewidth = 0.4) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_hline(yintercept = c(-1.96 / sqrt(n_obs), 1.96 / sqrt(n_obs)),
               linetype = 'dashed', color = 'blue') +
    labs(title = paste0('Pooled ', plot), x = 'Lag', y = plot) +
    theme_classic() +
    theme(
      axis.line = element_line(color = 'black'),
      axis.ticks = element_line(color = 'black'),
      axis.text = element_text(color = 'black'),
      plot.title = element_text(size = 14, face = 'bold', hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(-0.075, 0.1), xlim = c(0, 30)) +
    scale_x_continuous(breaks = round(seq(0, max(df$Lag), length.out = 6)))
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)


# 1. (Longterm) Test time df for capturing long-term/seasonal (time_df=14; AIC_season) ####

setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 1:10){
  test_result <- get_parameter_test(data=df, locations=locations, test='Longterm', time_df=i)
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_season <- plot_AIC_res(modelfitstats, best_fit, 'Pooled AIC: Long-term Seasonality'))


# 2. (Autocorr) Test autocorrelation of hospital admissions (AR=1; AIC_autocorr) #####

# Visualize PACF plot
test_result <- get_parameter_test(data=df, locations=locations, test='Autocorrelation', auto_lag=0)
acf_plot <- plot_pooled_acf(test_result, plot='ACF')
pacf_plot <- plot_pooled_acf(test_result, plot='PACF')
(AIC_ACF <- acf_plot + pacf_plot + plot_layout(ncol = 2))

# Assess AR lag term
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 0:14){
  test_result <- get_parameter_test(data=df, locations=locations, test='Autocorrelation', auto_lag=i)
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_autocorr <- plot_AIC_res(modelfitstats, best_fit, 'Pooled AIC: Autocorrelation lags'))

# Visualize PACF plot
test_result <- get_parameter_test(data=df, locations=locations, test='Autocorrelation', auto_lag=1)
acf_plot <- plot_pooled_acf(test_result, plot='ACF')
pacf_plot <- plot_pooled_acf(test_result, plot='PACF')
(AIC_ACF <- acf_plot + pacf_plot + plot_layout(ncol = 2))


# 3. (Lag) Test single weather variables for lag structures ####

# STATIC CHOICE: We are interested in two week periods (i.e. 14 days)

# TEMPERATURE (lag=6)
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 3:8){
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Lag' ,
                                    variables=c('t2m'), 
                                    # arg_var=list(fun='lin'), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_t2m_lag <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Temperature lags"))

# RAINFALL (lag=3)
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 3:8){
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Lag' ,
                                    variables=c('tp'), 
                                    # arg_var=list(fun='lin'), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_tp_lag <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Rainfall lags"))

# HUMIDITY (lag=5)
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 3:8){
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Lag' ,
                                    variables=c('rh'), 
                                    # arg_var=list(fun='lin'), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_rh_lag <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Relatve humidity lags"))

# WIND SPEED (lag=3)
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 3:8){
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Lag' ,
                                    variables=c('fg10'), 
                                    # arg_var=list(fun='lin'), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_fg10_lag <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Windspeed lags"))


# AIR POLLUTION (lag=3)
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 3:8){
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Lag' ,
                                    variables=c('pm25'), 
                                    # arg_var=list(fun='lin'), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_pm25_lag <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: PM2.5 lags"))


# Compile plots
(AIC_lags <- AIC_t2m_lag + AIC_tp_lag + AIC_rh_lag + AIC_fg10_lag + AIC_pm25_lag + plot_layout(ncol = 3))


# 4. (Shape) Test single weather variables for shape structures ####

# TEMPERATURE (ns, df=4)
i <- 6
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                  variables=c('t2m'), 
                                  arg_var=list(fun='lin'), 
                                  arg_lag=list(fun='ns', df=i))
modelfitstats <- rbind(modelfitstats, test_result$aic_result)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                    variables=c('t2m'), 
                                    arg_var=list(fun='ns', df=j), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_t2m_shp <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Temperature shape"))

# RAINFALL (ns, df=3)
i <- 3
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                  variables=c('tp'), 
                                  arg_var=list(fun='lin'), 
                                  arg_lag=list(fun='ns', df=i))
modelfitstats <- rbind(modelfitstats, test_result$aic_result)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                    variables=c('tp'), 
                                    arg_var=list(fun='ns', df=j), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_tp_shp <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Rainfall shape"))

# HUMIDITY (lin)
i <- 5
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                  variables=c('rh'), 
                                  arg_var=list(fun='lin'), 
                                  arg_lag=list(fun='ns', df=i))
modelfitstats <- rbind(modelfitstats, test_result$aic_result)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                    variables=c('rh'), 
                                    arg_var=list(fun='ns', df=j), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_rh_shp <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Relative humidity shape"))

# Wind speed (lin)
i <- 3
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                  variables=c('fg10'), 
                                  arg_var=list(fun='lin'), 
                                  arg_lag=list(fun='ns', df=i))
modelfitstats <- rbind(modelfitstats, test_result$aic_result)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                    variables=c('fg10'), 
                                    arg_var=list(fun='ns', df=j), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_fg10_shp <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: Windspeed shape"))


# AIR POLLUTION (lin)
i <- 3
setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                  variables=c('pm25'), 
                                  arg_var=list(fun='lin'), 
                                  arg_lag=list(fun='ns', df=i))
modelfitstats <- rbind(modelfitstats, test_result$aic_result)
for (j in 3:8){ # test exp-resp dfs (between 3 and 8 as per lit)
  test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                    variables=c('pm25'), 
                                    arg_var=list(fun='ns', df=j), 
                                    arg_lag=list(fun='ns', df=i))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats %>% arrange(AIC)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_pm25_shp <- plot_AIC_res(modelfitstats, best_fit, "Pooled AIC: PM2.5 shape"))


# Compile plots
(AIC_shape <- AIC_t2m_shp + AIC_tp_shp + AIC_rh_shp + AIC_fg10_shp + AIC_pm25_shp + plot_layout(ncol=3))


# 5. (Shape) Test shape of response (Quasi) #####

# Quasi poisson
primary <- 'total_patients'
quasi_test <- run_model(df, locations, primary, fam_type='quasi')
quasi_aic <- data.frame(model = 'quasi', AIC = quasi_test$meta_comp$sum_aic)
quasi_res <- format_results(quasi_test, primary)
qp_plot <- quasi_res$p_er_overall
# Negative binomial
binom_test <- run_model(df, locations, primary, fam_type='binom')
binom_aic <- data.frame(model = 'binom', AIC = binom_test$meta_comp$sum_aic)
binom_res <- format_results(binom_test, primary)
nb_plot <- binom_res$p_er_overall
# Checks
modelfitstats <- rbind(quasi_aic, binom_aic)
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)

# Stack and annotate plot
# qp_plot_titled <- qp_plot[[1]] + plot_annotation(title = "Quasipoisson")
# nb_plot_titled <- nb_plot[[1]] + plot_annotation(title = "Negative binomial")
qp_plot_titled <- qp_plot + plot_annotation(title = "Quasipoisson")
nb_plot_titled <- nb_plot + plot_annotation(title = "Negative binomial")
qp_wrapped <- wrap_elements(full = qp_plot_titled)
nb_wrapped <- wrap_elements(full = nb_plot_titled)
(AIC_model <- qp_wrapped + nb_wrapped)

#### EXPORT AIC plots ####
sens_analysis <- list(AIC_season=AIC_season, AIC_autocorr=AIC_autocorr, #AIC_ACF=AIC_ACF, 
                      AIC_lags=AIC_lags, AIC_shape=AIC_shape, AIC_model=AIC_model)

for (nm in names(sens_analysis)) {
  ggsave(filename = paste0(path_ch2, '/IV. Figures/SM_', nm, '.png'), plot = sens_analysis[[nm]],
         width = 10, height = 7, dpi = 300)
  ggsave(filename = paste0(path_report, '- figures/SM_', nm, '.png'), plot = sens_analysis[[nm]],
         width = 10, height = 7, dpi = 300)
}







test_result <- get_parameter_test(data=df, locations=locations, test='Variable Shape' ,
                                  variables=c('t2m'), 
                                  arg_var=list(fun='ns', df=4), 
                                  arg_lag=list(fun='ns', df=3))
main_res <- format_results(test_result)
(total_er_overall <- main_res$p_er_overall + plot_annotation(title = replace_exact_matches(primary, nm)))
(total_lr <- main_res$p_lr + plot_annotation(title = replace_exact_matches(primary, nm)))
(total_er <- main_res$p_er + plot_annotation(title = replace_exact_matches(primary, nm)))

