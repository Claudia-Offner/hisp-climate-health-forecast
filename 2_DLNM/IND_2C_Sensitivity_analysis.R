### ------------------------------------------------------------------------ ### 
###  DNLM: Multiple city Sensitivity Analysis
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()
setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/DLNM - Dependencies.R')

# Load packages
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

#### FUNCTIONS ####

# Function for parameter testing
get_parameter_test <- function(data, locations, outcome, pooled,
                               test='Variable Lag', variables=NULL, v_lag=10, 
                               longterm_df=8, autocorr_lag=14,
                               arg_lag=list(fun='ns', df=3),
                               arg_var=list(fun='ns', df=3)){
  # Build datalist
  datalist <- lapply(locations, function(region) data[data$location == region, ])
  names(datalist) <- locations
  
  # Drop empty datasets
  datalist <- Filter(function(x) nrow(x) > 0, datalist)
  
  # Use names that survived filtering
  locs_used <- names(datalist)
  
  m <- length(datalist)
  years <- length(unique(data$year))
  
  # Initialize output lists with correct names
  list_comp  <- vector('list', m)
  list_model <- vector('list', m)
  
  names(list_comp)  <- locs_used
  names(list_model) <- locs_used
  
  # IF variable=test parameters
  if(!is.null(variables)){
    
    reference <- ifelse(variables %in% c("tp", "fg10"), "10%", "50%")

    ranges <- t(sapply(datalist,function(x) range(x[[variables]],na.rm=T)))
    cen <- sapply(datalist,function(x) {
      quantile(x[[variables]],c(0:100)/100,na.rm=T)
    })
    df <- prod(arg_var$df, arg_lag$df)
    
    # Set storage matrix
    ymat <- matrix(NA,m,df,dimnames=list(locs_used,paste('spl',seq(df),sep='')))
    Slist <- vector('list',m)
    names(Slist) <- locs_used
    
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
        baseformula <- as.formula(paste0(outcome, ' ~  month + weekday + holiday + ns(time, df=', longterm_df,'*', years, ')'))
        model <- glm(formula=baseformula,family=quasipoisson(),data)
        list_comp[[i]] <- comp_metrics(model)
        list_model[[i]] <- model
      } else if (test == 'Autocorrelation') {
        
        if (autocorr_lag==0) {
          # Fit model
          baseformula <- as.formula(paste0(outcome, ' ~  month + weekday + holiday + ns(time, df=', longterm_df,'*', years, ')'))
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model)
          list_model[[i]] <- model
        } else{
          
          # Select lags
          lag_terms <- paste0(outcome, '_lag', 1:autocorr_lag)
          
          baseformula <- as.formula(
            paste0(
              outcome, ' ~ month + weekday + holiday + ns(time, df=', longterm_df, '*', years, ') + ',
              paste(lag_terms, collapse=' + ')
            )
          )
          
          # Fit model
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model)
          list_model[[i]] <- model
        }
      } else if (grepl('Variable', test)){

        # Select lags
        if(autocorr_lag==0){
          lags <- ''
        } else {
          lag_terms <- paste0(outcome, '_lag', 1:autocorr_lag)
          lags <- paste0('+ ', paste(lag_terms, collapse=' + '))
        }

        # Fit cross basis
        cb <- crossbasis(data[[variables]], lag=c(0,v_lag), argvar=arg_var, arglag=arg_lag, group=data$year)
        # Fit model
        baseformula <- as.formula(
          paste0(
            outcome, ' ~ month + weekday + holiday + ns(time, df=', longterm_df, '*', years, ') ', lags, ' + cb'
          )
        )
        model <- glm(formula=as.formula(baseformula),family=quasipoisson(),data)
        list_comp[[i]] <- comp_metrics(model, cb=cb, cen=cen[reference,i])
        list_model[[i]] <- model
        # Make predictions
        pred <- crosspred(cb, model, cen=cen[reference,i])
        ymat[i,] <- pred$coef
        Slist[[i]] <- pred$vcov
      }
    })
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['sum_aic']] <- sum(sapply(list_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(list_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']] <- rowMeans(sapply(list_comp, function(x) as.numeric(x$acf)))
  
  if (grepl('Variable', test)) {
    
    meta_comp[['avg']] <- mean(sapply(list_comp, function(x) x$esti))
    meta_coef <- get_meta_analysis(pooled=pooled, ymat=ymat,  Slist=Slist, arg_var=arg_var, arg_lag=arg_lag, ranges=ranges, cen=cen)
    
  } else {
    meta_coef <- NULL
  }
  
  # Format AIC
  if(test=='Longterm'){
    aic_result <- data.frame(paste0(test, ' (',longterm_df,')'), meta_comp[['sum_aic']])
  } else if (test == 'Autocorrelation') {
    aic_result <- data.frame(paste0(test, ' (',autocorr_lag,')'), meta_comp[['sum_aic']])
  } else if (test == 'Variable Lag'){
    aic_result <- data.frame(paste0(variables, 
                                    ' (arglag=', arg_lag$fun, '|', arg_lag$df,  ')'), meta_comp[['sum_aic']])
  } else if (test == 'Variable Shape'){
    aic_result <- data.frame(paste0(variables, 
                                    ' (argvar=', arg_var$fun, '|', arg_var$df,')'), meta_comp[['sum_aic']])
  }
  names(aic_result) <- c('model','AIC' )
  
  return(list(aic_result=aic_result, models=list_model, meta_coef=meta_coef,
              pacf_result=meta_comp[['sum_pacf']], 
              acf_result=meta_comp[['sum_acf']]))
  
}

# Function to plot AIC results (comparison)
plot_AIC_res <- function(modelfitstats, best_fit, main_title){
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, 'yes', 'no')
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels=modelfitstats$model)
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x=model, y=AIC, group=1)) +
    geom_line(color='grey20') +
    geom_point(aes(color=highlight), size=3) +
    scale_color_manual(values=c('no'='black', 'yes'='red')) +
    labs(title=main_title,
         y='AIC') +
    theme_classic() +
    theme(
      axis.line=element_line(color='gray20'),
      axis.text=element_text(color='gray20'),
      axis.text.x=element_text(angle=45, hjust=1),
      axis.title.x =element_blank(),
      axis.ticks=element_line(color='gray20'),
      plot.title=element_text(face='bold', hjust=0.5, size=12),
      legend.position='none'  # hide legend if not needed
    )
  return(p)
  
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
data <- df

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# Time and store results
sens_res <- list()
start <- Sys.time()
outcome <- 'total_patients'
var <- c('mx2t', 'tp', 'fg10')

# 1. (Longterm) Test time df for capturing long-term/seasonal ####

test_longterm <- function(df, outcome, pooled=TRUE){
  
  if(length(unique(df$location)) !=1 ){
    title <- 'Pooled AIC: Long-term Seasonality'
  } else{
    if(outcome == 'total_patients'){
      title <- as.character(unique(df$location))
    } else {
      title <- 'Emergency AIC: Long-term Seasonality'
    }
  }
  
  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  for(i in 1:10){
    test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled,
                                      test='Longterm', longterm_df=i)
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, best_fit=num))
}

# Pooled tests
sens_res$pooled$longterm$res <- test_longterm(df, outcome, pooled=TRUE)
sens_res$pooled$longterm$aic <- sens_res$pooled$longterm$res$aic_plot

# Site specific tests
city_longterm_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- test_longterm(sub, outcome, pooled=FALSE)
  city_longterm_res[[l]] <- res
  print(l)
}
city_longterm_aic <- wrap_plots(map(city_longterm_res, 'aic_plot'), nrow=2) &
    plot_annotation(title='City-specific AIC: Long-term seasonality',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
sens_res$city$longterm$res <- city_longterm_res
sens_res$city$longterm$aic <- city_longterm_aic

# Emergency ward tests
l <- 'Hamirpur'
sub <- df %>% filter(location==l)
res <- test_longterm(sub, outcome='emergency', pooled=FALSE)
sens_res$emergency$longterm$res <- res
sens_res$emergency$longterm$aic <- res$aic_plot

# EXPORT
sens_analysis <- list('1A_pooled_longterm_AIC'=list(p=sens_res$pooled$longterm$aic, w=6, h=4),
                      '1B_emg_longterm_AIC'=list(p=sens_res$emergency$longterm$aic, w=6, h=4),
                      '1C_city_longterm_AIC'=list(p=sens_res$city$longterm$aic, w=12, h=8))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'figures/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# Check results (update best fit params based on outputs)
sens_res$pooled$longterm$aic
sens_res$emergency$longterm$aic
sens_res$city$longterm$aic

# Updates
updates <- list('Chamba' = 7, 'Dharamshala' = 8, 'Hamirpur' = 9, 'Kullu' = 9,
                'Nahan' = 9, 'Reckong Peo' = 7, 'Solan' = 6, 'Una' = 8)
for(l in locations){
  sens_res$city$longterm$res[[l]]$best_fit <- updates[[l]]
}
sens_res$pooled$longterm$res$best_fit <- 8
sens_res$emergency$longterm$res$best_fit <- 8


# 2. (Autocorr) Test autocorrelation of hospital admissions #####

test_autocorr <- function(df, outcome, pooled=TRUE){
  
  # Function to plot pooled ACF/PACF results
  plot_pooled_acf <- function(test_result, plot='PACF', name=NULL){
    
    # Conditionally set plot
    value_col <- if(plot=='ACF') test_result$acf_result else test_result$pacf_result
    if(plot=='ACF') value_col <- value_col[-1]   # Skip lag 0 for ACF
    
    # Set data frame
    df <- data.frame(Lag=1:length(value_col), Value=value_col)
    n_obs <- sum(sapply(test_result$models, function(m) length(m$residuals)))
    
    # Title
    num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', test_result$aic_result$model))
    
    # Plot
    ggplot(df, aes(x=Lag, y=Value)) +
      geom_segment(aes(x=Lag, xend=Lag, y=0, yend=Value), linewidth=0.4) +
      geom_hline(yintercept=0, color='black') +
      geom_hline(yintercept=c(-1.96 / sqrt(n_obs), 1.96 / sqrt(n_obs)),
                 linetype='dashed', color='blue') +
      labs(title=paste0(name, ' (', plot, ' ', num, ')'), x='Lag', y=plot) +
      theme_classic() +
      theme(
        axis.line=element_line(color='black'),
        axis.ticks=element_line(color='black'),
        axis.text=element_text(color='black'),
        plot.title=element_text(size=12, face='bold', hjust=0.5)
      ) +
      coord_cartesian(ylim = c(-0.2, 0.2), xlim=c(0, 30)) +
      scale_x_continuous(breaks=round(seq(0, max(df$Lag), length.out=6)))
  }
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- 'Pooled AIC: Autocorrelation lags'
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    name <- 'Pooled'
  } else{
    if(outcome == 'total_patients'){
      title <- as.character(unique(df$location))
      longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
      name <- title
    } else {
      title <- 'Emergency AIC: Autocorrelation lags'
      longterm_df <- sens_res$emergency$longterm$res$best_fit
      name <- 'Emergency'
    }
  }
  
  # Visualize ACF/PACF plot
  test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled, 
                                    longterm_df=longterm_df,
                                    test='Autocorrelation', autocorr_lag=0)
  acf_plot_raw <- plot_pooled_acf(test_result, plot='ACF', name=name)
  pacf_plot_raw <- plot_pooled_acf(test_result, plot='PACF', name=name)
  
  # Assess 1st AR lag term
  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  for(i in 0:14){
    test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled, 
                                      longterm_df=longterm_df,
                                      test='Autocorrelation', autocorr_lag=i)
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  # Visualize final PACF plot
  test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled,
                                    longterm_df=longterm_df,
                                    test='Autocorrelation', autocorr_lag=num)
  acf_plot_fin <- plot_pooled_acf(test_result, plot='ACF', name=name)
  pacf_plot_fin <- plot_pooled_acf(test_result, plot='PACF', name=name)
  AIC_ACF <- acf_plot_raw + acf_plot_fin + plot_layout(ncol=2)
  AIC_PACF <- pacf_plot_raw + pacf_plot_fin + plot_layout(ncol=2)
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, best_fit=num, 
              AIC_ACF=AIC_ACF, AIC_PACF=AIC_PACF))
}

# Pooled tests
sens_res$pooled$autocorr$res <- test_autocorr(df, outcome, pooled=TRUE)
sens_res$pooled$autocorr$aic <- sens_res$pooled$autocorr$res$aic_plot
sens_res$pooled$autocorr$acf <- sens_res$pooled$autocorr$res$AIC_ACF

# Site specific tests
city_autocorr_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- test_autocorr(sub, outcome, pooled=FALSE)
  city_autocorr_res[[l]] <- res
  print(l)
}
city_autocorr_aic <- wrap_plots(map(city_autocorr_res, 'aic_plot'), nrow=2) &
    plot_annotation(title='City-specific AIC: Autocorrelation lags',
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
city_autocorr_acf <- (wrap_plots(map(city_autocorr_res, "AIC_ACF"), nrow = 4) & theme(axis.title.x = element_blank()))  +
  plot_annotation(
    title = "City-specific ACF Comparison",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )
sens_res$city$autocorr$res <- city_autocorr_res
sens_res$city$autocorr$aic <- city_autocorr_aic
sens_res$city$autocorr$acf <- city_autocorr_acf

# Emergency ward tests
l <- 'Hamirpur'
sub <- df %>% filter(location==l)
res <- test_autocorr(sub, outcome='emergency', pooled=FALSE)
sens_res$emergency$autocorr$res <- res
sens_res$emergency$autocorr$aic <- res$aic_plot
sens_res$emergency$autocorr$acf <- res$AIC_ACF

# EXPORT
sens_analysis <- list('2A_pooled_autocorr_AIC'=list(p=sens_res$pooled$autocorr$aic, w=6, h=4),
                      '2A_pooled_autocorr_ACF'=list(p=sens_res$pooled$autocorr$acf, w=6, h=4),
                      '2B_emg_autocorr_AIC'=list(p=sens_res$emergency$autocorr$aic, w=6, h=4),
                      '2B_emg_autocorr_ACF'=list(p=sens_res$emergency$autocorr$acf, w=6, h=4),
                      '2C_city_autocorr_AIC'=list(p=sens_res$city$autocorr$aic, w=12, h=8),
                      '2C_city_autocorr_ACF'=list(p=sens_res$city$autocorr$acf, w=12, h=12))

for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'figures/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# Check results (update best fit params based on outputs)
sens_res$pooled$autocorr$aic
sens_res$emergency$autocorr$aic
sens_res$city$autocorr$aic
sens_res$city$autocorr$res$Kullu$aic_plot
sens_res$city$autocorr$res$`Reckong Peo`$aic_plot
sens_res$city$autocorr$res$Solan$aic_plot
sens_res$city$autocorr$res$Una$aic_plot

# Updates
updates <- list('Chamba' = 13, 'Dharamshala' =13, 'Hamirpur' = 11, 'Kullu' = 14,
                'Nahan' = 8, 'Reckong Peo' = 10, 'Solan' = 1, 'Una' = 14)
for(l in locations){
  sens_res$city$autocorr$res[[l]]$best_fit <- updates[[l]]
}
sens_res$pooled$autocorr$res$best_fit <- 14
sens_res$emergency$autocorr$res$best_fit <- 14

# 3. (ArgLag) Test single weather variables for lag structures ####

test_arglag <- function(df, outcome, pooled, var){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- replace_exact_matches(var, nm)
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
  } else{
    if(outcome == 'total_patients'){
      title <- as.character(unique(df$location))
      longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
      autocorr_lag <- sens_res$city$autocorr$res[[title]]$best_fit
    } else {
      title <- replace_exact_matches(var, nm)
      longterm_df <- sens_res$emergency$longterm$res$best_fit
      autocorr_lag <- sens_res$emergency$autocorr$res$best_fit
    }
  }

  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  for(i in 3:6){
    print(i)
    test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled,
                                      longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                                      test='Variable Lag' , variables=c(var), 
                                      arg_lag=list(fun='ns', df=i),
                                      arg_var=list(fun='ns', df=4))
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub('.*\\|(\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, best_fit=num))
}

# Pooled tests
pooled_arglag_res <- list()
for (v in var){
  print(v)
  pooled_arglag_res[[v]] <- test_arglag(df, outcome, pooled=TRUE, v)
}
pooled_arglag_aic <- wrap_plots(map(pooled_arglag_res, 'aic_plot'), nrow=1) &
   plot_annotation(title='Pooled AIC: Variable lag | df',
                   theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
sens_res$pooled$arglag$res <- pooled_arglag_res
sens_res$pooled$arglag$aic <- pooled_arglag_aic

# Site specific tests
city_arglag_res <- list()
city_arglag_aic <- list()
for (v in var){
  print(v)
  for (l in locations){
    sub <- df %>% filter(location==l)
    city_arglag_res[[v]][[l]] <- test_arglag(sub, outcome, pooled=FALSE, v)
    print(l)
  }
  city_arglag_aic[[v]] <- wrap_plots(map(city_arglag_res[[v]], 'aic_plot'), nrow=2) &
    plot_annotation(title=paste0('City-specific AIC: ', replace_exact_matches(v, nm),' lag | df'),
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
}
sens_res$city$arglag$res <- city_arglag_res
sens_res$city$arglag$aic <- city_arglag_aic

# Emergency ward tests
l <- 'Hamirpur'
sub <- df %>% filter(location==l)
emg_arglag_res <- list()
for (v in var){
  print(v)
  emg_arglag_res[[v]] <- test_arglag(sub, outcome='emergency', pooled=FALSE, v)
}
emg_arglag_aic <- wrap_plots(map(emg_arglag_res, 'aic_plot'), nrow=1) &
  plot_annotation(title='Emergency AIC: Variable lag | df',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
sens_res$emergency$arglag$res <- emg_arglag_res
sens_res$emergency$arglag$aic <- emg_arglag_aic


# EXPORT
sens_analysis <- list('3A_pooled_arglag_AIC'=list(p=sens_res$pooled$arglag$aic, w=8, h=4),
                      '3B_emg_arglag_AIC'=list(p=sens_res$emergency$arglag$aic, w=8, h=4),
                      '3C_city_arglag1_mx2t_AIC'=list(p=sens_res$city$arglag$aic$mx2t, w=12, h=8),
                      '3C_city_arglag2_tp_AIC'=list(p=sens_res$city$arglag$aic$tp, w=12, h=8),
                      '3C_city_arglag3_fg10_AIC'=list(p=sens_res$city$arglag$aic$fg10, w=12, h=8))

for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'figures/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# Check results (update best fit params based on outputs)
sens_res$pooled$arglag$aic
sens_res$emergency$arglag$aic
sens_res$city$arglag$aic

# Updates
updates <- list(
  'mx2t' = list('pooled' = 4, 'emergency'=3, 'Chamba' = 4, 'Dharamshala' = 5, 'Hamirpur' = 4, 'Kullu' = 4, 'Nahan' = 4, 'Reckong Peo' = 3, 'Solan' = 4, 'Una' = 3),
  'tp'   = list('pooled' = 4, 'emergency'=3, 'Chamba' = 5, 'Dharamshala' = 4, 'Hamirpur' = 4, 'Kullu' = 4, 'Nahan' = 3, 'Reckong Peo' = 3, 'Solan' = 4, 'Una' = 3),
  'fg10' = list('pooled' = 4, 'emergency'=3, 'Chamba' = 3, 'Dharamshala' = 4, 'Hamirpur' = 3, 'Kullu' = 3, 'Nahan' = 4, 'Reckong Peo' = 3, 'Solan' = 4, 'Una' = 5)
)
for (v in var){
  sens_res$pooled$arglag$res[[v]]$best_fit <- updates[[v]]$pooled
  sens_res$emergency$arglag$res[[v]]$best_fit <- updates[[v]]$emergency
  for(l in locations){
    sens_res$city$arglag$res[[v]][[l]]$best_fit <- updates[[v]][[l]]
  }
}

# 4. (ArgVar) Test single weather variables for shape structures ####

test_argvar <- function(df, outcome, pooled, var){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    city <- 'ALL'
    p <- 'p_raw'
    title <- replace_exact_matches(var, nm)
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
    arglag_df <- sens_res$pooled$arglag$res[[var]]$best_fit
  } else{
    city <- as.character(unique(df$location))
    p <- 'p_raw'    
    if(outcome == 'total_patients'){
      title <- city
      longterm_df <- sens_res$city$longterm$res[[city]]$best_fit
      autocorr_lag <- sens_res$city$autocorr$res[[city]]$best_fit
      arglag_df <- sens_res$city$arglag$res[[var]][[city]]$best_fit
    } else {
      title <- replace_exact_matches(var, nm)
      longterm_df <- sens_res$emergency$longterm$res$best_fit
      autocorr_lag <- sens_res$emergency$autocorr$res$best_fit
      arglag_df <- sens_res$emergency$arglag$res[[var]]$best_fit
    }
  }
  
  reference <- ifelse(var %in% c("tp", "fg10"), "0%", "50%")
  

  color <- as.character(c(mx2t='#D74B3C', tp='#1D98EB', fg10='#8968CD')[var] %||% 'black')
  
  var_shp <- list()
  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled, 
                                    longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                                    arg_lag=list(fun='ns', df=arglag_df),
                                    test='Variable Shape', variables=c(var), 
                                    arg_var=list(fun='lin'))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  var_shp[['Linear']] <- city_exp_res(city=city, lag=NULL, meta=test_result$meta_coef, col=color, xlab=paste0(replace_exact_matches(var, nm), ' - %'), title='Linear', reference=reference)[[p]]
  for (j in 3:6){ # test exp-resp dfs (between 3 and 8 as per lit)
    print(j)
    name <- paste0('NS df=', j)
    test_result <- get_parameter_test(data=df, locations=locations, outcome=outcome, pooled=pooled,
                                      longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                                      arg_lag=list(fun='ns', df=arglag_df),
                                      test='Variable Shape' , variables=c(var), 
                                      arg_var=list(fun='ns', df=j))
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
    var_shp[[name]] <- city_exp_res(city=city, lag=NULL, meta=test_result$meta_coef, col=color, xlab=paste0(replace_exact_matches(var, nm), ' - %'), title=name, reference=reference)[[p]]
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  dis_plot <- (wrap_plots(var_shp, nrow=2)& theme(legend.position='none')) +
      plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub('.*\\|(\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, dis_plot=dis_plot, best_fit=num))
}

# Pooled tests
pooled_argvar_res <- list()
for (v in var){
  print(v)
  pooled_argvar_res[[v]] <- test_argvar(df, outcome, pooled=TRUE, var=v)
}
pooled_argvar_aic <- wrap_plots(map(pooled_argvar_res, 'aic_plot'), nrow=1) &
    plot_annotation(title='Pooled AIC: Variable shape | df',
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
pooled_argvar_dis <- wrap_plots(map(pooled_argvar_res, 'dis_plot'), nrow=3) &
    plot_annotation(title='Pooled Relationship: Variable exposure distribution, by df',
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
sens_res$pooled$argvar$res <- pooled_argvar_res
sens_res$pooled$argvar$aic <- pooled_argvar_aic
sens_res$pooled$argvar$dis <- pooled_argvar_dis

# Site specific tests
city_argvar_res <- list()
city_argvar_aic <- list()
city_argvar_dis <- list()
for (v in var){
  print(v)
  for (l in locations){
    sub <- df %>% filter(location==l)
    city_argvar_res[[v]][[l]] <- test_argvar(df=sub, outcome, pooled=FALSE, var=v)
    print(l)
  }
  city_argvar_aic[[v]] <- wrap_plots(map(city_argvar_res[[v]], 'aic_plot'), nrow=2) &
    plot_annotation(title=paste0('City-specific AIC: ', replace_exact_matches(v, nm),' shape | df'),
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
}

# Get city specific distributions
for (l in locations){
  out <- map(city_argvar_res, l)
  city_argvar_dis[[l]] <- wrap_plots(map(out, 'dis_plot'), ncol=1) &
    plot_annotation(title=paste0('City-specific shapes, by df: ', l),
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
  
}
sens_res$city$argvar$res <- city_argvar_res
sens_res$city$argvar$aic <- city_argvar_aic
sens_res$city$argvar$dis <- city_argvar_dis

# Emergency ward tests
l <- 'Hamirpur'
sub <- df %>% filter(location==l)
emg_argvar_res <- list()
for (v in var){
  print(v)
  emg_argvar_res[[v]] <- test_argvar(sub, outcome='emergency', pooled=FALSE, v)
}
emg_argvar_aic <- wrap_plots(map(emg_argvar_res, 'aic_plot'), nrow=1) &
  plot_annotation(title='Emergency AIC: Variable shape | df',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
emg_argvar_dis <- wrap_plots(map(emg_argvar_res, 'dis_plot'), nrow=3) &
  plot_annotation(title='Emergency Relationship: Variable exposure distribution, by df',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
sens_res$emergency$argvar$res <- emg_argvar_res
sens_res$emergency$argvar$aic <- emg_argvar_aic
sens_res$emergency$argvar$dis <- emg_argvar_dis

# EXPORT
sens_analysis <- list('4A_pooled_argvar_AIC'=list(p=sens_res$pooled$argvar$aic, w=8, h=4),
                      '4A_pooled_argvar_DIS'=list(p=sens_res$pooled$argvar$dis, w=8, h=12),
                      '4B_emg_argvar_AIC'=list(p=sens_res$emergency$argvar$aic, w=8, h=4),
                      '4B_emg_argvar_DIS'=list(p=sens_res$emergency$argvar$dis, w=8, h=12),
                      '4C_city_argvar1_AIC_mx2t'=list(p=sens_res$city$argvar$aic$mx2t, w=12, h=8),
                      '4C_city_argvar2_AIC_tp'=list(p=sens_res$city$argvar$aic$tp, w=12, h=8),
                      '4C_city_argvar3_AIC_fg10'=list(p=sens_res$city$argvar$aic$fg10, w=12, h=8))

# AIC
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'figures/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}
# City distributions only
for (n in names(city_argvar_dis)) {
  title <- paste0('4C_city_argvar_DIS_',n)
  ggsave(filename=paste0(path_report, 'figures/sensitivity/SM_', title, '.png'), 
         plot=city_argvar_dis[[n]], w=8, h=12, dpi=300)
  print(n)
}

# Check results (update best fit params based on outputs)
sens_res$pooled$argvar$aic
sens_res$pooled$argvar$dis
sens_res$emergency$argvar$aic
sens_res$emergency$argvar$dis
sens_res$city$argvar$aic
sens_res$city$argvar$dis

# Updates
updates <- list(
  'mx2t' = list('pooled' = 3, 'emergency'=1, 'Chamba' = 3, 'Dharamshala' = 4, 'Hamirpur' = 3, 'Kullu' = 3, 'Nahan' = 3, 'Reckong Peo' = 3, 'Solan' = 3, 'Una' = 1),
  'tp'   = list('pooled' = 1, 'emergency'=1, 'Chamba' = 1, 'Dharamshala' = 3, 'Hamirpur' = 1, 'Kullu' = 1, 'Nahan' = 1, 'Reckong Peo' = 1, 'Solan' = 5, 'Una' = 3),
  'fg10' = list('pooled' = 1, 'emergency'=1, 'Chamba' = 1, 'Dharamshala' = 1, 'Hamirpur' = 1, 'Kullu' = 1, 'Nahan' = 1, 'Reckong Peo' = 1, 'Solan' = 1, 'Una' = 1)
)

for (v in var){
  sens_res$pooled$argvar$res[[v]]$best_fit <- updates[[v]]$pooled
  sens_res$emergency$argvar$res[[v]]$best_fit <- updates[[v]]$emergency
  for(l in locations){
    sens_res$city$argvar$res[[v]][[l]]$best_fit <- updates[[v]][[l]]
  }
}


# -- Store parameters -- ####

# Get pooled params
final_params <- list()
final_params$pooled$longterm <- sens_res$pooled$longterm$res$best_fit
final_params$pooled$autocorr <- sens_res$pooled$autocorr$res$best_fit
for(v in var){
  final_params$pooled$arglag[[v]] <- sens_res$pooled$arglag$res[[v]]$best_fit
  final_params$pooled$argvar[[v]] <- sens_res$pooled$argvar$res[[v]]$best_fit
}

# Get city params
for (l in locations){
  final_params[[l]]$longterm <- sens_res$city$longterm$res[[l]]$best_fit
  final_params[[l]]$autocorr <- sens_res$city$autocorr$res[[l]]$best_fit
  for(v in var){
    final_params[[l]]$arglag[[v]] <- sens_res$city$arglag$res[[v]][[l]]$best_fit
    final_params[[l]]$argvar[[v]] <- sens_res$city$argvar$res[[v]][[l]]$best_fit
  }
}

# Get emergency params
final_params$emergency$longterm <- sens_res$emergency$longterm$res$best_fit
final_params$emergency$autocorr <- sens_res$emergency$autocorr$res$best_fit
for(v in var){
  final_params$emergency$arglag[[v]] <- sens_res$emergency$arglag$res[[v]]$best_fit
  final_params$emergency$argvar[[v]] <- sens_res$emergency$argvar$res[[v]]$best_fit
}

# 5. (Var) Test added value of each weather variable within model ####

# Function for parameter testing
get_compvar_test <- function(data, locations, outcome, pooled=TRUE, variables=NULL, v_lag=14, 
                               longterm_df=8, autocorr_lag=14,
                               arg_lag=list('mx2t'=3,'tp'=1,'fg10'=3),
                               arg_var=list('mx2t'=3,'tp'=1,'fg10'=3)){
  # Build datalist
  datalist <- lapply(locations, function(region) data[data$location == region, ])
  names(datalist) <- locations
  
  # Drop empty datasets
  datalist <- Filter(function(x) nrow(x) > 0, datalist)
  
  # Use names that survived filtering
  locs_used <- names(datalist)
  
  m <- length(datalist)
  years <- length(unique(data$year))
  
  # Initialize output lists with correct names
  list_comp  <- vector('list', m)
  list_model <- vector('list', m)
  
  names(list_comp)  <- locs_used
  names(list_model) <- locs_used
  
  # PER-VARIABLE reference, ranges, centering, df, and storage objects

  # Reference percentile per variable
  reference <- setNames(
    lapply(variables, function(v) {
      if (v %in% c("tp", "fg10")) "10%" else "50%"
    }),
    variables
  )
  
  # Ranges per variable (locations × 2)
  ranges <- setNames(
    lapply(variables, function(v) {
      t(sapply(datalist, function(x) range(x[[v]], na.rm = TRUE)))
    }),
    variables
  )
  
  # Centering quantiles per variable (101 × locations)
  cen <- setNames(
    lapply(variables, function(v) {
      sapply(datalist, function(x) {
        quantile(x[[v]], c(0:100) / 100, na.rm = TRUE)
      })
    }),
    variables
  )
  
  arg_lag_list <- lapply(arg_lag, function(x) {
    if (x == 1) list(fun = "lin") else list(fun = "ns", df = x)
  })
  
  # Build arglag list per variable (always natural spline)
  arg_var_list <- lapply(arg_var, function(x) {
    list(fun = "ns", df = x)
  })
  
  # Degrees of freedom per variable (arg_var × arg_lag)
  df <- setNames(
    lapply(variables, function(v) {
      prod(arg_var[[v]], arg_lag[[v]])
    }),
    variables
  )
  
  # Storage matrices per variable
  ymat <- setNames(
    lapply(variables, function(v) {
      matrix(
        NA,
        m,
        df[[v]],
        dimnames = list(locs_used, paste0("spl", seq(df[[v]])))
      )
    }),
    variables
  )
  
  # Variance–covariance storage per variable
  Slist <- setNames(
    lapply(variables, function(v) vector("list", m)),
    variables
  )
  
  for (v in variables) {
    names(Slist[[v]]) <- locs_used
  }
  
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  for(i in seq(m)) {
  
    # PRINT ITERATION
    cat(i,'')
    
    # LOAD
    data <- datalist[[i]]
    
    # RUN THE MODEL & EXTRACT COMPARATIVE METRICS

      # Select lags
      if(autocorr_lag==0){
        lags <- ''
      } else {
        lag_terms <- paste0(outcome, '_lag', 1:autocorr_lag)
        lags <- paste0('+ ', paste(lag_terms, collapse=' + '))
      }
      
      # Fit cross basis
    nvar <- length(variables)
    
    if (nvar == 1) {
      
      cb1 <- crossbasis(data[[variables[1]]], c(0, v_lag), arg_var_list[[1]], arg_lag_list[[1]], data$year)
      
      data$cb1 <- cb1
      cb_terms <- "cb1"
      cb_main  <- cb1
      
    } else if (nvar == 2) {
      
      cb1 <- crossbasis(data[[variables[1]]], c(0, v_lag), arg_var_list[[1]], arg_lag_list[[1]], data$year)
      cb2 <- crossbasis(data[[variables[2]]], c(0, v_lag), arg_var_list[[2]], arg_lag_list[[2]], data$year)
      
      data$cb1 <- cb1
      data$cb2 <- cb2
      cb_terms <- "cb1 + cb2"
      cb_main  <- cb1   # anchor for metrics / prediction
      
    } else if (nvar == 3) {
      
      cb1 <- crossbasis(data[[variables[1]]], c(0, v_lag), arg_var_list[[1]], arg_lag_list[[1]], data$year)
      cb2 <- crossbasis(data[[variables[2]]], c(0, v_lag), arg_var_list[[2]], arg_lag_list[[2]], data$year)
      cb3 <- crossbasis(data[[variables[3]]], c(0, v_lag), arg_var_list[[3]], arg_lag_list[[3]], data$year)
      
      data$cb1 <- cb1; data$cb2 <- cb2; data$cb3 <- cb3
      cb_terms <- "cb1 + cb2 + cb3"
      cb_main  <- cb1
      
    } else if (nvar == 4) {
      
      cb1 <- crossbasis(data[[variables[1]]], c(0, v_lag), arg_var_list[[1]], arg_lag_list[[1]], data$year)
      cb2 <- crossbasis(data[[variables[2]]], c(0, v_lag), arg_var_list[[2]], arg_lag_list[[2]], data$year)
      cb3 <- crossbasis(data[[variables[3]]], c(0, v_lag), arg_var_list[[3]], arg_lag_list[[3]], data$year)
      cb4 <- crossbasis(data[[variables[4]]], c(0, v_lag), arg_var_list[[4]], arg_lag_list[[4]], data$year)
      
      data$cb1 <- cb1; data$cb2 <- cb2; data$cb3 <- cb3; data$cb4 <- cb4
      cb_terms <- "cb1 + cb2 + cb3 + cb4"
      cb_main  <- cb1
      
    } else if (nvar == 5) {
      
      cb1 <- crossbasis(data[[variables[1]]], c(0, v_lag), arg_var_list[[1]], arg_lag_list[[1]], data$year)
      cb2 <- crossbasis(data[[variables[2]]], c(0, v_lag), arg_var_list[[2]], arg_lag_list[[2]], data$year)
      cb3 <- crossbasis(data[[variables[3]]], c(0, v_lag), arg_var_list[[3]], arg_lag_list[[3]], data$year)
      cb4 <- crossbasis(data[[variables[4]]], c(0, v_lag), arg_var_list[[4]], arg_lag_list[[4]], data$year)
      cb5 <- crossbasis(data[[variables[5]]], c(0, v_lag), arg_var_list[[5]], arg_lag_list[[5]], data$year)
      
      data$cb1 <- cb1; data$cb2 <- cb2; data$cb3 <- cb3; data$cb4 <- cb4; data$cb5 <- cb5
      cb_terms <- "cb1 + cb2 + cb3 + cb4 + cb5"
      cb_main  <- cb1
      
    } else {
      stop("Function supports 1–5 variables only.")
    }
    
    baseformula <- as.formula(
      paste0(
        outcome,
        " ~ month + weekday + holiday +  ns(time, df=",
        longterm_df, "*", years, ") ",
        lags,
        " + ",
        cb_terms, '+ rh + pm25'
      )
    )
    
  
    model <- glm(formula=as.formula(baseformula),family=quasipoisson(),data)
    list_comp[[i]] <- comp_metrics(model)
    list_model[[i]] <- model
    
    # Make predictions
    if (nvar == 1) {
      
      pred1 <- crosspred(cb1, model, cen = cen[[variables[1]]][reference[[variables[1]]], i])
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      
    } else if (nvar == 2) {
      
      pred1 <- crosspred(cb1, model, cen = cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen = cen[[variables[2]]][reference[[variables[2]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      
    } else if (nvar == 3) {
      
      pred1 <- crosspred(cb1, model, cen = cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen = cen[[variables[2]]][reference[[variables[2]]], i])
      pred3 <- crosspred(cb3, model, cen = cen[[variables[3]]][reference[[variables[3]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      ymat[[variables[3]]][i, ]  <- pred3$coef
      Slist[[variables[3]]][[i]] <- pred3$vcov
      
    } else if (nvar == 4) {
      
      pred1 <- crosspred(cb1, model, cen = cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen = cen[[variables[2]]][reference[[variables[2]]], i])
      pred3 <- crosspred(cb3, model, cen = cen[[variables[3]]][reference[[variables[3]]], i])
      pred4 <- crosspred(cb4, model, cen = cen[[variables[4]]][reference[[variables[4]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      ymat[[variables[3]]][i, ]  <- pred3$coef
      Slist[[variables[3]]][[i]] <- pred3$vcov
      ymat[[variables[4]]][i, ]  <- pred4$coef
      Slist[[variables[4]]][[i]] <- pred4$vcov
      
    } else if (nvar == 5) {
      
      pred1 <- crosspred(cb1, model, cen = cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen = cen[[variables[2]]][reference[[variables[2]]], i])
      pred3 <- crosspred(cb3, model, cen = cen[[variables[3]]][reference[[variables[3]]], i])
      pred4 <- crosspred(cb4, model, cen = cen[[variables[4]]][reference[[variables[4]]], i])
      pred5 <- crosspred(cb5, model, cen = cen[[variables[5]]][reference[[variables[5]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      ymat[[variables[3]]][i, ]  <- pred3$coef
      Slist[[variables[3]]][[i]] <- pred3$vcov
      ymat[[variables[4]]][i, ]  <- pred4$coef
      Slist[[variables[4]]][[i]] <- pred4$vcov
      ymat[[variables[5]]][i, ]  <- pred5$coef
      Slist[[variables[5]]][[i]] <- pred5$vcov
      
    } else {
      stop("Function supports 1–5 variables only.")
    }
    
  }
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['sum_aic']] <- sum(sapply(list_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(list_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']] <- rowMeans(sapply(list_comp, function(x) as.numeric(x$acf)))
  meta_comp[['avg']] <- mean(sapply(list_comp, function(x) x$esti))

  meta_coef <- list()
  for (v in variables) {
    meta_coef[[v]] <- get_meta_analysis(pooled=pooled, ymat=ymat[[v]],  Slist=Slist[[v]], 
                                   arg_var=arg_var_list[[v]], arg_lag=arg_lag_list[[v]], ranges=ranges[[v]], cen=cen[[v]])
  }
  
  # Format AIC
  aic_result <- data.frame(paste(variables, collapse = "|"), meta_comp[['sum_aic']])
  names(aic_result) <- c('model','AIC' )
  
  return(list(aic_result=aic_result, models=list_model, meta_coef=meta_coef,
              pacf_result=meta_comp[['sum_pacf']], 
              acf_result=meta_comp[['sum_acf']]))
  
}

arglag <- final_params$pooled$arglag
argvar <- final_params$pooled$argvar

tmp <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('mx2t'),
                        arg_lag=list('mx2t'=arglag$mx2t), 
                        arg_var=list('mx2t'=argvar$mx2t))
tp <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('tp'),
                        arg_lag=list('tp'=arglag$tp), 
                        arg_var=list('tp'=argvar$tp))
fg10 <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('fg10'),
                        arg_lag=list('fg10'=arglag$fg10), 
                        arg_var=list('fg10'=argvar$fg10))
tmp_tp <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('mx2t', 'tp'),
                           arg_lag=list('mx2t'=arglag$mx2t,'tp'=arglag$tp),
                           arg_var=list('mx2t'=argvar$mx2t,'tp'=argvar$tp))
tmp_fg10 <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('mx2t', 'fg10'),
                             arg_lag=list('mx2t'=arglag$mx2t,'fg10'=arglag$fg10),
                             arg_var=list('mx2t'=argvar$mx2t,'fg10'=argvar$fg10))
tp_fg10 <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('tp', 'fg10'),
                             arg_lag=list('tp'=arglag$tp,'fg10'=arglag$fg10),
                             arg_var=list('tp'=argvar$tp,'fg10'=argvar$fg10))
tmp_tp_fg10 <- get_compvar_test(df, locations, outcome, pooled=TRUE, variables=c('mx2t', 'tp', 'fg10'),
                           arg_lag=list('mx2t'=arglag$mx2t,'tp'=arglag$tp, 'fg10'=arglag$fg10),
                           arg_var=list('mx2t'=argvar$mx2t,'tp'=argvar$tp, 'fg10'=argvar$fg10))

pooled_var_res <- rbind(tmp$aic_result, tp$aic_result, fg10$aic_result, tmp_tp$aic_result, 
                        tmp_fg10$aic_result, tp_fg10$aic_result, tmp_tp_fg10$aic_result)

# # Check tmp dis
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp$meta_coef$mx2t, col='red', xlab=NULL, title='tmp', reference='50%')[['p_raw']])
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp_tp$meta_coef$mx2t, col='red', xlab=NULL, title='tmp_tp', reference='50%')[['p_raw']])
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp_fg10$meta_coef$mx2t, col='red', xlab=NULL, title='tmp_fg10', reference='50%')[['p_raw']])
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp_tp_fg10$meta_coef$mx2t, col='red', xlab=NULL, title='tmp_tp_fg10', reference='50%')[['p_raw']])

# Check AIC: Suggests removing windspeed from analysis #278079.9
sens_res$pooled$var$aic <- pooled_var_res[order(pooled_var_res$AIC), ]



# 6. (Shape) Test shape of response #####

test_shape <- function(df, params, outcome, pooled=TRUE, city=NULL){
  
  # Get title
  title <- if (outcome == "total_patients") {
    if (pooled) "Pooled: " else paste0(unique(df$location), ": ")
  } else {
    paste0(replace_exact_matches(outcome, nm), ": ")
  }
  
  # Get params
  longterm_df <- params$longterm
  autocorr_lag <- params$autocorr
  arglag_tmp <- params$arglag$mx2t
  arglag_tp <- params$arglag$tp
  arglag_fg10 <- params$arglag$fg10
  argvar_tmp <- params$argvar$mx2t
  argvar_tp <- params$argvar$tp
  argvar_fg10 <- params$argvar$fg10
  
  # Run models
  models <- c('quasi', 'binom')
  test <- list()
  for(d in models){
    output <- run_model(df, locations, outcome, fam_type=d, pooled=pooled,
                            longterm_df=longterm_df, autocorr_lag=autocorr_lag, 
                            arglag_tmp=arglag_tmp, arglag_tp=arglag_tp, arglag_fg10=arglag_fg10,
                            argvar_tmp=argvar_tmp, argvar_tp=argvar_tp, argvar_fg10=argvar_fg10)
    test[[d]]$res <- format_results(output, outcome, city=city)
    test[[d]]$aic <- data.frame(model=d, AIC=output$meta_comp$sum_aic)
    test[[d]]$dis <-  test[[d]]$res$p_er_overall
    print(d)
  }
  
  # Checks
  modelfitstats <- rbind(test$quasi$aic, test$binom$aic)
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  
  # Plot dist
  qp_plot_titled <- test$quasi$dis + plot_annotation(title=paste0(title, 'Quasipoisson')) & theme(legend.position='none')
  nb_plot_titled <- test$binom$dis + plot_annotation(title=paste0(title, 'Negative binomial')) & theme(legend.position='none')
  dis_plot <- wrap_elements(full=qp_plot_titled) / wrap_elements(full=nb_plot_titled)
  
  return(list(best_fit=best_fit, dis_plot=dis_plot))
}

# Get pooled results
shape_res <- list()
shape_res$pooled <- test_shape(df=df, params=final_params$pooled, outcome=outcome, pooled=TRUE, city=NULL)
final_params$pooled$shape <- shape_res$pooled$best_fit
shape_res$pooled$dis_plot

# Get city-specific results
for (l in locations){
  sub <- df %>% filter(location==l)
  shape_res[[l]] <- test_shape(df=sub, params=final_params[[l]], outcome=outcome, pooled=FALSE, city=l)
  final_params[[l]]$shape <-  shape_res[[l]]$best_fit
  print(l)
}

# Emergency ward tests
l <- 'Hamirpur'
sub <- df %>% filter(location==l)
shape_res$emergency <- test_shape(df=sub, params=final_params$emergency, outcome='emergency', pooled=FALSE, city=l)
final_params$pooled$shape <- shape_res$emergency$best_fit

# EXPORT
for (n in names(shape_res)) {
  if(n=='pooled'){
    title <- paste0('5A_pooled_shape')
  } else if (n=='emergency'){
    title <- paste0('5B_emergency_shape')
  } else{
    title <- paste0('5C_city_shape_', n)
  }
  ggsave(filename=paste0(path_report, 'figures/sensitivity/SM_', title, '.png'),
         plot=shape_res[[n]]$dis_plot, width=8,height=6, dpi=300)
  print(n)
}

# EXPORT PARAMS ####

save(sens_res, final_params, file=paste0(path_project, 'ind_sens_params.RData'))

# End timer: 1 hours 31 minutes 10 seconds (20/01/2026) ####
end <- Sys.time()
dt <- round(as.numeric(end - start, units = 'secs'))
hours   <- dt %/% 3600
minutes <- (dt %% 3600) %/% 60
seconds <- dt %% 60
sprintf('%d hours %d minutes %.d seconds', hours, minutes, seconds)

