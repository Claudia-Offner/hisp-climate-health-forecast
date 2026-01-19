### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city Sensitivity  Analysis                                 ###
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()
setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# Load packages
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

#### FUNCTIONS ####

# Function for parameter testing
get_parameter_test <- function(data, locations, pooled,
                               test='Variable Lag', variables=NULL, v_lag=14, 
                               longterm_df=8, 
                               autocorr_lag=28,
                               arg_lag=list(fun='ns', df=3),
                               arg_var=list(fun='ns', df=4)){
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
  list_comp  <- vector("list", m)
  list_model <- vector("list", m)
  
  names(list_comp)  <- locs_used
  names(list_model) <- locs_used
  
  # IF variable=test parameters
  if(!is.null(variables)){
    
    if (variables %in% c("tp", "fg10")) {
      reference <- '10%'
    } else {
      reference <- '50%'
    }
    
    ranges <- t(sapply(datalist,function(x) range(x[[variables]],na.rm=T)))
    cen <- sapply(datalist,function(x) {
      quantile(x[[variables]],c(0:100)/100,na.rm=T)
    })
    df <- prod(arg_var$df, arg_lag$df)
    
    # Set storage matrix
    ymat <- matrix(NA,m,df,dimnames=list(locs_used,paste("spl",seq(df),sep="")))
    Slist <- vector("list",m)
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
        baseformula <- as.formula(paste0('total_patients ~  season + weekday +  ns(time, df=', longterm_df,'*', years, ')'))
        model <- glm(formula=baseformula,family=quasipoisson(),data)
        list_comp[[i]] <- comp_metrics(model)
        list_model[[i]] <- model
      } else if (test == 'Autocorrelation') {
        
        if (autocorr_lag==0) {
          # Fit model
          baseformula <- as.formula(paste0('total_patients ~  season + weekday + ns(time, df=', longterm_df,'*', years, ')'))
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model)
          list_model[[i]] <- model
        } else{
          
          # Select lags
          lag_terms <- paste0("total_patients_lag", 1:autocorr_lag)
          
          baseformula <- as.formula(
            paste0(
              "total_patients ~ season + weekday + ns(time, df=", longterm_df, "*", years, ") + ",
              paste(lag_terms, collapse=" + ")
            )
          )
          
          # Fit model
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model)
          list_model[[i]] <- model
        }
      } else if (grepl('Variable', test)){

        # Select lags
        lag_terms <- paste0("total_patients_lag", 1:autocorr_lag)
        
        # Fit cross basis
        cb <- crossbasis(data[[variables]], lag=c(0,v_lag), argvar=arg_var, arglag=arg_lag, group=data$year)
        # Fit model
        baseformula <- as.formula(
          paste0(
            "total_patients ~ season + weekday + ns(time, df=", longterm_df, "*", years, ") + ",
            paste(lag_terms, collapse=" + "), '+ cb'
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
         y='Pooled AIC') +
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

# Function to plot pooled ACF/PACF results
plot_pooled_acf <- function(test_result, plot='PACF'){
  
  # Conditionally set plot
  value_col <- if(plot=='ACF') test_result$acf_result else test_result$pacf_result
  if(plot=='ACF') value_col <- value_col[-1]   # Skip lag 0 for ACF
  
  # Set data frame
  df <- data.frame(Lag=1:length(value_col), Value=value_col)
  n_obs <- sum(sapply(test_result$models, function(m) length(m$residuals)))

  # Plot
  ggplot(df, aes(x=Lag, y=Value)) +
    geom_segment(aes(x=Lag, xend=Lag, y=0, yend=Value), linewidth=0.4) +
    geom_hline(yintercept=0, color='black') +
    geom_hline(yintercept=c(-1.96 / sqrt(n_obs), 1.96 / sqrt(n_obs)),
               linetype='dashed', color='blue') +
    labs(title=paste0(plot, ': ', test_result$aic_result$model), x='Lag', y=plot) +
    theme_classic() +
    theme(
      axis.line=element_line(color='black'),
      axis.ticks=element_line(color='black'),
      axis.text=element_text(color='black'),
      plot.title=element_text(size=14, face='bold', hjust=0.5)
    ) +
    coord_cartesian(ylim=c(-0.075, 0.1), xlim=c(0, 30)) +
    scale_x_continuous(breaks=round(seq(0, max(df$Lag), length.out=6)))
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# Store results
sens_res <- list()

# 1. (Longterm) Test time df for capturing long-term/seasonal ####

test_longterm <- function(df, pooled=TRUE){
  
  if(length(unique(df$location)) !=1 ){
    title <- 'Pooled AIC: Long-term Seasonality'
  } else{
    title <- unique(df$location)
  }
  
  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  for(i in 1:10){
    test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled,
                                      test='Longterm', longterm_df=i)
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub(".*\\((\\d+)\\).*", "\\1", best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, best_fit=num))
}

# Pooled tests
sens_res$pooled$longterm$res <- test_longterm(df, pooled=TRUE)
sens_res$pooled$longterm$aic <- sens_res$pooled$longterm$res$aic_plot

# Site specific tests
city_longterm_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- test_longterm(sub, pooled=FALSE)
  city_longterm_res[[l]] <- res
  print(l)
}
city_longterm_aic <- wrap_plots(map(city_longterm_res, "aic_plot"), nrow=2) &
    plot_annotation(title="City-specific AIC: Long-term seasonality",
    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
sens_res$city$longterm$res <- city_longterm_res
sens_res$city$longterm$aic <- city_longterm_aic


# # Check results (update best fit params based on outputs)
# sens_res$pooled$longterm$res
# sens_res$pooled$longterm$aic
# sens_res$city$longterm$res
# sens_res$city$longterm$aic

# EXPORT
sens_analysis <- list('pooled_1_longterm_AIC'=list(p=sens_res$pooled$longterm$aic, w=7, h=6),
                      'city_1_longterm_AIC'=list(p=sens_res$city$longterm$aic, w=12, h=8))
for (nm in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, '- figures/sensitivity/SM_', nm, '.png'), 
         plot=sens_analysis[[nm]]$p, width=sens_analysis[[nm]]$w,height=sens_analysis[[nm]]$h, dpi=300)
  print(nm)
}

# 2. (Autocorr) Test autocorrelation of hospital admissions #####

test_autocorr <- function(df, pooled=TRUE){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- 'Pooled AIC: Autocorrelation lags'
    longterm_df <- sens_res$pooled$longterm$res$best_fit
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
  }
  
  # Visualize ACF/PACF plot
  test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled, 
                                    longterm_df=longterm_df,
                                    test='Autocorrelation', autocorr_lag=0)
  acf_plot_raw <- plot_pooled_acf(test_result, plot='ACF')
  pacf_plot_raw <- plot_pooled_acf(test_result, plot='PACF')
  
  # Assess 1st AR lag term
  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  for(i in 0:7){
    test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled, 
                                      longterm_df=longterm_df,
                                      test='Autocorrelation', autocorr_lag=i)
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub(".*\\((\\d+)\\).*", "\\1", best_fit))
  num[is.na(num)] <- 1
  
  # Visualize final PACF plot
  test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled,
                                    longterm_df=longterm_df,
                                    test='Autocorrelation', autocorr_lag=num)
  acf_plot_fin <- plot_pooled_acf(test_result, plot='ACF')
  pacf_plot_fin <- plot_pooled_acf(test_result, plot='PACF')
  AIC_ACF <- acf_plot_raw + acf_plot_fin + plot_layout(ncol=2)
  AIC_PACF <- pacf_plot_raw + pacf_plot_fin + plot_layout(ncol=2)
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, best_fit=num, 
              AIC_ACF=AIC_ACF, AIC_PACF=AIC_PACF))
}

# Pooled tests
sens_res$pooled$autocorr$res <- test_autocorr(df, pooled=TRUE)
sens_res$pooled$autocorr$aic <- sens_res$pooled$autocorr$res$aic_plot

# Site specific tests
city_autocorr_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- test_autocorr(sub, pooled=FALSE)
  city_autocorr_res[[l]] <- res
  print(l)
}
city_autocorr_aic <- wrap_plots(map(city_autocorr_res, "aic_plot"), nrow=2) &
    plot_annotation(title="City-specific AIC: Autocorrelation lags",
                    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
sens_res$city$autocorr$res <- city_autocorr_res
sens_res$city$autocorr$aic <- city_autocorr_aic

# # Check results (update best fit params based on outputs)
# sens_res$pooled$autocorr$res
# sens_res$pooled$autocorr$aic
# sens_res$city$autocorr$res
# sens_res$city$autocorr$aic

# EXPORT
sens_analysis <- list('pooled_2_autocorr_AIC'=list(p=sens_res$pooled$autocorr$aic, w=7, h=6),
                      'city_2_autocorr_AIC'=list(p=sens_res$city$autocorr$aic, w=12, h=8))
for (nm in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, '- figures/sensitivity/SM_', nm, '.png'), 
         plot=sens_analysis[[nm]]$p, width=sens_analysis[[nm]]$w,height=sens_analysis[[nm]]$h, dpi=300)
  print(nm)
}

# 3. (Lag) Test single weather variables for lag structures ####

test_arglag <- function(df, var, pooled){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- replace_exact_matches(var, nm)
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
    autocorr_lag <- sens_res$city$autocorr$res[[title]]$best_fit
  }

  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  for(i in 3:4){
    test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled,
                                      longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                                      test='Variable Lag' , variables=c(var), 
                                      arg_lag=list(fun='ns', df=i),
                                      arg_var=list(fun='ns', df=4))
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub(".*\\|(\\d+)\\).*", "\\1", best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, best_fit=num))
}

var <- c('mx2t', 'tp', 'fg10')

# Pooled tests
pooled_arglag_res <- list()
for (v in var){
  pooled_arglag_res[[v]] <- test_arglag(df, v, pooled=TRUE)
}
pooled_arglag_aic <- wrap_plots(map(pooled_arglag_res, "aic_plot"), nrow=1) &
   plot_annotation(title="Pooled AIC: Variable lag df",
                   theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
sens_res$pooled$arglag$res <- pooled_arglag_res
sens_res$pooled$arglag$aic <- pooled_arglag_aic

# Site specific tests
city_arglag_res <- list()
city_arglag_aic <- list()
for (v in var){
  print(v)
  for (l in locations){
    sub <- df %>% filter(location==l)
    city_arglag_res[[v]][[l]] <- test_arglag(sub, v, pooled=FALSE)
    print(l)
  }
  city_arglag_aic[[v]] <- wrap_plots(map(city_arglag_res[[v]], "aic_plot"), nrow=2) &
    plot_annotation(title=paste0('City-specific AIC: ', replace_exact_matches(v, nm),' lag df'),
                    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
}
sens_res$city$arglag$res <- city_arglag_res
sens_res$city$arglag$aic <- city_arglag_aic

# # Check results (update best fit params based on outputs)
# sens_res$pooled$arglag$res
# sens_res$pooled$arglag$aic
# sens_res$city$arglag$res
# sens_res$city$arglag$aic

# EXPORT
sens_analysis <- list('pooled_3_arglag_AIC'=list(p=sens_res$pooled$arglag$aic, w=7, h=6),
                      'city_3_arglag_mx2t_AIC'=list(p=sens_res$city$arglag$aic$mx2t, w=12, h=8),
                      'city_3_arglag_tp_AIC'=list(p=sens_res$city$arglag$aic$tp, w=12, h=8),
                      'city_3_arglag_fg10_AIC'=list(p=sens_res$city$arglag$aic$fg10, w=12, h=8))

for (nm in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, '- figures/sensitivity/SM_', nm, '.png'), 
         plot=sens_analysis[[nm]]$p, width=sens_analysis[[nm]]$w,height=sens_analysis[[nm]]$h, dpi=300)
  print(nm)
}


# 4. (Shape) Test single weather variables for shape structures ####

test_argvar <- function(df, var, pooled){
  
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
    title <- paste0(city, ': Distribution (', var, ')')
    longterm_df <- sens_res$city$longterm$res[[city]]$best_fit
    autocorr_lag <- sens_res$city$autocorr$res[[city]]$best_fit
    arglag_df <- sens_res$city$arglag$res[[var]][[city]]$best_fit
  }

  color <- as.character(c(mx2t="#D74B3C", tp="#1D98EB", fg10="#8968CD")[var] %||% "black")
  
  var_shp <- list()
  setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
  test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled, 
                                    longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                                    arg_lag=list(fun='ns', df=arglag_df),
                                    test='Variable Shape', variables=c(var), 
                                    arg_var=list(fun='lin'))
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
  var_shp[['Linear']] <- city_exp_res(city=city, lag=NULL, meta=test_result$meta_coef, col=color, xlab=paste0(var, ' (%)'), title='Linear')[[p]]
  for (j in 3:4){ # test exp-resp dfs (between 3 and 8 as per lit)
    name <- paste0('NS df=', j)
    test_result <- get_parameter_test(data=df, locations=locations, pooled=pooled,
                                      longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                                      arg_lag=list(fun='ns', df=arglag_df),
                                      test='Variable Shape' , variables=c(var), 
                                      arg_var=list(fun='ns', df=j))
    modelfitstats <- rbind(modelfitstats, test_result$aic_result)
    var_shp[[name]] <- city_exp_res(city=city, lag=NULL, meta=test_result$meta_coef, col=color, xlab=paste0(var, ' (%)'), title=name)[[p]]
  }
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, best_fit, title)
  dis_plot <- (wrap_plots(var_shp, nrow=2)& theme(legend.position='none')) +
      plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub(".*\\|(\\d+)\\).*", "\\1", best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, aic_plot=aic_plot, dis_plot=dis_plot, best_fit=num))
}

var <- c('mx2t', 'tp', 'fg10')

# Pooled tests
pooled_argvar_res <- list()
for (v in var){
  pooled_argvar_res[[v]] <- test_argvar(df, var=v, pooled=TRUE)
}
pooled_argvar_aic <- wrap_plots(map(pooled_argvar_res, "aic_plot"), nrow=1) &
    plot_annotation(title="Pooled AIC: Variable exposure df",
                    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
pooled_argvar_dis <- wrap_plots(map(pooled_argvar_res, "dis_plot"), nrow=3) &
    plot_annotation(title="Pooled Relationship: Variable exposure distribution, by df",
                    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
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
    city_argvar_res[[v]][[l]] <- test_argvar(df=sub, var=v, pooled=FALSE)
    print(l)
  }
  city_argvar_aic[[v]] <- wrap_plots(map(city_argvar_res[[v]], "aic_plot"), nrow=2) &
    plot_annotation(title=paste0('City-specific AIC: ', replace_exact_matches(v, nm),' shape | df'),
                    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
}

# Get city specific distributions
for (l in locations){
  out <- map(city_argvar_res, l)
  city_argvar_dis[[l]] <- wrap_plots(map(out, 'dis_plot'), ncol=1) &
    plot_annotation(title=paste0('City-specific shapes, by df: ', l),
                    theme=theme(plot.title=element_text(hjust=0.5, face="bold")))
  
}
sens_res$city$argvar$res <- city_argvar_res
sens_res$city$argvar$aic <- city_argvar_aic
sens_res$city$argvar$dis <- city_argvar_dis

# # Check results (update best fit params based on outputs)
# sens_res$pooled$argvar$res
# sens_res$pooled$argvar$aic
# sens_res$pooled$argvar$dis
# sens_res$city$argvar$res
# sens_res$city$argvar$aic
# sens_res$city$argvar$dis

# EXPORT
sens_analysis <- list('pooled_4_argvar_AIC'=list(p=sens_res$pooled$argvar$aic, w=7, h=6),
                      'pooled_4_argvar_DIS'=list(p=sens_res$pooled$argvar$dis, w=15, h=15),
                      'city_4_argvar_AIC_mx2t'=list(p=sens_res$city$argvar$aic$mx2t, w=12, h=8),
                      'city_4_argvar_AIC_tp'=list(p=sens_res$city$argvar$aic$tp, w=12, h=8),
                      'city_4_argvar_AIC_fg10'=list(p=sens_res$city$argvar$aic$fg10, w=12, h=8))
# AIC
for (nm in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, '- figures/sensitivity/SM_', nm, '.png'), 
         plot=sens_analysis[[nm]]$p, width=sens_analysis[[nm]]$w,height=sens_analysis[[nm]]$h, dpi=300)
  print(nm)
}
# City dis only
for (nm in names(city_argvar_dis)) {
  title <- paste0('city_4_argvar_DIS_',nm)
  ggsave(filename=paste0(path_report, '- figures/sensitivity/SM_', title, '.png'), 
         plot=city_argvar_dis[[nm]], width=15,height=15, dpi=300)
  print(nm)
}


# EXPORT PARAMS ####

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

save(sens_res, final_params, file=paste0(path_project, 'sens_params.RData'))

# 5. (Shape) Test shape of response (Quasi) #####

load(paste0(path_project, '/sens_params.RData'))

test_shape <- function(df, params, pooled=TRUE, city=NULL){
  
  # Get title
  if (pooled==TRUE){
    title <- 'Pooled: ' 
  } else {
    title <- paste0(unique(df$location), ': ')
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
  outcome <- 'total_patients'
  models <- c('quasi', 'binom')
  test <- list()
  for(d in models){
    output <- run_model(df, locations, primary, fam_type=d, pooled=pooled,
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
  qp_plot_titled <- test$quasi$dis + plot_annotation(title=paste0(title, "Quasipoisson")) & theme(legend.position="none")
  nb_plot_titled <- test$binom$dis + plot_annotation(title=paste0(title, "Negative binomial")) & theme(legend.position="none")
  dis_plot <- wrap_elements(full=qp_plot_titled) / wrap_elements(full=nb_plot_titled)
  
  return(list(best_fit=best_fit, dis_plot=dis_plot))
}

# Get pooled results
shape_res <- list()
shape_res$pooled <- test_shape(df=df, params=final_params$pooled, pooled=TRUE, city=NULL)

# Get city-specific results
for (l in locations){
  sub <- df %>% filter(location==l)
  shape_res[[l]] <- test_shape(df=sub, params=final_params[[l]], pooled=FALSE, city=l)
}


# EXPORT
for (nm in names(shape_res)) {
  if(nm=='pooled'){
    title <- paste0('pooled_5_shape')
  } else{
    title <- paste0('city_5_shape_', nm)
  }
  ggsave(filename=paste0(path_report, '- figures/sensitivity/SM_', title, '.png'), 
         plot=city_argvar_dis[[nm]], width=12,height=10, dpi=300)
  print(nm)
}
