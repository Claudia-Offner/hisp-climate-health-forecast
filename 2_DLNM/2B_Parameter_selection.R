### ------------------------------------------------------------------------ ### 
### DLNM - Sensitivity Analysis (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

####  DEPENDENCIES  ####
source('2_DLNM/2_Utilities.R')

# Load packages
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork); library(TTR)

# Function for parameter testing
run_param <- function(data, locations, outcome, fam_type='quasi', pooled=TRUE, 
                    control_test=NULL, var_test=NULL, 
                    longterm_df=10, autocorr_lag=0, lag_window=c(0,10),
                    arglag_var=4, argvar_var=1){
  
  # 1. Data and dependency prep ----
  
  # LIST OF DATAFRAMES FOR ALL SITES
  datalist <- lapply(locations, function(region) data[data$location == region, ])
  names(datalist) <- locations
  
  # Drop empty datasets
  datalist <- Filter(function(x) nrow(x) > 0, datalist)
  
  # Use names that survived filtering
  locs_used <- names(datalist)
  m <- length(datalist)
  years <- length(unique(data$year))
  
  # Optional: Test controls
  if(!is.null(var_test)){
    ranges_var <- t(sapply(datalist,function(x) range(x[[var_test]],na.rm=T)))
    cen_var <- sapply(datalist,function(x) {
      quantile(x[[var_test]],c(0:100)/100,na.rm=T)
    })
    arglag_var <- list(fun='ns', df=arglag_var)
    argvar_var <- if (argvar_var == 1) list(fun="lin") else list(fun="ns", df=argvar_var)
    df_var <- prod(arglag_var$df, argvar_var$df)
    
  }

  # 2. Stage 1: Fit the DLNM ----
  
  # BUILT OBJECTS WHERE outputs WILL BE STORED:
  
  Slist_comp <- vector('list',m)
  Slist_model <- vector('list', m)
  names(Slist_comp) <- locs_used
  names(Slist_model) <- locs_used
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  system.time(
    for(i in seq(m)) {
      
      # PRINT ITERATION
      cat(i,'')
      
      # LOAD
      sub <- datalist[[i]]
      
      # SELECT AUTOCORR LAGS
      if(autocorr_lag==0){
        lags <- ''
      } else {
        lag_terms <- paste0(outcome, '_lag', 1:autocorr_lag)
        lags <- paste0('+ ', paste(lag_terms, collapse=' + '))
      }
      
      if(!is.null(var_test)){
        
        if (control_test=='None'){
          x_rep <- NULL
        } else if (control_test == 'Single Lag') {
          x_rep <- dplyr::lag(sub[[var_test]], n = 1)
          
        } else if (control_test == 'Moving Average') {
          x_rep <- TTR::SMA(na.approx(sub[[var_test]], na.rm = FALSE), n = lag_window[2])
          
        } else if (control_test == 'Exponential Moving Average') {
          x_rep <- TTR::EMA(na.approx(sub[[var_test]], na.rm = FALSE), n = lag_window[2], wilder = FALSE)
          
        } else if (control_test == 'Distributed Lag'){
          x_rep <- sub[[var_test]]
        }
        
        if(!is.null(x_rep)){
          
          # Store this location's range/percentiles for later pooling
          ranges_var[i, ] <- range(x_rep, na.rm = TRUE)
          cen_var[, i] <- quantile(x_rep, probs = c(0:100) / 100, na.rm = TRUE)
          
          # Build the basis (no lag dimension -- already collapsed above)
          if (control_test=='Distributed Lag'){
            
            cb <- crossbasis(x_rep, lag=lag_window, argvar=argvar_var, arglag=arglag_var)
            
          } else{
            
            if (argvar_var$fun == 'lin') {
              cb <- onebasis(x_rep, fun = 'lin')
            } else {
              cb <- onebasis(x_rep, fun = argvar_var$fun, df = argvar_var$df)
            }
          }
          
          
          controls <- '+ weekday + holiday + cb'
          
        } else {
          controls <- '+ weekday + holiday'
        }
        
      } else{
        controls <- ''
      }

      
      # Set formula
      form <- as.formula(
        paste0(outcome, 
               " ~ ns(time, df=", longterm_df, "*", years, ')', 
               controls, lags)
      )
      
      # RUN THE MODEL
      if(fam_type=='quasi'){
        model <- glm(formula=form, family=quasipoisson(), data=sub)
      } else if (fam_type=='binom'){
        model <- glm.nb(formula=form,data=sub)
      } else {
        stop('Unknown model type: ', fam_type)
      }
      
      # EXTRACT AND SAVE THE RELATED COEF AND VCOV
      
      # EXTRACT COMPARATIVE METRICS
      Slist_comp[[i]] <- comp_metrics(model)
      Slist_model[[i]] <- model
      
    })
  
  
  # 3. Stage 2: Multivariate Meta-analysis ----
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['sum_aic']] <- sum(sapply(Slist_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$acf)))
  
  # RESET WARNING
  options(warn=0)
  
  return(list(models=Slist_model, meta_comp=meta_comp))
  
}


# Function to plot AIC results (comparison)
plot_AIC_res <- function(modelfitstats, title, best_fit, xlab='Degrees of freedom'){
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, 'yes', 'no')
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels=modelfitstats$model)
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x=model, y=AIC, group=1)) +
    geom_line(color='grey20') +
    geom_point(aes(color=highlight), size=3) +
    scale_color_manual(values=c('no'='black', 'yes'='red')) +
    labs(title=title,
         y='Sum of AIC', x=xlab) +
    theme_classic() +
    theme(
      axis.line=element_line(color='gray20'),
      axis.text=element_text(color='gray20'),
      # axis.text.x=element_text(angle=45, hjust=1),
      # axis.title.x =element_blank(),
      axis.ticks=element_line(color='gray20'),
      plot.title=element_text(face='bold', hjust=0.5, size=12),
      legend.position='none'  # hide legend if not needed
    )
  return(p)
  
}

# Function to plot average estimates results (comparison)
plot_ESTI_res <- function(modelfitstats, exp){
  
  if(exp=='tmp'){
    y <- 'Average temperature estimates'
    modelfitstats$exp <- modelfitstats$avg_tmp
  } else {
    y <- 'Average precipitation estimates'
    modelfitstats$exp <- modelfitstats$avg_tp
  }
  
  best_fit <- if (mean(modelfitstats$AIC) < 0) {
    modelfitstats[which.max(modelfitstats$AIC), ]$model
  } else {
    modelfitstats[which.min(modelfitstats$AIC), ]$model
  }
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, 'yes', 'no')
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels=modelfitstats$model)
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x=model, y=exp, group=1)) +
    geom_line(color='grey20') +
    geom_point(aes(color=highlight), size=3) +
    scale_color_manual(values=c('no'='black', 'yes'='red')) +
    labs(title=NULL, y=y, x='Degree of freedom') +
    theme_classic() +
    theme(
      axis.line=element_line(color='gray20'),
      axis.text=element_text(color='gray20'),
      # axis.text.x=element_text(angle=45, hjust=1),
      axis.ticks=element_line(color='gray20'),
      plot.title=element_text(face='bold', hjust=0.5, size=12),
      legend.position='none'  # hide legend if not needed
    )
  return(p)
  
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData'))
df <- df_dlnm; rm(df_dlnm)

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# Time and store results
sens_res <- list()
start <- Sys.time()
outcome <- 'total_patients'

# 1. (Longterm) Test time df for capturing long-term/seasonal ####

test_longterm <- function(df, outcome, pooled=TRUE){
  
  if(length(unique(df$location)) !=1 ){
    title <-  'Pooled AIC & Estimates: Long-term seasonality'
  } else{
    title <- as.character(unique(df$location))
  }
  
  modelfitstats <- data.frame()
  for (i in 1:10) {
    test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                             pooled=pooled, longterm_df=i, autocorr_lag=0,
                             control_test='None')
    modelfitstats <- rbind(modelfitstats, data.frame(
      model  =i,
      AIC    =test_result$meta_comp$sum_aic
      # avg_tmp=test_result$meta_comp$avg_tmp,
      # avg_tp =test_result$meta_comp$avg_tp
    ))
  }
  
  # Get AIC
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, title, best_fit)
  # # Get avg tmp
  # esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
  # esti_tp_plot <- plot_ESTI_res(modelfitstats, exp='tp')
  # Combine plots
  # all_plot <- wrap_plots(list(aic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1)  &
  #   plot_annotation(title=title,theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
  
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, best_fit=num, aic_plot=aic_plot))
}

# Pooled tests
sens_res$pooled$longterm$res <- test_longterm(df, outcome, pooled=TRUE)
pool_longterm_aic <- sens_res$pooled$longterm$res$aic_plot

# Site specific tests
city_longterm_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- test_longterm(sub, outcome, pooled=FALSE)
  city_longterm_res[[l]] <- res
  print(l)
}
sens_res$city$longterm$res <- city_longterm_res
city_longterm_aic <- wrap_plots(map(city_longterm_res, 'aic_plot'), nrow=2) &
    plot_annotation(title='City-specific AIC & Estimates: Long-term seasonality',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))


# EXPORT
sens_analysis <- list('1A_pooled_longterm_AIC'=list(p=pool_longterm_aic, w=6, h=4),
                      '1B_city_longterm_AIC'=list(p=city_longterm_aic, w=15, h=8))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# # Updates
# sens_res$pooled$longterm$res$best_fit <- 8
# sens_res$emergency$longterm$res$best_fit <- 9
# updates <- list('Una'=10, 'Kangra'=8, 'Nahan'=8, 'Hamirpur'=8, 'Mandi'=8,
#                 'Chamba'=7, 'Dharamshala'=8, 'Solan'=6, 'Kullu'=8, 'Reckong Peo'=8)
# for(l in locations){
#   sens_res$city$longterm$res[[l]]$best_fit <- updates[[l]]
# }

# 2. (Autocorr) Test autocorrelation of hospital admissions #####

test_autocorr <- function(df, outcome, pooled=TRUE){
  
  # Function to plot pooled ACF/PACF results
  plot_pooled_acf <- function(test_result, plot='PACF', name=NULL, num=NULL){
    
    # Conditionally set plot
    value_col <- if(plot=='ACF') test_result$meta_comp$sum_acf else test_result$meta_comp$sum_pacf
    if(plot=='ACF') value_col <- value_col[-1]   # Skip lag 0 for ACF
    
    # Set data frame
    df <- data.frame(Lag=1:length(value_col), Value=value_col)
    n_obs <- sum(sapply(test_result$models, function(m) length(m$residuals)))
    
    # Title
    num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', num))
    
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
      coord_cartesian(ylim=c(-0.2, 0.2), xlim=c(0, 30)) +
      scale_x_continuous(breaks=round(seq(0, max(df$Lag), length.out=6)))
  }
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- 'Pooled AIC'
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    name <- 'Pooled'
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
    name <- title
  }
  
  # Visualize ACF/PACF plot
  test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                         pooled=pooled, longterm_df=longterm_df, autocorr_lag=0)
  
  acf_plot_raw <- plot_pooled_acf(test_result, plot='ACF', name=name, num=0)
  pacf_plot_raw <- plot_pooled_acf(test_result, plot='PACF', name=name, num=0)
  
  # Assess 1st AR lag term
  modelfitstats <- data.frame()
  for (i in 1:14) {
    test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                             pooled=pooled, longterm_df=longterm_df, autocorr_lag=i)
    
    modelfitstats <- rbind(modelfitstats, data.frame(
      model  =i,
      AIC    =test_result$meta_comp$sum_aic
      # avg_tmp=test_result$meta_comp$avg_tmp,
      # avg_tp =test_result$meta_comp$avg_tp
    ))
  }
  
  # Get AIC
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, title, best_fit, 'Autocorrelation Lag')
  # Get avg tmp
  # esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
  # esti_tp_plot <- plot_ESTI_res(modelfitstats, exp='tp')
  # # Combine plots
  # all_plot <- wrap_plots(list(aic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1)  &
  #   plot_annotation(title=title,theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
  
  modelfitstats <- modelfitstats %>% arrange(AIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  # Visualize final PACF plot
  test_result <- run_param(data=df, locations=locations, outcome=outcome, pooled=pooled, autocorr_lag=num)
  acf_plot_fin <- plot_pooled_acf(test_result, plot='ACF', name=name, num)
  pacf_plot_fin <- plot_pooled_acf(test_result, plot='PACF', name=name, num)
  
  AIC_ACF <- list(); AIC_PACF <- list()
  AIC_ACF$raw <- acf_plot_raw
  AIC_ACF$fin <- acf_plot_fin
  AIC_PACF$raw <- pacf_plot_raw
  AIC_PACF$fin <- pacf_plot_fin
  AIC_ACF$all <- acf_plot_raw + acf_plot_fin + plot_layout(ncol=2)
  AIC_PACF$all <- pacf_plot_raw + pacf_plot_fin + plot_layout(ncol=2)
  
  return(list(modelfitstats=modelfitstats,  best_fit=num, aic_plot=aic_plot,
              AIC_ACF=AIC_ACF, AIC_PACF=AIC_PACF))
}

# Pooled tests
sens_res$pooled$autocorr$res <- test_autocorr(df, outcome, pooled=TRUE)

# Site specific tests
city_autocorr_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- test_autocorr(sub, outcome, pooled=FALSE)
  city_autocorr_res[[l]] <- res
  print(l)
}
sens_res$city$autocorr$res <- city_autocorr_res
city_autocorr_aic <- wrap_plots(map(city_autocorr_res, 'aic_plot'), nrow=2) &
  plot_annotation(title='City-specific AIC: Autocorrelation lags',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
city_autocorr_acf <- (wrap_plots(map(map(city_autocorr_res, "AIC_ACF"), 'all'), nrow=4) & theme(axis.title.x=element_blank()))  +
  plot_annotation(
    title="City-specific ACF Comparison",
    theme=theme(plot.title=element_text(hjust=0.5, face="bold"))
  )

# ALL plots for appendix
all <- list(sens_res$pooled$autocorr$res$aic_plot, sens_res$pooled$autocorr$res$AIC_ACF$raw, sens_res$pooled$autocorr$res$AIC_ACF$fin)
pooled_autocorr <- (wrap_plots(all, nrow=1))  +
  plot_annotation(
    tag_levels = 'A',
    theme=theme(plot.title=element_text(hjust=0.5, face="bold"))
  )


# EXPORT
sens_analysis <- list('2A_pooled_autocorr_ALL'=list(p=pooled_autocorr, w=12, h=4),
                      '2A_pooled_autocorr_AIC'=list(p=sens_res$pooled$autocorr$res$aic_plot, w=6, h=4),
                      '2A_pooled_autocorr_ACF'=list(p=sens_res$pooled$autocorr$res$AIC_ACF$all, w=6, h=4),
                      '2B_city_autocorr_AIC'=list(p=city_autocorr_aic, w=12, h=8),
                      '2B_city_autocorr_ACF'=list(p=city_autocorr_acf, w=12, h=12))

for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# # Updates
# sens_res$pooled$autocorr$res$best_fit <- 14
# sens_res$emergency$autocorr$res$best_fit <- 14
# updates <- list('Una'=14, 'Kangra'=14, 'Nahan'=14, 'Hamirpur'=14, 'Mandi'=10,
#                 'Chamba'=14, 'Dharamshala'=1, 'Solan'=11, 'Kullu'=14, 'Reckong Peo'=10)
# for(l in locations){
#   sens_res$city$autocorr$res[[l]]$best_fit <- updates[[l]]
# }

# 3. (Var) Test added value of each weather variable within model ####

test_var <- function (df, outcome, exp, pooled=TRUE){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <-  paste0('Pooled AIC & Estimates: Time-varying confounders (', 
                     replace_exact_matches(exp, nm), ')')
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
    autocorr_lag <- sens_res$city$longterm$res[[title]]$best_fit
  }

  lag_df <- 4; var_df <- 3
  settings <- list(
    C0=list(test='None'),
    C1=list(test='Single Lag',                 arg_var=1),
    C2=list(test='Moving Average',             arg_var=1),
    C3=list(test='Exponential Moving Average', arg_var=1),
    C4=list(test='Distributed Lag',            arg_var=1, arg_lag=lag_df),
    C5=list(test='Single Lag',                 arg_var=var_df),
    C6=list(test='Moving Average',             arg_var=var_df),
    C7=list(test='Exponential Moving Average', arg_var=var_df),
    C8=list(test='Distributed Lag',            arg_var=var_df, arg_lag=lag_df)
  )
  
  modelfitstats <- data.frame()
  for (cname in names(settings)) {
    
    s <- settings[[cname]]
    cat('\n[', cname, '-', s$test, ']\n')
    
    if (cname == 'C0') {

      test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                               pooled=pooled, longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                               control_test=s$test, var_test=exp)
      
    } else {
      test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                               pooled=pooled, longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                               control_test=s$test, var_test=exp, 
                               argvar_var=s$arg_var, arglag_var=s$arg_lag)
    }
    
    modelfitstats <- rbind(modelfitstats, data.frame(
      model  =cname,
      AIC    =test_result$meta_comp$sum_aic
      # avg_tmp=test_result$meta_comp$avg_tmp,
      # avg_tp =test_result$meta_comp$avg_tp
    )) 
    
  }

  # Get AIC
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  aic_plot <- plot_AIC_res(modelfitstats, title, best_fit, xlab='Model')
  # Get avg tmp
  # esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
  # esti_tp_plot <- plot_ESTI_res(modelfitstats, exp='tp')
  # Combine plots
  # all_plot <- wrap_plots(list(aic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1)  &
  #   plot_annotation(title=title,theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

  modelfitstats <- modelfitstats %>% arrange(AIC)

  return(list(modelfitstats=modelfitstats, best_fit=best_fit, aic_plot=aic_plot))

}

controls <- c('mx2t', 'tp', 'pm25', 'fg10', 'rh')
results <- list()
for (c in controls){
  
  print(c)
  results[[c]] <- test_var(df, outcome, exp=c, pooled=TRUE)
  
}

sens_res$pooled$confound$res <- results

plots <- map(results, 'aic_plot')
titles <- replace_exact_matches(controls, nm)
pool_confound_aic <- wrap_plots(map2(plots, titles, \(p, t) p + ggtitle(t)), nrow=2) &
  plot_annotation(title='Pooled AIC: Time-varying exposures',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

# EXPORT
sens_analysis <- list('3A_pooled_confound_ALL'=list(p=pool_confound_aic, w=15, h=5))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}


# 3. (Var) Check added confounder combinations #####

# Function for parameter testing
get_compvar_test <- function(data, locations, outcome, pooled=TRUE, variables=NULL, v_lag=10,
                             longterm_df=8, autocorr_lag=14, arg_lag=4, arg_var=3){
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
      if (v %in% c('mx2t', 'rh')) "50%" else "0%"
    }),
    variables
  )
  
  # Ranges per variable (locations × 2)
  ranges <- setNames(
    lapply(variables, function(v) {
      t(sapply(datalist, function(x) range(x[[v]], na.rm=TRUE)))
    }),
    variables
  )
  
  # Centering quantiles per variable (101 × locations)
  cen <- setNames(
    lapply(variables, function(v) {
      sapply(datalist, function(x) {
        quantile(x[[v]], c(0:100) / 100, na.rm=TRUE)
      })
    }),
    variables
  )
  
  arg_lag_list <- lapply(arg_lag, function(x) {
    if (x == 1) list(fun="lin") else list(fun="ns", df=x)
  })
  
  # Build arglag list per variable (always natural spline)
  arg_var_list <- lapply(arg_var, function(x) {
    if (x == 1) list(fun="lin") else list(fun="ns", df=x)
  })
  
  # Degrees of freedom per variable (arg_var × arg_lag)
  df <- setNames(
    lapply(variables, function(v) {
      var_df  <- if (arg_var[[v]] == 1) 1 else arg_var[[v]]
      lag_df  <- if (arg_lag[[v]] == 1) 1 else arg_lag[[v]]
      if (var_df * lag_df == 1) 2 else var_df * lag_df
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
        dimnames=list(locs_used, paste0("spl", seq(df[[v]])))
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
        cb_terms
      )
    )
    
    
    model <- glm(formula=as.formula(baseformula),family=quasipoisson(),data)
    list_comp[[i]] <- comp_metrics(model)
    list_model[[i]] <- model
    
    # Make predictions
    if (nvar == 1) {
      
      pred1 <- crosspred(cb1, model, cen=cen[[variables[1]]][reference[[variables[1]]], i])
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      
    } else if (nvar == 2) {
      
      pred1 <- crosspred(cb1, model, cen=cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen=cen[[variables[2]]][reference[[variables[2]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      
    } else if (nvar == 3) {
      
      pred1 <- crosspred(cb1, model, cen=cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen=cen[[variables[2]]][reference[[variables[2]]], i])
      pred3 <- crosspred(cb3, model, cen=cen[[variables[3]]][reference[[variables[3]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      ymat[[variables[3]]][i, ]  <- pred3$coef
      Slist[[variables[3]]][[i]] <- pred3$vcov
      
    } else if (nvar == 4) {
      
      pred1 <- crosspred(cb1, model, cen=cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen=cen[[variables[2]]][reference[[variables[2]]], i])
      pred3 <- crosspred(cb3, model, cen=cen[[variables[3]]][reference[[variables[3]]], i])
      pred4 <- crosspred(cb4, model, cen=cen[[variables[4]]][reference[[variables[4]]], i])
      
      ymat[[variables[1]]][i, ]  <- pred1$coef
      Slist[[variables[1]]][[i]] <- pred1$vcov
      ymat[[variables[2]]][i, ]  <- pred2$coef
      Slist[[variables[2]]][[i]] <- pred2$vcov
      ymat[[variables[3]]][i, ]  <- pred3$coef
      Slist[[variables[3]]][[i]] <- pred3$vcov
      ymat[[variables[4]]][i, ]  <- pred4$coef
      Slist[[variables[4]]][[i]] <- pred4$vcov
      
    } else if (nvar == 5) {
      
      pred1 <- crosspred(cb1, model, cen=cen[[variables[1]]][reference[[variables[1]]], i])
      pred2 <- crosspred(cb2, model, cen=cen[[variables[2]]][reference[[variables[2]]], i])
      pred3 <- crosspred(cb3, model, cen=cen[[variables[3]]][reference[[variables[3]]], i])
      pred4 <- crosspred(cb4, model, cen=cen[[variables[4]]][reference[[variables[4]]], i])
      pred5 <- crosspred(cb5, model, cen=cen[[variables[5]]][reference[[variables[5]]], i])
      
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
  aic_result <- data.frame(paste(variables, collapse="|"), meta_comp[['sum_aic']])
  names(aic_result) <- c('model','AIC' )
  
  return(list(aic_result=aic_result, models=list_model, meta_coef=meta_coef,
              pacf_result=meta_comp[['sum_pacf']],
              acf_result=meta_comp[['sum_acf']]))
  
}


# # Pre-compute EMA columns for any _ema variables
# df <- df %>%
#   group_by(location) %>%
#   mutate(across(all_of(c('pm25', 'rh', 'fg10')), 
#                 ~ TTR::EMA(zoo::na.approx(.x, na.rm=FALSE), n=14, wilder=FALSE),
#                 .names = "{.col}_ema")) %>%
#   mutate(across(all_of(c('rh', 'fg10')), 
#                 ~ dplyr::lag(.x, n=1),
#                 .names = "{.col}_lag1")) %>%
#   ungroup()
# 
# 
# variables_ema<- c('mx2t', 'tp', 'rh_ema', 'fg10_ema', 'pm25_ema')
# lin_arglag_ema <- list('mx2t'=4, 'tp'=4, 'rh_ema'=1, 'fg10_ema'=1, 'pm25_ema'=1)
# lin_argvar_ema <- list('mx2t'=3, 'tp'=1, 'rh_ema'=1, 'fg10_ema'=1, 'pm25_ema'=1)
# nonlin_arglag_ema <- list('mx2t'=4, 'tp'=4, 'rh_ema'=4, 'fg10_ema'=4, 'pm25_ema'=1)
# nonlin_argvar_ema <- list('mx2t'=3, 'tp'=1, 'rh_ema'=3, 'fg10_ema'=1, 'pm25_ema'=1)
# 
# variables_lag <- c('mx2t', 'tp', 'rh_lag1', 'fg10_lag1', 'pm25_ema')
# lin_arglag_lag <- list('mx2t'=4, 'tp'=4, 'rh_lag1'=1, 'fg10_lag1'=1, 'pm25_ema'=1)
# lin_argvar_lag <- list('mx2t'=3, 'tp'=1, 'rh_lag1'=1, 'fg10_lag1'=1, 'pm25_ema'=1)
# nonlin_arglag_lag <- list('mx2t'=4, 'tp'=4, 'rh_lag1'=4, 'fg10_lag1'=4, 'pm25_ema'=1)
# nonlin_argvar_lag <- list('mx2t'=3, 'tp'=1, 'rh_lag1'=3, 'fg10_lag1'=1, 'pm25_ema'=1)
# 
# variables_all <- c('mx2t', 'tp', 'rh', 'fg10', 'pm25_ema')
# lin_arglag <- list('mx2t'=4, 'tp'=4, 'rh'=1, 'fg10'=1, 'pm25_ema'=1)
# lin_argvar <- list('mx2t'=3, 'tp'=1, 'rh'=1, 'fg10'=1, 'pm25_ema'=1)
# nonlin_arglag <- list('mx2t'=4, 'tp'=4, 'rh'=4, 'fg10'=4, 'pm25_ema'=1)
# nonlin_argvar <- list('mx2t'=3, 'tp'=1, 'rh'=3, 'fg10'=1, 'pm25_ema'=1)
# 
# 
# scenarios <- list(
#   linear     = list(vars=variables_all, arglag=lin_arglag,    argvar=lin_argvar),
#   nonlinear  = list(vars=variables_all, arglag=nonlin_arglag, argvar=nonlin_argvar),
#   linear_lag = list(vars=variables_lag, arglag=lin_arglag_lag, argvar=lin_argvar_lag),
#   nonlin_lag = list(vars=variables_lag, arglag=nonlin_arglag_lag, argvar=nonlin_argvar_lag),
#   linear_ema = list(vars=variables_ema, arglag=lin_arglag_ema, argvar=lin_argvar_ema),
#   nonlin_ema = list(vars=variables_ema, arglag=nonlin_arglag_ema, argvar=nonlin_argvar_ema)
#   )
# 
# all_results <- lapply(scenarios, function(s) {
#   combos <- unlist(lapply(seq_along(s$vars), function(k) combn(s$vars, k, simplify=FALSE)), recursive=FALSE)
#   
#   res <- lapply(combos, function(vars) {
#     print(vars)
#     get_compvar_test(df, locations, outcome, pooled=TRUE,
#                      variables=vars, arg_lag=s$arglag[vars], arg_var=s$argvar[vars])
#   })
#   names(res) <- sapply(combos, paste, collapse="|")
#   
#   out <- do.call(rbind, lapply(res, `[[`, "aic_result"))
#   row.names(out) <- NULL
#   out
# })
# 
# pooled_var_res <- do.call(rbind, Map(function(d, nm) { d$type <- nm; d }, all_results, names(all_results)))
# pooled_var_res <- pooled_var_res[, c('type', setdiff(names(pooled_var_res), 'type'))]
# row.names(pooled_var_res) <- NULL
# pooled_var_res[order(pooled_var_res$AIC), ]


# # Check tmp dis
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp$meta_coef$mx2t, col='red', xlab=NULL, title='tmp', reference='50%')[['p_raw']])
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp_tp$meta_coef$mx2t, col='red', xlab=NULL, title='tmp_tp', reference='50%')[['p_raw']])
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp_fg10$meta_coef$mx2t, col='red', xlab=NULL, title='tmp_fg10', reference='50%')[['p_raw']])
# (x <- city_exp_res(city='ALL', lag=NULL, meta=tmp_tp_fg10$meta_coef$mx2t, col='red', xlab=NULL, title='tmp_tp_fg10', reference='50%')[['p_raw']])

# Check AIC: Suggests removing windspeed from analysis #278079.9
# sens_res$pooled$var$aic <- pooled_var_res[order(pooled_var_res$AIC), ]




# 4. (ArgLag) Test single weather variables for lag structures ####

test_arglag <- function(df, outcome, pooled){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- 'Pooled AIC: Variable lag'
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
    autocorr_lag <- sens_res$city$autocorr$res[[title]]$best_fit
  }
  
  # Get tmp and tp lag params
  fit_arglag <- function(param, ex) {
    modelfitstats <- do.call(rbind, lapply(3:6, \(i) {
      args <- list(data=df, locations=locations, outcome=outcome,
                   pooled=pooled, longterm_df=longterm_df, autocorr_lag=autocorr_lag)
      args[[param]] <- i
      test_result <- do.call(run_model, args)
      data.frame(model=i, AIC=test_result$meta_comp$sum_aic,
                 avg_tmp=test_result$meta_comp$avg_tmp, avg_tp=test_result$meta_comp$avg_tp)
    }))
    
    best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
    aic_plot <- plot_AIC_res(modelfitstats, ex, best_fit)
    # esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
    # esti_tp_plot  <- plot_ESTI_res(modelfitstats, exp='tp')
    # all_plot <- wrap_plots(list(aic_plot + theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1) &
    #   plot_annotation(title=paste0(title, '(', ex, ')'), theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
    
    modelfitstats <- modelfitstats %>% arrange(AIC)
    num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
    num[is.na(num)] <- 1
    
    return(list(modelfitstats=modelfitstats, best_fit=num, aic_plot=aic_plot))
  }
  
  modelfitstats_tmp <- fit_arglag('arglag_tmp', replace_exact_matches('mx2t',nm))
  modelfitstats_tp  <- fit_arglag('arglag_tp',  replace_exact_matches('tp',nm))
  modelfitstats_fg10  <- fit_arglag('arglag_fg10',  replace_exact_matches('fg10',nm))
  
  return(list(tmp=modelfitstats_tmp, tp=modelfitstats_tp, fg10=modelfitstats_fg10))
}

# Pooled tests
pooled_arglag_res <- test_arglag(df, outcome, pooled=TRUE)
sens_res$pooled$arglag$res <- pooled_arglag_res
pooled_arglag_aic <- wrap_plots(map(pooled_arglag_res, 'aic_plot'), nrow=1) &
   plot_annotation(title='Pooled AIC: Variable lag',
                   theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

# Site specific tests
city_arglag_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  city_arglag_res[[l]] <- test_arglag(sub, outcome, pooled=FALSE)
  print(l)
}
sens_res$city$arglag$res <- city_arglag_res
city_arglag_aic_tmp <- wrap_plots(
  imap(city_arglag_res, \(city, locations) city$tmp$aic_plot + ggtitle(locations)),
  nrow = 2
) &
  plot_annotation(title = paste0('City-specific AIC: ', replace_exact_matches('mx2t', nm), ' lag'),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold')))
city_arglag_aic_tp <- wrap_plots(
  imap(city_arglag_res, \(city, locations) city$tp$aic_plot + ggtitle(locations)),
  nrow = 2
) &
  plot_annotation(title = paste0('City-specific AIC: ', replace_exact_matches('tp', nm), ' lag'),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold')))

# EXPORT
sens_analysis <- list('4A_pooled_arglag_AIC'=list(p=pooled_arglag_aic, w=10, h=4),
                      '4B_city_arglag_mx2t_AIC'=list(p=city_arglag_aic_tmp, w=12, h=8),
                      '4B_city_arglag_tp_AIC'=list(p=city_arglag_aic_tp, w=12, h=8))

for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}


# # Updates
# updates <- list(
#   'mx2t'=list('pooled'=4, 'emergency'=3, 'Una'=3, 'Kangra'=6, 'Nahan'=4, 'Hamirpur'=4, 'Mandi'=4, 'Chamba'=4, 'Dharamshala'=6, 'Solan'=4, 'Kullu'=4, 'Reckong Peo'=3),
#   'tp'  =list('pooled'=4, 'emergency'=3, 'Una'=3, 'Kangra'=3, 'Nahan'=4, 'Hamirpur'=4, 'Mandi'=4, 'Chamba'=4, 'Dharamshala'=4, 'Solan'=4, 'Kullu'=4, 'Reckong Peo'=3),
#   'fg10'=list('pooled'=4, 'emergency'=3, 'Una'=3, 'Kangra'=3, 'Nahan'=4, 'Hamirpur'=3, 'Mandi'=4, 'Chamba'=4, 'Dharamshala'=6, 'Solan'=4, 'Kullu'=4, 'Reckong Peo'=3)
# )
# 
# for (v in var){
#   sens_res$pooled$arglag$res[[v]]$best_fit <- updates[[v]]$pooled
#   sens_res$emergency$arglag$res[[v]]$best_fit <- updates[[v]]$emergency
#   for(l in locations){
#     sens_res$city$arglag$res[[v]][[l]]$best_fit <- updates[[v]][[l]]
#   }
# }

# 5. (ArgVar) Test single weather variables for shape structures ####

test_argvar <- function(df, outcome, pooled){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    city <- 'ALL'
    p <- 'p_raw'
    title <- NULL #replace_exact_matches(var, nm)
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
    arglag_tmp_df <- sens_res$pooled$arglag$res$tmp$best_fit
    arglag_tp_df <- sens_res$pooled$arglag$res$tp$best_fit
  } else{
    city <- as.character(unique(df$location))
    p <- 'p_raw'    
    title <- city
    longterm_df <- sens_res$city$longterm$res[[city]]$best_fit
    autocorr_lag <- sens_res$city$autocorr$res[[city]]$best_fit
    arglag_tmp_df <- sens_res$city$arglag$res[[city]]$tmp$best_fit
    arglag_tp_df <- sens_res$city$arglag$res[[city]]$tp$best_fit
  }
  
  fit_argvar <- function(param, meta_key, ref, col, title) {
    
    var_shp <- list()
    modelfitstats <- data.frame()
    
    # Linear
    args <- list(data=df, locations=locations, outcome=outcome,
                 pooled=pooled, longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                 arglag_tmp=arglag_tmp_df, arglag_tp=arglag_tp_df)
    args[[param]] <- 1
    test_result <- do.call(run_model, args)
    modelfitstats <- rbind(modelfitstats, data.frame(
      model=  'Linear',
      AIC    = test_result$meta_comp$sum_aic,
      avg_tmp= test_result$meta_comp$avg_tmp,
      avg_tp = test_result$meta_comp$avg_tp
    ))
    var_shp[['Linear']] <- city_exp_res(city=city, lag=NULL, meta=test_result[[meta_key]], 
                                        col=col, xlab='Percentiles (%)', title='Linear', reference=ref)[[p]]
    
    # NS df=3:6
    for (i in 3:6) {
      name <- paste0('NS df=', i)
      args[[param]] <- i
      test_result <- do.call(run_model, args)
      modelfitstats <- rbind(modelfitstats, data.frame(
        model  = name,
        AIC    = test_result$meta_comp$sum_aic,
        avg_tmp= test_result$meta_comp$avg_tmp,
        avg_tp = test_result$meta_comp$avg_tp
      ))
      var_shp[[name]] <- city_exp_res(city=city, lag=NULL, meta=test_result[[meta_key]], 
                                      col=col, xlab='Percentiles (%)', title=name, reference=ref)[[p]]
    }
    
    best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
    aic_plot <- plot_AIC_res(modelfitstats, title, best_fit)
    esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
    esti_tp_plot  <- plot_ESTI_res(modelfitstats, exp='tp')
    all_plot <- wrap_plots(list(aic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1) &
      plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
    var_shp <- wrap_plots(var_shp, ncol=3) &
      plot_annotation(title='Pooled Relationship: Variable exposure distribution, by df',
                      theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
    
    modelfitstats <- modelfitstats %>% arrange(AIC)
    num <- as.numeric(gsub("[^0-9.]", "", best_fit))
    num[is.na(num)] <- 1
    
    return(list(modelfitstats=modelfitstats, best_fit=num, aic_plot=aic_plot, var_shp=var_shp))
  }

  res_argvar_tmp <- fit_argvar('argvar_tmp', 'meta_tmp', ref='50%', col='#D74B3C', replace_exact_matches('mx2t', nm))
  # res_argvar_tp  <- fit_argvar('argvar_tp',  'meta_tp', ref='0%', col='#1D98EB', replace_exact_matches('tp',   nm))
  
  return(list(tmp=res_argvar_tmp))
}

# Pooled tests
pooled_argvar_res <- test_argvar(df, outcome, pooled=TRUE)
sens_res$pooled$argvar$res <- pooled_argvar_res
pooled_argvar_aic <- wrap_plots(map(pooled_argvar_res, 'aic_plot'), nrow=1) &
    plot_annotation(title='Pooled AIC: Variable shape',
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
pooled_argvar_dis <- wrap_plots(map(pooled_argvar_res, 'var_shp'), ncol=1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  plot_annotation(title='Pooled AIC: Variable distribution',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

# Site specific tests
city_argvar_res <- list()
for (l in locations){
  sub <- df %>% filter(location==l)
  city_argvar_res[[l]] <- test_argvar(df=sub, outcome, pooled=FALSE)
  print(l)
}
sens_res$city$argvar$res <- city_argvar_res
city_argvar_aic_tmp <- wrap_plots(
  imap(city_argvar_res, \(city, locations) city$tmp$aic_plot + ggtitle(locations)),
  nrow = 2
) &
  plot_annotation(title = paste0('City-specific AIC: ', replace_exact_matches('mx2t', nm), ' shape'),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold')))

# Get city specific distributions
city_argvar_dis <- list()
for (l in locations){
  out <- city_argvar_res[[l]]
  city_argvar_dis[[l]] <- wrap_plots(map(out, 'var_shp'), ncol=1) &
    plot_annotation(title=paste0('City-specific shapes, by df: ', l),
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
  
}

# EXPORT
sens_analysis <- list('5A_pooled_argvar_AIC'=list(p=pooled_argvar_aic, w=6, h=4),
                      '5A_pooled_argvar_DIS'=list(p=pooled_argvar_dis, w=8, h=6),
                      '5B_city_argvar_AIC_mx2t'=list(p=city_argvar_aic_tmp, w=12, h=8))

# AIC
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}
# City distributions only
for (n in names(city_argvar_dis)) {
  title <- paste0('5C_city_argvar_DIS_',n)
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', title, '.png'), 
         plot=city_argvar_dis[[n]], w=8, h=6, dpi=300)
  print(n)
}

# # Updates
# updates <- list(
#   'mx2t'=list('pooled'=3, 'emergency'=1, 'Una'=6, 'Kangra'=3, 'Nahan'=3, 'Hamirpur'=3, 'Mandi'=5, 'Chamba'=4, 'Dharamshala'=5, 'Solan'=3, 'Kullu'=3, 'Reckong Peo'=3),
#   'tp'  =list('pooled'=1, 'emergency'=1, 'Una'=1, 'Kangra'=3, 'Nahan'=1, 'Hamirpur'=1, 'Mandi'=1, 'Chamba'=4, 'Dharamshala'=4, 'Solan'=3, 'Kullu'=1, 'Reckong Peo'=3),
#   'fg10'=list('pooled'=1, 'emergency'=1, 'Una'=1, 'Kangra'=1, 'Nahan'=1, 'Hamirpur'=1, 'Mandi'=1, 'Chamba'=3, 'Dharamshala'=3, 'Solan'=3, 'Kullu'=1, 'Reckong Peo'=1)
# )
# 
# for (v in var){
#   sens_res$pooled$argvar$res[[v]]$best_fit <- updates[[v]]$pooled
#   sens_res$emergency$argvar$res[[v]]$best_fit <- updates[[v]]$emergency
#   for(l in locations){
#     sens_res$city$argvar$res[[v]][[l]]$best_fit <- updates[[v]][[l]]
#   }
# }

# -- Store parameters -- ####

# Get pooled params
final_params <- list()
final_params$pooled$longterm <- sens_res$pooled$longterm$res$best_fit
final_params$pooled$autocorr <- sens_res$pooled$autocorr$res$best_fit
final_params$pooled$arglag_tmp <- sens_res$pooled$arglag$res$tmp$best_fit
final_params$pooled$argvar_tmp <- sens_res$pooled$argvar$res$tmp$best_fit
final_params$pooled$arglag_tp <- sens_res$pooled$arglag$res$tp$best_fit
final_params$pooled$argvar_tp <- sens_res$pooled$argvar$res$tp$best_fit
final_params$pooled$arglag_fg10 <- sens_res$pooled$arglag$res$fg10$best_fit

# Get city params
for (l in locations){
  final_params[[l]]$longterm <- sens_res$city$longterm$res[[l]]$best_fit
  final_params[[l]]$autocorr <- sens_res$city$autocorr$res[[l]]$best_fit
  final_params[[l]]$arglag_tmp <- sens_res$city$arglag$res[[l]]$tmp$best_fit
  final_params[[l]]$argvar_tmp <- sens_res$city$argvar$res[[l]]$tmp$best_fit
  final_params[[l]]$arglag_tp <- sens_res$city$arglag$res[[l]]$tp$best_fit
  final_params[[l]]$argvar_tp <- sens_res$city$argvar$res[[l]]$tp$best_fit
  final_params[[l]]$arglag_fg10 <- sens_res$city$arglag$res[[l]]$fg10$best_fit
  
}

# 6. (Shape) Test shape of response #####

test_shape <- function(df, params, outcome, pooled=TRUE, city='ALL'){
  
  # Get title
  title <- if (pooled) "Pooled: " else paste0(unique(df$location), ": ")
  
  # Get params
  longterm_df <- params$longterm
  autocorr_lag <- params$autocorr
  arglag_tmp <- params$arglag_tmp
  argvar_tmp <- params$argvar_tmp
  arglag_tp <- params$arglag_tp
  arglag_fg10 <- params$arglag_fg10
  
  # Run models
  models <- c('quasi', 'binom')
  test <- list()
  for(d in models){
    output <- run_model(df, locations, outcome, fam_type=d, pooled=pooled,
                            longterm_df=longterm_df, autocorr_lag=autocorr_lag, 
                            arglag_tmp=arglag_tmp, arglag_tp=arglag_tp, arglag_fg10=arglag_fg10,
                            argvar_tmp=argvar_tmp)
    test[[d]]$aic <- data.frame(model=d, AIC=output$meta_comp$sum_aic)
    test[[d]]$dis <- get_results(output, city=city, lags=NULL, outcome=outcome) #$p_er_tmp$p_blup

    print(d)
  }
  
  # Checks
  modelfitstats <- rbind(test$quasi$aic, test$binom$aic)
  best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model
  
  # Plot dist
  qp_plot_titled <- test$quasi$dis + plot_annotation(title=paste0(title, 'Quasipoisson'), theme=theme(plot.title=element_text(hjust=0.5))) + theme(legend.position='none')
  nb_plot_titled <- test$binom$dis + plot_annotation(title=paste0(title, 'Negative binomial'), theme=theme(plot.title=element_text(hjust=0.5))) + theme(legend.position='none')
  dis_plot <- wrap_elements(full=qp_plot_titled) + wrap_elements(full=nb_plot_titled) 
  
  return(list(best_fit=best_fit, dis_plot=dis_plot))
}

# Get pooled results
shape_res <- list()
shape_res$pooled <- test_shape(df=df, params=final_params$pooled, outcome=outcome, pooled=TRUE, city='ALL')
shape_res$pooled$dis_plot
final_params$pooled$shape <- shape_res$pooled$best_fit

# # Get city-specific results
# for (l in locations){
#   sub <- df %>% filter(location==l)
#   shape_res[[l]] <- test_shape(df=sub, params=final_params[[l]], outcome=outcome, pooled=FALSE, city=l)
#   final_params[[l]]$shape <-  shape_res[[l]]$best_fit
#   print(l)
# }

# EXPORT
for (n in names(shape_res)) {
  if(n=='pooled'){
    title <- paste0('6A_pooled_shape')
    } else{
    title <- paste0('6C_city_shape_', n)
    }
  ggsave(filename=paste0(path_report, 'II. dlnm/sensitivity/SM_', title, '.png'),
         plot=shape_res[[n]]$dis_plot, width=8,height=4, dpi=300)
  print(n)
}

# EXPORT PARAMS ####

save(sens_res, file=paste0(path_project, 'Results/ind_sens_params.RData'))
save(final_params, file=paste0(path_project, 'Results/ind_final_params.RData'))

# End timer: 35 minutes 10 seconds (19/06/2026) ####
end <- Sys.time()
dt <- round(as.numeric(end - start, units='secs'))
hours   <- dt %/% 3600
minutes <- (dt %% 3600) %/% 60
seconds <- dt %% 60
sprintf('%d hours %d minutes %.d seconds', hours, minutes, seconds)




#########
