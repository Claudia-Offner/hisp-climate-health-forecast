### ------------------------------------------------------------------------ ### 
### DLNM - Sensitivity Analysis
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES  ####

# Function for univariate parameter testing
run_param <- function(data, locations, outcome, fam_type='binom', pooled=TRUE, 
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
    if(var_test == 'tp'){
      ref <- get_tp_ref(data, locations)
    } else{
      ref <- '50%'
    }
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
      
      # EXTRACT COMPARATIVE METRICS
      Slist_comp[[i]] <- comp_metrics(model)
      Slist_model[[i]] <- model
      if(!is.null(var_test) && control_test == 'Distributed Lag'){
        pred<- crosspred(cb, model, cen=cen_var[ref,i])
        Slist_comp[[i]]$esti_var <- unname(pred$coefficients[1])
      }
      
    })
  
  # 3. Stage 2: Multivariate Meta-analysis ----

  # Meta metrics
  meta_comp <- list()
  meta_comp[['sum_qaic']] <- sum(sapply(Slist_comp, function(x) x$qaic))
  meta_comp[['sum_acf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$acf)))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_esti']] <- mean(sapply(Slist_comp, function(x) as.numeric(x$esti_var)))
  
  # RESET WARNING
  options(warn=0)
  
  return(list(models=Slist_model, meta_comp=meta_comp))
  
}

# Function to plot QAIC results (comparison)
plot_QAIC_res <- function(modelfitstats, title, best_fit, xlab='Degrees of freedom'){
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, 'yes', 'no')
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels=modelfitstats$model)
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x=model, y=QAIC, group=1)) +
    geom_line(color='grey20') +
    geom_point(aes(color=highlight), size=3) +
    scale_color_manual(values=c('no'='black', 'yes'='red')) +
    labs(title=title,
         y='Sum of QAIC', x=xlab) +
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
  
  best_fit <- if (mean(modelfitstats$QAIC) < 0) {
    modelfitstats[which.max(modelfitstats$QAIC), ]$model
  } else {
    modelfitstats[which.min(modelfitstats$QAIC), ]$model
  }
  
  modelfitstats$highlight <- ifelse(modelfitstats$model == best_fit, 'yes', 'no')
  
  # Preserve original order
  modelfitstats$model <- factor(modelfitstats$model, levels=modelfitstats$model)
  
  
  # Plot with highlighted point
  p <- ggplot(modelfitstats, aes(x=model, y=esti, group=1)) +
    geom_hline(yintercept=mean(modelfitstats$esti), color='red') +
    geom_line(color='grey20') +
    geom_point(size=3) +
    # geom_point(aes(color=highlight), size=3) +
    # scale_color_manual(values=c('no'='black', 'yes'='red')) +
    labs(title=replace_exact_matches(exp, nm), y='Average of estimates', x='Days') +
    theme_classic() +
    theme(
      axis.line=element_line(color='gray20'),
      axis.text=element_text(color='gray20'),
      # axis.text.x=element_text(angle=45, hjust=1),
      axis.ticks=element_line(color='gray20'),
      plot.title=element_text(face='bold', hjust=0.5, size=12),
      legend.position='none'
    )
  return(p)
  
}

#### MAIN CODE ####

export <- paste0(a, '_', model)
df$groups <- df$location

# Store results
sens_res <- list()

# 1.  (Longterm) Test time df for capturing long-term/seasonal ####

test_longterm <- function(df, outcome, pooled=TRUE){
  
  if(length(unique(df$location)) !=1 ){
    title <-  'Pooled QAIC & Estimates: Long-term seasonality'
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
      QAIC   =test_result$meta_comp$sum_qaic
      # avg_tmp=test_result$meta_comp$avg_tmp,
      # avg_tp =test_result$meta_comp$avg_tp
    ))
  }
  
  # Get QAIC
  best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
  qaic_plot <- plot_QAIC_res(modelfitstats, title, best_fit)
  # # Get avg tmp
  # esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
  # esti_tp_plot <- plot_ESTI_res(modelfitstats, exp='tp')
  # Combine plots
  # all_plot <- wrap_plots(list(qaic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1)  &
  #   plot_annotation(title=title,theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
  
  modelfitstats <- modelfitstats %>% arrange(QAIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, best_fit=num, qaic_plot=qaic_plot))
}

# Pooled tests
sens_res$pooled$longterm$res <- test_longterm(df, outcome, pooled=TRUE)
pool_longterm_qaic <- sens_res$pooled$longterm$res$qaic_plot

# EXPORT
sens_analysis <- list('1A'=list(p=pool_longterm_qaic, w=6, h=4))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export,'/S3_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# 2.  (Autocorr) Test autocorrelation of hospital admissions #####

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
    title <- 'Pooled QAIC'
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
      QAIC    =test_result$meta_comp$sum_qaic
      # avg_tmp=test_result$meta_comp$avg_tmp,
      # avg_tp =test_result$meta_comp$avg_tp
    ))
  }
  
  # Get QAIC
  best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
  qaic_plot <- plot_QAIC_res(modelfitstats, title, best_fit, 'Autocorrelation Lag')
  # Get avg tmp
  # esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
  # esti_tp_plot <- plot_ESTI_res(modelfitstats, exp='tp')
  # # Combine plots
  # all_plot <- wrap_plots(list(qaic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1)  &
  #   plot_annotation(title=title,theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
  
  modelfitstats <- modelfitstats %>% arrange(QAIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  # Visualize final PACF plot
  test_result <- run_param(data=df, locations=locations, outcome=outcome, pooled=pooled, autocorr_lag=num)
  acf_plot_fin <- plot_pooled_acf(test_result, plot='ACF', name=name, num)
  pacf_plot_fin <- plot_pooled_acf(test_result, plot='PACF', name=name, num)
  
  QAIC_ACF <- list(); QAIC_PACF <- list()
  QAIC_ACF$raw <- acf_plot_raw
  QAIC_ACF$fin <- acf_plot_fin
  QAIC_PACF$raw <- pacf_plot_raw
  QAIC_PACF$fin <- pacf_plot_fin
  QAIC_ACF$all <- acf_plot_raw + acf_plot_fin + plot_layout(ncol=2)
  QAIC_PACF$all <- pacf_plot_raw + pacf_plot_fin + plot_layout(ncol=2)
  
  return(list(modelfitstats=modelfitstats,  best_fit=num, qaic_plot=qaic_plot,
              QAIC_ACF=QAIC_ACF, QAIC_PACF=QAIC_PACF))
}

# Pooled tests
sens_res$pooled$autocorr$res <- test_autocorr(df, outcome, pooled=TRUE)

# ALL plots for appendix
all <- list(sens_res$pooled$autocorr$res$QAIC_ACF$raw, sens_res$pooled$autocorr$res$qaic_plot, sens_res$pooled$autocorr$res$QAIC_ACF$fin)
pooled_autocorr <- (wrap_plots(all, nrow=1))  +
  plot_annotation(
    tag_levels = 'A',
    theme=theme(plot.title=element_text(hjust=0.5, face="bold"))
  )

# EXPORT
sens_analysis <- list('2A'=list(p=pooled_autocorr, w=12, h=4))

for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export,'/S3_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# 3.  (Var) Test added value of each weather variable within model ####

test_var <- function (df, outcome, exp, pooled=TRUE, tmp_title='tmp'){

  if(pooled==TRUE){
    title <-  paste0('Pooled QAIC & Estimates: Time-varying confounders (', 
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
    C0=list(vartest='None'),
    C1=list(vartest='Single Lag',                 arg_var=1),
    C2=list(vartest='Moving Average',             arg_var=1),
    C3=list(vartest='Exponential Moving Average', arg_var=1),
    C4=list(vartest='Distributed Lag',            arg_var=1, arg_lag=lag_df),
    C5=list(vartest='Single Lag',                 arg_var=var_df),
    C6=list(vartest='Moving Average',             arg_var=var_df),
    C7=list(vartest='Exponential Moving Average', arg_var=var_df),
    C8=list(vartest='Distributed Lag',            arg_var=var_df, arg_lag=lag_df)
  )
  
  modelfitstats <- data.frame()
  for (cname in names(settings)) {
    
    s <- settings[[cname]]
    cat('\n[', cname, '-', s$vartest, ']\n')
    
    if (cname == 'C0') {

      test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                               pooled=pooled, longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                               control_test=s$vartest, var_test=exp)
      
    } else {
      test_result <- run_param(data=df, locations=locations, outcome=outcome, 
                               pooled=pooled, longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                               control_test=s$vartest, var_test=exp, 
                               argvar_var=s$arg_var, arglag_var=s$arg_lag)
    }
    
    modelfitstats <- rbind(modelfitstats, data.frame(
      model  =cname,
      QAIC    =test_result$meta_comp$sum_qaic
    )) 
    
  }

  # Get QAIC
  best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
  qaic_plot <- plot_QAIC_res(modelfitstats, title, best_fit, xlab='Model')
  modelfitstats <- modelfitstats %>% arrange(QAIC)

  return(list(modelfitstats=modelfitstats, best_fit=best_fit, qaic_plot=qaic_plot))

}

controls <- c('tmp', 'tp')
results <- list()
for (c in controls){
  
  print(c)
  results[[c]] <- test_var(df, outcome, exp=c, pooled=TRUE, tmp_title=a)
  
}

sens_res$pooled$confound$res <- results
plots <- map(results, 'qaic_plot')
controls[1] <- a
titles <- replace_exact_matches(controls, nm)
pool_confound_qaic <- wrap_plots(map2(plots, titles, \(p, t) p + ggtitle(t)), nrow=1) &
  plot_annotation(title='Pooled QAIC: Time-varying exposures',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

# EXPORT
sens_analysis <- list('3B'=list(p=pool_confound_qaic, w=9, h=3))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export,'/S3_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# 4a. (Lag_window) Test feasible lag windows ####

test_window <- function(df, outcome, pooled=TRUE, exp='t2m'){
  
  if(pooled==TRUE){
    if(exp=='tmp'){t <- a} else {t <- exp}
    title <- replace_exact_matches(t, nm)
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
    autocorr_lag <- sens_res$city$longterm$res[[title]]$best_fit
  }
  
  modelfitstats <- data.frame()
  for (i in 3:28){
    test_result <- run_param(data=df, locations=locations, outcome=outcome, pooled=pooled, 
                             longterm_df=longterm_df, autocorr_lag=autocorr_lag,
                             control_test='Distributed Lag', var_test=exp,
                             lag_window=c(0,i))
    modelfitstats <- rbind(modelfitstats, data.frame(
      model  =i,
      QAIC    =test_result$meta_comp$sum_qaic,
      esti   =test_result$meta_comp$sum_esti
    ))
  }
  
  # Get QAIC
  best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
  qaic_plot <- plot_QAIC_res(modelfitstats, title=title, best_fit)
  esti_plot <- plot_ESTI_res(modelfitstats, exp)

  modelfitstats <- modelfitstats %>% arrange(QAIC)
  num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
  num[is.na(num)] <- 1
  
  return(list(modelfitstats=modelfitstats, best_fit=num, qaic_plot=qaic_plot, esti_plot=esti_plot))
}

exp <- c('tmp', 'tp')
pooled_window_lag <- list()
for (i in exp){
  pooled_window_lag[[i]] <- test_window(df, outcome, pooled=TRUE, exp=i)
}
pooled_window_esti <- wrap_plots(imap(pooled_window_lag, ~ .x$esti), nrow=1) &
  plot_annotation(title='Pooled Estimations: Lag window',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
pooled_window_qaic <- wrap_plots(imap(pooled_window_lag, ~ .x$qaic), nrow=1) &
  plot_annotation(title='Pooled QAIC: Lag window',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

sens_analysis <- list('4A'=list(p=pooled_window_esti, w=12, h=4))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export,'/S3_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

sens_res$pooled$lag_window$res$best_fit <- 10


# 4b. (ArgLag) Test single weather variables for lag structures ####

test_arglag <- function(df, outcome, pooled, tmp_title='tmp'){
  
  # Select best params from previous tests
  if(pooled==TRUE){
    title <- 'Pooled QAIC: Variable lag'
    longterm_df <- sens_res$pooled$longterm$res$best_fit
    autocorr_lag <- sens_res$pooled$autocorr$res$best_fit
  } else{
    title <- as.character(unique(df$location))
    longterm_df <- sens_res$city$longterm$res[[title]]$best_fit
    autocorr_lag <- sens_res$city$autocorr$res[[title]]$best_fit
  }
  
  # Get tmp and tp lag params
  fit_arglag <- function(param, ex) {
    modelfitstats <- purrr::map_dfr(3:6, \(i) {
      test_result <- run_model(
        df, groups=locations, outcome = outcome, fam_type = 'binom',
        longterm_df = longterm_df, autocorr_lag = autocorr_lag, lag_window = c(0, 10),
        arglag_tmp = if (param == 'arglag_tmp') i else 3, argvar_tmp = 1,
        arglag_tp  = if (param == 'arglag_tp')  i else 3, argvar_tp  = 1,
        group_meta = NULL, subgroup_var = NULL
      )
      data.frame(model = i, QAIC = test_result$meta_comp$sum_qaic,
                 avg_tmp = test_result$meta_comp$avg_tmp, avg_tp = test_result$meta_comp$avg_tp)
    })
    
    best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
    qaic_plot <- plot_QAIC_res(modelfitstats, ex, best_fit)
    modelfitstats <- modelfitstats %>% arrange(QAIC)
    num <- as.numeric(sub('.*\\((\\d+)\\).*', '\\1', best_fit))
    num[is.na(num)] <- 1
    
    return(list(modelfitstats=modelfitstats, best_fit=num, qaic_plot=qaic_plot))
  }
  
  modelfitstats_tmp <- fit_arglag('arglag_tmp', replace_exact_matches(tmp_title,nm))
  modelfitstats_tp  <- fit_arglag('arglag_tp',  replace_exact_matches('tp',nm))

  return(list(tmp=modelfitstats_tmp, tp=modelfitstats_tp))
}

# Pooled tests
pooled_arglag_res <- test_arglag(df, outcome, pooled=TRUE, tmp_title=a)
sens_res$pooled$arglag$res <- pooled_arglag_res
pooled_arglag_qaic <- wrap_plots(map(pooled_arglag_res, 'qaic_plot'), nrow=1) &
   plot_annotation(title='Pooled QAIC: Variable lag',
                   theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

# EXPORT
sens_analysis <- list('4B'=list(p=pooled_arglag_qaic, w=9, h=3))
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export,'/S3_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# 5.  (ArgVar) Test single weather variables for shape structures ####

test_argvar <- function(df, outcome, pooled, tmp_title='tmp'){
  
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
  
  fit_argvar <- function(param, meta_key, ref, col, title, xlab) {
    
    run_one <- function(i) {
      run_model(
        df, groups=locations, outcome = outcome, fam_type = 'binom',
        longterm_df = longterm_df, autocorr_lag = autocorr_lag, lag_window = c(0, 10),
        arglag_tmp = arglag_tmp_df, arglag_tp = arglag_tp_df,
        argvar_tmp = if (param == 'argvar_tmp') i else 1,
        argvar_tp  = if (param == 'argvar_tp')  i else 1
      )
    }
    
    dfs <- c(Linear = 1, setNames(3:6, paste0('NS', 3:6)))
    results <- purrr::imap(dfs, \(i, name) {
      res <- run_one(i)
      list(
        stats = data.frame(model = name, QAIC = res$meta_comp$sum_qaic,
                           avg_tmp = res$meta_comp$avg_tmp, avg_tp = res$meta_comp$avg_tp),
        shp = city_exp_res(city = city, lag = NULL, meta = res[[meta_key]], col = col,
                           xlab = if (name == 'Linear') xlab else 'Percentiles (%)',
                           title = name, reference = ref)[[p]]
      )
    })
    
    modelfitstats <- purrr::map_dfr(results, 'stats')
    var_shp <- purrr::map(results, 'shp')
    
    best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
    qaic_plot <- plot_QAIC_res(modelfitstats, title, best_fit)
    esti_tmp_plot <- plot_ESTI_res(modelfitstats, exp='tmp')
    esti_tp_plot  <- plot_ESTI_res(modelfitstats, exp='tp')
    all_plot <- wrap_plots(list(qaic_plot+theme(plot.title=element_blank()), esti_tmp_plot, esti_tp_plot), nrow=1) &
      plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
    var_shp <- wrap_plots(var_shp, nrow=1) &
      plot_annotation(title='Pooled Relationship: Variable exposure distribution, by df',
                      theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
    
    modelfitstats <- modelfitstats %>% arrange(QAIC)
    num <- as.numeric(gsub("[^0-9.]", "", best_fit))
    num[is.na(num)] <- 1
    
    return(list(modelfitstats=modelfitstats, best_fit=num, qaic_plot=qaic_plot, var_shp=var_shp))
  }

  res_argvar_tmp <- fit_argvar('argvar_tmp', 'meta_tmp', ref='50%', col='#D74B3C', replace_exact_matches(tmp_title, nm), 'Temperature (%)')
  res_argvar_tp  <- fit_argvar('argvar_tp',  'meta_tp', ref=tp_ref, col='#1D98EB', replace_exact_matches('tp',   nm), 'Percipitation (%)')
  
  return(list(tmp=res_argvar_tmp, tp=res_argvar_tp))
}

# Pooled tests
pooled_argvar_res <- test_argvar(df, outcome, pooled=TRUE, tmp_title=a)
sens_res$pooled$argvar$res <- pooled_argvar_res
pooled_argvar_qaic <- wrap_plots(map(pooled_argvar_res, 'qaic_plot'), nrow=1) &
    plot_annotation(title='Pooled QAIC: Variable shape',
                    theme=theme(plot.title=element_text(hjust=0.5, face='bold')))
pooled_argvar_dis <- wrap_plots(map(pooled_argvar_res, 'var_shp'), ncol=1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none") &
  plot_annotation(title='Pooled QAIC: Variable distribution',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold')))

# EXPORT
sens_analysis <- list('5A'=list(p=pooled_argvar_qaic, w=9, h=3),
                      '5B'=list(p=pooled_argvar_dis, w=15, h=5))
# QAIC
for (n in names(sens_analysis)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export,'/S3_', n, '.png'), 
         plot=sens_analysis[[n]]$p, width=sens_analysis[[n]]$w,height=sens_analysis[[n]]$h, dpi=300)
  print(n)
}

# --- Store parameters --- ####

# Get pooled params
final_params <- list()
final_params$pooled$longterm <- sens_res$pooled$longterm$res$best_fit
final_params$pooled$autocorr <- sens_res$pooled$autocorr$res$best_fit
final_params$pooled$lag_window  <- sens_res$pooled$lag_window$res$best_fit
final_params$pooled$arglag_tmp <- sens_res$pooled$arglag$res$tmp$best_fit
final_params$pooled$argvar_tmp <- 3
final_params$pooled$arglag_tp <- sens_res$pooled$arglag$res$tp$best_fit
final_params$pooled$argvar_tp <- 1

# 7.  (Shape) Test shape of response #####

test_shape <- function(df, params, outcome, pooled=TRUE, city='ALL'){
  
  # Get title
  title <- if (pooled) "Pooled: " else paste0(unique(df$location), ": ")
  
  # Get params
  longterm_df <- params$longterm
  autocorr_lag <- params$autocorr
  arglag_tmp <- params$arglag_tmp
  argvar_tmp <- params$argvar_tmp
  arglag_tp <- params$arglag_tp

    # Run models
  models <- c('quasi', 'binom')
  shapetest <- list()
  for(d in models){
    output <- run_model(df, groups=locations, outcome = outcome, fam_type=d,
                        longterm_df = longterm_df, autocorr_lag = autocorr_lag, lag_window = c(0, 10),
                        arglag_tmp=arglag_tmp, arglag_tp=arglag_tp,
                        argvar_tmp=argvar_tmp)
    shapetest[[d]]$qaic <- data.frame(model=d, QAIC=output$meta_comp$sum_qaic)
    
    # Exposure response (with city specs) by lag
    shapetest[[d]]$dis$p_er_tmp <- city_exp_res(city=city, lag=NULL, meta=output$meta_tmp, col='#D74B3C', xlab='Temperature (%)', title=NULL, reference='50%')$p_blup
    shapetest[[d]]$dis$p_er_tp <- city_exp_res(city=city, lag=NULL, meta=output$meta_tp, col='#1D98EB', xlab='Precipitation (%)', title=NULL, reference=tp_ref)$p_blup

    print(d)
  }
  
  # Checks
  modelfitstats <- rbind(shapetest$quasi$qaic, shapetest$binom$qaic)
  best_fit <- modelfitstats[which.min(modelfitstats$QAIC), ]$model
  
  qp_plot_titled <- wrap_plots(shapetest$quasi$dis, nrow=1) + 
    plot_annotation(title=paste0(title, 'Quasipoisson'), theme=theme(plot.title=element_text(hjust=0.5, face='bold'))) &
    theme(legend.position='none')
  nb_plot_titled <- wrap_plots(shapetest$binom$dis, nrow=1) + 
    plot_annotation(title=paste0(title, 'Negative binomial'), theme=theme(plot.title=element_text(hjust=0.5, face='bold'))) &
    theme(legend.position='none')
  dis_plot <- wrap_elements(full=qp_plot_titled) / wrap_elements(full=nb_plot_titled) 
  
  return(list(best_fit=best_fit, dis_plot=dis_plot))
}

# Get pooled results
shape_res <- list()
shape_res$pooled <- test_shape(df=df, params=final_params$pooled, outcome=outcome, pooled=TRUE, city='ALL')
shape_res$pooled$dis_plot
final_params$pooled$shape <- 'binom'

# EXPORT
for (n in names(shape_res)) {
  ggsave(filename=paste0(path_report, 'sensitivity/', export, '/S3_6A.png'), 
         plot=shape_res[[n]]$dis_plot, width=9,height=6, dpi=300)
  print(n)
}

#### EXPORT PARAMS ####

save(final_params, file=paste0(path_git, 'outputs/', export, '_final_params.RData'))


#### End timer: ~35 minutes (19/06/2026) ####
end <- Sys.time()
dt <- round(as.numeric(end - start, units='secs'))
hours   <- dt %/% 3600
minutes <- (dt %% 3600) %/% 60
seconds <- dt %% 60
sprintf('%d hours %d minutes %.d seconds', hours, minutes, seconds)



