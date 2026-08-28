### ------------------------------------------------------------------------ ### 
###  DNLM – Pooled and subgroup models with meta predictors (HISP)
### ------------------------------------------------------------------------ ### 

# rm(list=ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

####  DEPENDENCIES  ####
source('2_DLNM/2_Utilities.R')

library(dlnm); library(splines); library(mvmeta); library(MASS); library(TTR)
library(ggplot2); library(dplyr); library(purrr); library(patchwork)
library(openxlsx); library(colorspace); library(grid); library(reshape2)
library(cowplot); library(lmtest); library(flextable); library(officer)


# Function to run model with static formula
run_model <- function(data, locations, outcome, fam_type='quasi', pooled=TRUE,
                      longterm_df=10, autocorr_lag=0, lag_window=c(0,10),
                      arglag_tmp=4, arglag_tp=4, arglag_tp_roll=4,
                      argvar_tmp=1, argvar_tp=1, argvar_tp_roll=1){
  
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
  
  # TMP
  ranges_tmp <- t(sapply(datalist,function(x) range(x$mx2t,na.rm=T)))
  cen_tmp <- sapply(datalist,function(x) {
    quantile(x$mx2t,c(0:100)/100,na.rm=T)
  })
  arglag_tmp <- list(fun='ns', df=arglag_tmp)
  argvar_tmp <- if (argvar_tmp == 1) list(fun="lin") else list(fun="ns", df=argvar_tmp)
  df_tmp <- prod(arglag_tmp$df, argvar_tmp$df)
  
  # TP - DAILY
  ranges_tp <- t(sapply(datalist,function(x) range(x$tp,na.rm=T)))
  cen_tp <- sapply(datalist,function(x) {
    quantile(x$tp,c(0:100)/100,na.rm=T)
  })
  arglag_tp <- list(fun='ns', df=arglag_tp)
  argvar_tp <- if (argvar_tp == 1) list(fun="lin") else list(fun="ns", df=argvar_tp)
  df_tp <- prod(arglag_tp$df, argvar_tp$df)
  
  # TP - WEEKLY
  ranges_tp_roll <- t(sapply(datalist,function(x) range(x$tp_roll,na.rm=T)))
  cen_tp_roll <- sapply(datalist,function(x) {
    quantile(x$tp_roll,c(0:100)/100,na.rm=T)
  })
  arglag_tp_roll <- list(fun='ns', df=arglag_tp_roll)
  argvar_tp_roll <- if (argvar_tp_roll == 1) list(fun="lin") else list(fun="ns", df=argvar_tp_roll)
  df_tp_roll <- prod(arglag_tp_roll$df, argvar_tp_roll$df)
  
  
  # 2. Stage 1: Fit the DLNM ----
  
  # BUILT OBJECTS WHERE outputs WILL BE STORED:
  
  Slist_comp <- vector('list',m)
  Slist_model <- vector('list', m)
  names(Slist_comp) <- locs_used
  names(Slist_model) <- locs_used
  
  ymat_tmp <- matrix(NA,m,df_tmp,dimnames=list(locs_used,paste('spl',seq(df_tmp),sep='')))
  Slist_tmp <- vector('list',m)
  names(Slist_tmp) <- locs_used
  
  ymat_tp <- matrix(NA,m,df_tp,dimnames=list(locs_used,paste('spl',seq(df_tp),sep='')))
  Slist_tp <- vector('list',m)
  names(Slist_tp) <- locs_used
  
  ymat_tp_roll <- matrix(NA,m,df_tp_roll,dimnames=list(locs_used,paste('spl',seq(df_tp_roll),sep='')))
  Slist_tp_roll <- vector('list',m)
  names(Slist_tp_roll) <- locs_used
  
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
      
      # CREATE THE SPLINES - CUMULATIVE
      cb.tmp <- crossbasis(sub$mx2t, lag=lag_window, argvar=argvar_tmp, arglag=arglag_tmp, group=sub$year)
      cb.tpD <- crossbasis(sub$tp, lag=c(0,1), argvar=argvar_tp, arglag=arglag_tp, group=sub$year)
      cb.tpW <- crossbasis(sub$tp_roll, lag=c(2,10), argvar=argvar_tp_roll, arglag=arglag_tp_roll, group=sub$year)
      controls <- '+ weekday + holiday + rh + fg10 + pm25'
      
      # Set formula
      form <- as.formula(
        paste0(outcome,
               " ~ cb.tmp + cb.tpD + cb.tpW + ns(time, df=", longterm_df, "*", years, ')',
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
      
      # TMP
      predtmp <- crosspred(cb.tmp, model, cen=cen_tmp['50%',i])
      ymat_tmp[i,] <- predtmp$coef
      Slist_tmp[[i]] <- predtmp$vcov
      
      # TP - DAILY
      predtp <- crosspred(cb.tpD, model, cen=cen_tp['0%',i])
      ymat_tp[i,] <- predtp$coef
      Slist_tp[[i]] <- predtp$vcov
      
      # TP - WEEKLY
      predtp_roll <- crosspred(cb.tpW, model, cen=cen_tp_roll['0%',i])
      ymat_tp_roll[i,] <- predtp_roll $coef
      Slist_tp_roll[[i]] <- predtp_roll$vcov
      
      # EXTRACT COMPARATIVE METRICS
      Slist_comp[[i]] <- comp_metrics(model)
      Slist_model[[i]] <- model
      
      Slist_comp[[i]]$esti_tmp <- unname(predtmp$coefficients[1])
      Slist_comp[[i]]$esti_tp <- unname(predtp$coefficients[1])
      Slist_comp[[i]]$esti_tp_roll <- unname(predtp_roll$coefficients[1])
      
    })
  
  
  # 3. Stage 2: Multivariate Meta-analysis ----
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['avg_tmp']] <- mean(sapply(Slist_comp, function(x) x$esti_tmp))
  meta_comp[['avg_tp']] <- mean(sapply(Slist_comp, function(x) x$esti_tp))
  meta_comp[['avg_tp_roll']] <- mean(sapply(Slist_comp, function(x) x$esti_tp_roll))
  meta_comp[['sum_aic']] <- sum(sapply(Slist_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$acf)))
  
  # Meta coefs
  meta_tmp <- get_meta_analysis(pooled=pooled, ymat=ymat_tmp,  Slist=Slist_tmp, arg_var=argvar_tmp, arg_lag=arglag_tmp, ranges=ranges_tmp, cen=cen_tmp, lag_window=lag_window)
  meta_tp <- get_meta_analysis(pooled=pooled, ymat=ymat_tp,  Slist=Slist_tp, arg_var=argvar_tp, arg_lag=arglag_tp, ranges=ranges_tp, cen=cen_tp, lag_window=c(0,3))
  meta_tp_roll <- get_meta_analysis(pooled=pooled, ymat=ymat_tp_roll,  Slist=Slist_tp_roll, arg_var=argvar_tp_roll, arg_lag=arglag_tp_roll, ranges=ranges_tp_roll, cen=cen_tp_roll, lag_window=c(4,14))
  
  # RESET WARNING
  options(warn=0)
  
  return(list(models=Slist_model, meta_comp=meta_comp, meta_tmp=meta_tmp, meta_tp=meta_tp, meta_tp_roll=meta_tp_roll))
  
}

# Function to get usable tables and figures
get_results <- function(output, city='ALL', lags=NULL, probs='50%', outcome){
  
  # Function to extract model comparison criteria
  ftab <- function(m, mref=NULL) {
    
    # RETURN NOTHING IF m IS NULL
    if (is.null(m)) return(rep(NA_real_, 12))
    
    # HETEROGENEITY AND IC STATS
    q <- qtest(m)
    het <- c(q$Q[1],q$df[1],q$pvalue[1],(q$Q[1]-q$df[1])/q$Q[1]*100,AIC(m),BIC(m))
    # LR TEST (ONLY FOR META-REGRESSION)
    if(!is.null(mref)) {
      lrstat <- -2*(logLik(mref)-logLik(m))[1]
      df <- attr(logLik(m),'df')-attr(logLik(mref),'df')
      pvalue <- 1-pchisq(lrstat,df)
      lr <- c(lrstat,df,pvalue)
    }
    # WALD TEST (ONLY FOR META-REGRESSION)
    if(!is.null(mref)) {
      coef <- coef(m)[-grep('Int',names(coef(m)))]
      vcov <- vcov(m)[-grep('Int',names(coef(m))),-grep('Int',names(coef(m)))]
      waldstat <- coef%*%solve(vcov)%*%coef
      df <- length(coef)
      pvalue <- 1-pchisq(waldstat,df)
      wald <- c(waldstat,df,pvalue)
    }
    # outputS
    if(!is.null(mref)) {
      return(c(het,lr,wald))
    } else return(c(het,rep(NA,6)))
  }
  
  meta_tmp <- output$meta_tmp
  meta_tp <- output$meta_tp
  meta_tp_roll <- output$meta_tp_roll
  
  # Model selection criteria
  tab <- matrix(NA,3,12)
  colnames(tab) <- c('Q','df','p','I-square','AIC','BIC','stat','df','p','stat','df','p')
  rownames(tab) <- c('tmp','tp', 'tp_roll') #, 'fg10', 'rh', 'fg10', 'pm25')
  tab[1,] <- ftab(m=meta_tmp$mv)
  tab[2,] <- ftab(m=meta_tp$mv)
  tab[3,] <- ftab(m=meta_tp_roll$mv)
  
  # Exposure response (with city specs) by lag
  p_er_tmp <- city_exp_res(city=city, lag=lags, meta=meta_tmp, col='#D74B3C', xlab='Max. Temperature (%)', title=NULL, reference='50%')
  p_er_tp <- city_exp_res(city=city, lag=lags, meta=meta_tp, col='#1D98EB', xlab='Precipitation - Daily (%)', title=NULL, reference='60%')
  p_er_tp_roll <- city_exp_res(city=city, lag=lags, meta=meta_tp_roll , col='#1D98EB', xlab='Precipitation - Weekly (%)', title=NULL, reference='20%')
  
  # Lag response (with city specs) by percentile
  q_tmp <- rowMeans(meta_tmp$cen)[probs]
  q_tp <- rowMeans(meta_tp$cen)[probs]
  q_tp_roll <- rowMeans(meta_tp_roll$cen)[probs]
  p_lr_tmp <- city_lag_res(city=city, quant_v=q_tmp, meta=meta_tmp, col='#D74B3C', xlab='Max. Temperature (%)', title=NULL)
  p_lr_tp <- city_lag_res(city=city, quant_v=q_tp, meta=meta_tp, col='#1D98EB', xlab='Precipitation - Daily (%)', title=NULL, reference='0%')
  p_lr_tp_roll <- city_lag_res(city=city, quant_v=q_tp_roll, meta=meta_tp_roll, col='#1D98EB', xlab='Precipitation - Weekly (%)', title=NULL, reference='0%')
  
  return(list(tab=tab,
              p_er_tmp=p_er_tmp, p_er_tp=p_er_tp, p_er_tp_roll=p_er_tp_roll, #p_er_rh=p_er_rh, p_er_pm25=p_er_pm25,
              p_lr_tmp=p_lr_tmp, p_lr_tp=p_lr_tp, p_lr_tp_roll=p_lr_tp_roll)) #p_lr_rh=p_lr_rh, p_lr_pm25=p_lr_pm25))
  
}

format_results <- function(output, outcome, city=NULL){
  
  # OVERALL
  overall <- get_results(output, city='ALL', lags=NULL, outcome=outcome)
  # LAG
  quant_test <- c('10%', '25%', '75%', '90%')
  quant <- list()
  for (q in quant_test){
    quant[[q]] <- get_results(output, city='ALL', lags=NULL, probs=q, outcome=outcome)
  }
  # EXPOSURE
  lag_test <- c(1, 2)
  lags <- list()
  for (l in lag_test){
    lags[[as.character(l)]] <- get_results(output, city='ALL', lags=l, outcome=outcome)
  }
  
  if(!is.null(output$meta_tmp$mv)){
    
    # Format overall plots
    all <- overall[!names(overall) %in% 'tab']
    p_blup_all <- lapply(all, function(x) x[['p_blup']])
    p_blup_er <- p_blup_all[grepl('er', names(p_blup_all))] #|> head(-2)
    p_blup_lr <- p_blup_all[grepl('lr', names(p_blup_all))] #|> head(-2)
    p_er_overall <- wrap_plots(p_blup_er, nrow=1)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    # Format lag plots
    p_lr_tmp <- wrap_plots(lapply(quant, function(x) x$p_lr_tmp$p_blup), ncol=4) & theme(legend.position='none')
    p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    p_lr_tp_roll <- wrap_plots(lapply(quant, function(x) x$p_lr_tp_roll$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    p_lr <- wrap_plots(p_lr_tmp, p_lr_tp, p_lr_tp_roll, ncol=1)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    # Format exposure plots
    p_er_tmp <- wrap_plots(lapply(lags, function(x) x$p_er_tmp$p_blup), ncol=4) & theme(legend.position='none')
    p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    p_er_tp_roll <- wrap_plots(lapply(lags, function(x) x$p_er_tp_roll$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    p_er <- wrap_plots(p_er_tmp, p_er_tp, p_er_tp_roll, ncol=1)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    
    return(list(overall=overall, p_er_overall=p_er_overall,
                quant=quant, p_lr=p_lr, lags=lags, p_er=p_er,
                city=city, city_er=city_er, city_lr=city_lr))
    
  } else {
    
    # Format overall plots
    all <- overall[!names(overall) %in% 'tab']
    p_raw_all <- lapply(all, function(x) x[['p_raw']])
    p_raw_er <- p_raw_all[grepl('er', names(p_raw_all))] #|> head(-2)
    p_raw_lr <- p_raw_all[grepl('lr', names(p_raw_all))] #|> head(-2)
    p_er_overall <- wrap_plots(p_raw_er, nrow=1)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    # Format lag plots
    p_lr_tmp <- wrap_plots(lapply(quant, function(x) x$p_lr_tmp$p_raw), nrow=1) & theme(legend.position='none')
    p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_raw), nrow=1) & theme(legend.position='none', plot.title=element_blank())
    # p_lr_tp_roll <- wrap_plots(lapply(quant, function(x) x$p_lr_tp_roll$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    # p_lr_rh <- wrap_plots(lapply(quant, function(x) x$p_lr_rh$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    # p_lr_pm25 <- wrap_plots(lapply(quant, function(x) x$p_lr_pm25$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    p_lr <- wrap_plots(p_lr_tmp, p_lr_tp, nrow=2)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    # Format exposure plots
    p_er_tmp <- wrap_plots(lapply(lags, function(x) x$p_er_tmp$p_raw), nrow=1) & theme(legend.position='none')
    p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_raw), nrow=1) & theme(legend.position='none', plot.title=element_blank())
    # p_er_tp_roll <- wrap_plots(lapply(lags, function(x) x$p_er_tp_roll$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    # p_er_rh <- wrap_plots(lapply(lags, function(x) x$p_er_rh$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    # p_er_pm25 <- wrap_plots(lapply(lags, function(x) x$p_er_pm25$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
    p_er <- wrap_plots(p_er_tmp, p_er_tp, nrow=2)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
    
    return(list(overall=overall, quant=quant, lags=lags,
                p_er_overall=p_er_overall, p_lr=p_lr, p_er=p_er))
  }
}

run_bestfit <- function(df, locations, outcome, params, pooled, city) {
  
  output  <- run_model(
    df, locations, outcome,
    fam_type    ='binom',
    pooled      =pooled,
    longterm_df =params$longterm,
    autocorr_lag=params$autocorr,
    arglag_tmp  =params$arglag_tmp,
    argvar_tmp  =params$argvar_tmp,
    arglag_tp   =params$arglag_tp,
    arglag_tp_roll =params$arglag_tp_roll
  )
  results <- format_results(output, outcome, city=NULL)
  message('Model AIC: ', output$meta_comp$sum_aic)
  list(output=output, results=results)
}


#### MAIN CODE ####

# Load data
load(paste0(path_project, 'Results/ind_final_params.RData'))
load('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData')
df        <- df_dlnm; rm(df_dlnm)

locations <- df |> distinct(location, elev) |> arrange(elev) |> pull(location) |> as.character()

m <- length(locations)
outcome   <- 'total_patients'

# Lag windows used for exposure plots
tmp_lags <- list(
  ALL         =c(0, 3, 4),    Una        =c(0, 3),
  Kangra      = c(0, 3),      Nahan      =c(1, 8),
  Hamirpur    =c(1, 7),       Mandi= c(1, 4),
  Chamba      =c(1, 3),       Kullu      =c(1, 3),
  Dharamshala =c(0, 3),       Solan      =c(0, 3),
  'Reckong Peo'= c(0, 3)
)
tp_lags <- list(
  ALL         =c(0, 3, 9),    Una        =c(0, 7),
  Kangra      =c(0, 7),       Nahan      =c(0, 7),
  Hamirpur    =c(0, 7),       Mandi      =c(0, 7),
  Chamba      =c(0, 7),       Kullu      =c(0, 7),
  Dharamshala =c(0, 7),       Solan      =c(0, 7),
  'Reckong Peo'= c(0, 7)
)

# Part 1 – Pooled model ----
final_params <- final_params$pooled; gc()
final_params$shape <- 'quasi'
final_params$arglag_fg10 <- NULL
final_params$arglag_tp_roll <- 4

df$tp_roll <- df$tp_roll_log
df$tp <- df$tp_log

res      <- run_bestfit(df, locations, outcome, params=final_params, pooled=TRUE, city=NULL)
output   <- res$output
main_res <- res$results


## Shared helper: build contour plots for all cities + pooled
make_contour_plots <- function(meta, xlab, reference, max_dev=0.20, sig=FALSE) {
  shared_lims <- c(0.999 - max_dev, 1 + max_dev) * 1.001
  shared_brks <- pretty(shared_lims, n=8)
  
  args <- list(meta=meta, xlab=xlab, reference=reference, sig=sig,
               shared_lims=shared_lims, shared_brks=shared_brks)
  
  setNames(
    lapply(c("ALL", locations), function(l)
      do.call(city_con_res, c(args, list(city=l, title=if (l == "ALL") "Pooled" else l)))$p),
    c("ALL", locations)
  )
}

cont_sig <- list(
  tmp=make_contour_plots(output$meta_tmp, "Temperature (%)",   "50%",  sig=TRUE),
  tp =make_contour_plots(output$meta_tp,  "Precipitation - Daily (%)", "60%",  sig=TRUE),
  tp_roll =make_contour_plots(output$meta_tp_roll,  "Precipitation - Weekly (%)", "20%",  sig=TRUE)
)

# Inspect significant lags
wrap_plots(cont_sig$tmp, nrow=3) +
  plot_layout(guides='collect') +
  plot_annotation(
    title='Significant exposure-lag-response',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
  ) & theme(legend.position='bottom')

wrap_plots(cont_sig$tp, nrow=3) +
  plot_layout(guides='collect') +
  plot_annotation(
    title='Significant exposure-lag-response',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
  ) & theme(legend.position='bottom')

wrap_plots(cont_sig$tp_roll, nrow=3) +
  plot_layout(guides='collect') +
  plot_annotation(
    title='Significant exposure-lag-response',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
  ) & theme(legend.position='bottom')
