### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city Sensitivity  Analysis                                 ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta)
library(ggplot2); library(dplyr); library(purrr); library(patchwork)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

#### FUNCTIONS ####
library('MuMIn')



# Function to get reference meta-analysis
get_meta_analysis <- function(ymat, Slist, arg, ranges, cen){
  
  # Meta-analysis w/ prediction basis
  mv <- mvmeta(ymat, Slist, method="ml")
  bound <- c(colMeans(ranges))
  mean <- seq(bound[1],bound[2],length=30)
  b_mean <- crossbasis(mean, lag=c(0,14),argvar=list(fun='ns', df=arg$argvar), arglag=list(fun='ns', df=arg$arglag))
  
  return(list(ymat=ymat, Slist=Slist, arg=arg, ranges=ranges, mv=mv, bound=bound, b_mean=b_mean, cen=cen))
  
}

# Function to calculate comparative metrics for meta-DLNM
comp_metrics <- function(model, base=T) {
  
  res <- list()
  
  # Get sum aic
  res[['aic']] <- summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion
  
  # Get sum of PACF
  res[['pacf']] <- pacf(residuals(model, type = "pearson"), plot = FALSE, lag.max = 28)$acf
  
  if(base==F){
    
    # Get average of estimates (lag 4)
    res[['esti_t2m']] <- predt2m$coef[4]  
    res[['esti_tp']]  <- predtp$coef[4]  
    res[['esti_rh']] <- predrh$coef[4]  
    res[['esti_pm25']]  <- predpm25$coef[4] 
    
  }
  
  return(res) 
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

# Function for parameter testing
get_parameter_test <- function(data, locations, test='Longterm', variables, 
                               auto_lag=0, v_dist=NULL, v_lag=14, df_lag=4, time_df=10){
  
  # SET PARAMETERS
  datalist <- lapply(locations, function(region) data[data$location==region,])
  names(datalist) <- locations
  years <- length(unique(data$year))
  
  # EMPTY OBJECTS TO STORE OUTPUTS
  list_comp <- vector("list",m)
  names(list_comp) <- locations
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  system.time(
    for(i in seq(m)) {
      
      # PRINT ITERATION
      cat(i,"")
      
      # LOAD
      data <- datalist[[i]]
      
      # RUN THE MODEL & EXTRACT COMPARATIVE METRICS
      if(test=='Longterm'){
        baseformula <- paste0('total_patients ~  season + weekday + holiday + ns(time, df=', time_df,'*years)')
        model <- glm(formula=baseformula,family=quasipoisson(),data)
        list_comp[[i]] <- comp_metrics(model, base=T)
      } else if (test == 'Autocorrelation') {
        
        if (auto_lag==0) {
          baseformula <- paste0('total_patients ~  season + weekday + holiday + ns(time, df=', time_df, '*years)')
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model, base=T)
        } else{
          # Add lags to df
          data <- data %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, auto_lag))
          # Fit multivariate model
          baseformula <- paste0('total_patients ~  season + weekday + holiday + ns(time, df=', time_df, '*years) + patients_lag')
          model <- glm(formula=baseformula,family=quasipoisson(),data)
          list_comp[[i]] <- comp_metrics(model, base=T)
        }
      } 
    })
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['sum_aic']] <- sum(sapply(list_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowSums(sapply(list_comp, function(x) as.numeric(x$pacf)))
  
  # Format AIC
  if(test=='Longterm'){
    aic_result <- data.frame(paste0(test, ' (',time_df,')'), meta_comp[['sum_aic']])
  } else if (test == 'Autocorrelation') {
    aic_result <- data.frame(paste0(test, ' (',auto_lag,')'), meta_comp[['sum_aic']])
  }
  names(aic_result) <- c('model','AIC' )
  
  
  # # LONG-TERM TESTS
  # if(baseline==TRUE){
  #   # Only fit the baseline model
  #   model <- glm(formula=as.formula(baseformula), data=df,  family = quasipoisson(link="log"), na.action='na.exclude') 
  #   test_result <- get_model_fit(model,paste0('baseline (',time_df,')'))
  # } else{
  #   
  #   # LAG/SHAPE TESTS
  #   if(length(variables)==1){
  #     # Fit univariate model
  #     cb <- crossbasis(df[[variables]], lag=c(0,v_lag),argvar=v_dist, arglag=list(fun='ns', df=df_lag))
  #     model <- glm(formula=as.formula(paste(baseformula,'+ cb')),data=df, family = quasipoisson(link="log"), na.action='na.exclude') 
  #     test_result <- get_model_fit(model,paste0(variables,'(',v_dist$fun,'[',v_dist$df , '], ', df_lag,')'))
  #     
  #     # AUTOCORR TESTS
  #   } else{
  #     # Fit multivariate model
  #     cb.t2m <- crossbasis(df$t2m, lag=c(0,v_lag),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=df_lag))
  #     cb.tp <- crossbasis(df$tp, lag=c(0,v_lag),argvar=list(fun='lin'), arglag=list(fun='ns', df=df_lag))
  #     cb.fg10 <- crossbasis(df$fg10, lag=c(0,v_lag),argvar=list(fun='lin'), arglag=list(fun='ns', df=df_lag))
  #     model <- glm(formula=as.formula(paste(baseformula,'+ cb.t2m + cb.tp + cb.fg10')), data=df, family = quasipoisson(link="log"), na.action='na.exclude') 
  #     test_result <- get_model_fit(model,paste0('autocorr(', 1 ,')'))
  #   }
  
  # }
  
  return(list(aic_result=aic_result, pacf_result=meta_comp[['sum_pacf']], 
              model=model))
  
}



#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)


# 1. (Baseline) Test time df for capturing long-term/seasonal (time_df=9) ####

setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 1:15){
  test_result <- get_parameter_test(data=df, locations=locations, test='Longterm', time_df=i)
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_season <- plot_AIC_res(modelfitstats, best_fit, "Model AIC (seasonality)"))


# 2. (Autocorr) Test autocorrelation of hospital admissions (AR=1) #####

setNames(modelfitstats <- data.frame(matrix(NA,0,2)),c('model', 'AIC' ))
for(i in 0:14){
  test_result <- get_parameter_test(data=df, locations=locations, test='Autocorrelation', auto_lag=i)
  modelfitstats <- rbind(modelfitstats, test_result$aic_result)
}
modelfitstats
(best_fit <- modelfitstats[which.min(modelfitstats$AIC), ]$model)
(AIC_season <- plot_AIC_res(modelfitstats, best_fit, "Model AIC (seasonality)"))

# Visualize PACF plot
test_result <- get_parameter_test(data=df, locations=locations, test='Autocorrelation', auto_lag=0)
pacf_plot <- data.frame(
  lag = 1:length(test_result$pacf_result),
  pacf = test_result$pacf_result
)
ggplot(pacf_plot, aes(x = lag, y = pacf)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
  scale_y_continuous(limits = c(min(test_result$pacf_result, -0.2), max(test_result$pacf_result, 1))) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "gray20"),
    axis.text = element_text(color = "gray20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "gray20"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  ggtitle("Pooled PACF Across Cities")




# 3. (Lag) Test single weather variables for lag structures ####

# 4. (Shape) Test single weather variables for shape structures ####

# 5. (Shape) Test shape of response (negative binom) #####

