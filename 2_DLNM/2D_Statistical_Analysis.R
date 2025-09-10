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
# install.packages('dlnm')
library(dlnm); library(splines); library(dplyr); library(lubridate); 
library(MASS); library(patchwork); library(ggplot2); library(ggfortify);
library(gridExtra)

# Variables
nm <- c(total_patients='Total patients', M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        general='General medicine', emergency='Emergency', ophthalmology='Ophthalmology', 
        pediatrics='Paediatric', orthopedics='Orthopedics', ent='ENT', 
        gynecology='Gynaecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        t2m='T2m (\u00B0C)', mx2t='Mx2t (\u00B0C)', mn2t='Mn2t (\u00B0C)', tp='Tp (mm)', fg10='Fg10 (km/h)', 
        closed='Closed', holiday='Holiday')

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

# Function to retrieve results from each predictor, with visuals and tables
get_results <- function(df, outcome, quant_v=0.95, lag_v=c(NULL, NULL, NULL)){
  
  # Fit cross-basis
  cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df$year)
  # cb.t2m <- crossbasis(df$mx2t, lag=c(0,14), argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=8))
  # cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=4))
  # cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=5))
  
  # Fit model
  model <- glm(formula=as.formula(paste0(outcome, ' ~  season + weekday + holiday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag')), data = df, family = quasipoisson(link="log"))
  
  cen_tem <- 12.5
  cen_tp <- 30
  cen_fg10 <- 40
  
  if(is.null(quant_v)){
    lr_t2m <- crosspred(cb.t2m, model, cen=cen_tem)
    lr_tp <- crosspred(cb.tp, model, cen=cen_tp)
    lr_fg10 <- crosspred(cb.fg10, model, cen=cen_fg10)
    
  } else{
    # Select quantile
    q_t2m <- quantile(df$t2m, probs = quant_v, na.rm = TRUE)
    q_tp <- quantile(df$tp, probs = quant_v, na.rm = TRUE)
    q_fg10 <- quantile(df$fg10, probs = quant_v, na.rm = TRUE)
    
    # DISPLAY OVERALL ESTIMATE OF RR 
    lr_t2m <- crosspred(cb.t2m, model, at=q_t2m, cen=cen_tem,)
    lr_tp <- crosspred(cb.tp, model, at=q_tp, cen=cen_tp)
    lr_fg10 <- crosspred(cb.fg10, model, at=q_fg10, cen=cen_fg10)
  }
  
  # Get overall values
  cp <- data.frame(rbind(lr_t2m[c('predvar', 'allRRfit','allRRlow','allRRhigh')], 
                         lr_tp[c('predvar', 'allRRfit','allRRlow','allRRhigh')], 
                         lr_fg10[c('predvar', 'allRRfit','allRRlow','allRRhigh')]))
  df_cp <- data.frame(variables = rownames(cp), cp, row.names = NULL)
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK T2M value
  p_lr_t2m <- get_visuals(lr_t2m, lag_rep=TRUE, col="#D74B3C", title='t2m')
  lr_t2m <- data.frame(Variable= 't2m', Lag  = c(colnames(lr_t2m$matRRfit)), RR  = as.numeric(t(lr_t2m$matRRfit)),
                       Lower  = as.numeric(t(lr_t2m$matRRlow)), Upper = as.numeric(t(lr_t2m$matRRhigh)))
  # GRAPH CUMULATIVE T2M-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & FG10) 
  er_t2m <- crosspred(cb.t2m, model, at=seq(-10, 35, by=1), cen=cen_tem)
  p_er_t2m <- get_visuals(er_t2m, lag_rep=FALSE, lag=lag_v[1], col="#D74B3C", xlab="Temperature (C)", title='t2m')
  l_er_t2m <- data.frame(Variable= 't2m', Exposure = er_t2m$predvar, RR = er_t2m$matRRfit[, lag_v[1] + 1], Lower = er_t2m$matRRlow[, lag_v[1] + 1], Upper = er_t2m$matRRhigh[, lag_v[1] + 1], row.names = NULL)
  o_er_t2m <- data.frame(Variable= 't2m', Lag  = names(er_t2m$allRRfit), RR  = as.numeric(t(er_t2m$allRRfit)),
                         Lower  = as.numeric(t(er_t2m$allRRlow)), Upper = as.numeric(t(er_t2m$allRRhigh)))
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK TP value
  p_lr_tp <- get_visuals(lr_tp, lag_rep=TRUE, col="#1D98EB", title='tp')
  lr_tp <- data.frame(Variable= 'tp', Lag  = c(colnames(lr_tp$matRRfit)), RR  = as.numeric(t(lr_tp$matRRfit)),
                      Lower  = as.numeric(t(lr_tp$matRRlow)), Upper = as.numeric(t(lr_tp$matRRhigh)))
  # GRAPH CUMULATIVE TP-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR T2M & FG10) 
  er_tp <- crosspred(cb.tp, model, at=seq(0, 50, by=1), cen=cen_tp)
  p_er_tp <- get_visuals(er_tp, lag_rep=FALSE, lag=lag_v[2], col="#1D98EB", xlab="Rainfall (mm)", title='tp')
  l_er_tp <- data.frame(Variable= 'tp', Exposure = er_tp$predvar, RR = er_tp$matRRfit[, lag_v[2] + 1], Lower = er_tp$matRRlow[, lag_v[2] + 1], Upper = er_tp$matRRhigh[, lag_v[2] + 1], row.names = NULL)
  o_er_tp <- data.frame(Variable= 'tp', Lag  = names(er_tp$allRRfit), RR  = as.numeric(t(er_tp$allRRfit)),
                        Lower  = as.numeric(t(er_tp$allRRlow)), Upper = as.numeric(t(er_tp$allRRhigh)))
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK FG10 value
  p_lr_fg10 <- get_visuals(lr_fg10, lag_rep=TRUE, col="#A6B727", title='fg10')
  lr_fg10 <- data.frame(Variable= 'fg10', Lag  = c(colnames(lr_fg10$matRRfit)), RR  = as.numeric(t(lr_fg10$matRRfit)),
                        Lower  = as.numeric(t(lr_fg10$matRRlow)), Upper = as.numeric(t(lr_fg10$matRRhigh)))
  # GRAPH CUMULATIVE FG10-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & T2M) 
  er_fg10 <- crosspred(cb.fg10, model, at=seq(0, 60, by=1), cen=cen_fg10)
  p_er_fg10 <- get_visuals(er_fg10, lag_rep=FALSE, lag=lag_v[3], col="#A6B727", xlab="Wind speed (km/hr)",  title='fg10')
  l_er_fg10 <- data.frame(Variable= 'fg10', Exposure = er_fg10$predvar, RR = er_fg10$matRRfit[, lag_v[3] + 1], Lower = er_fg10$matRRlow[, lag_v[3] + 1], Upper = er_fg10$matRRhigh[, lag_v[3] + 1], row.names = NULL)
  o_er_fg10 <- data.frame(Variable= 'fg10', Lag  = names(er_fg10$allRRfit), RR  = as.numeric(t(er_fg10$allRRfit)),
                          Lower  = as.numeric(t(er_fg10$allRRlow)), Upper = as.numeric(t(er_fg10$allRRhigh)))
  
  # Store values in tables
  lr <-rbind(lr_t2m, lr_tp, lr_fg10)
  er <- rbind(l_er_t2m, l_er_tp, l_er_fg10)
  oer <- rbind(o_er_t2m, o_er_tp, o_er_fg10)
  
  result <- list(model=model, overall=list(overall=df_cp, er=oer), 
                 exp_rep=list(all=er, t2m=p_er_t2m, tp=p_er_tp, fg10=p_er_fg10),
                 lag_rep=list(all=lr, t2m=p_lr_t2m, tp=p_lr_tp, fg10=p_lr_fg10))
  
  return(result)
  
}

# Function to get cumulative exposure-response and lag-response effects
get_visuals <- function(cp, xlab, title, col, lag_rep=FALSE, lag=NULL){
  
  
  if(lag_rep==FALSE){
    if(!is.null(lag)){
      cum <- data.frame(x = cp$predvar, RR = cp$matRRfit[, lag + 1], lower = cp$matRRlow[, lag + 1], upper = cp$matRRhigh[, lag + 1], row.names = NULL)
    } else {
      cum <- data.frame(x = cp$predvar, RR = cp$allRRfit, lower = cp$allRRlow, upper = cp$allRRhigh, row.names = NULL)
    }
    p <- ggplot(cum, aes(x = x, y = RR)) +
      geom_line(color = col) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = col) +
      geom_hline(yintercept = 1, color = "gray20") +
      scale_x_continuous(breaks = pretty(cum$x, n = 6)) +
      # scale_y_continuous(limits = c(0.8, 1.4)) +
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
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,  size  = 0.2) +
      geom_point(color = col,  size = 2.5) +
      scale_x_continuous(breaks = seq(min(slice$lag), max(slice$lag), by = 2)) +  # more breaks
      # scale_y_continuous(limits = c(0.8, 1.4)) +
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

# Function: Format plots for visualizing in grid
multi_plot_format <- function(plot_list, rows = 4, title_map = nm){
  
  p <- plot_list
  
  # Loop to replace titles using title_map if provided
  p <- lapply(names(p), function(n) {
    plot <- p[[n]]
    
    # Determine new title
    new_title <- if(!is.null(title_map) && n %in% names(title_map)) {
      title_map[n]
    } else {
      n  # fallback to the list name
    }
    
    plot + labs(title = new_title)
  })
  
  # Keep names after lapply
  names(p) <- names(plot_list)
  
  # Combine plots
  result <- wrap_plots(p, nrow = rows)
  
  return(result)
}

# Function to format and organizing visuals 
get_all_plots <- function(plot_list, single_outcome=NULL, rows=3){
  
  # Visualise all outputs
  if(!is.null(single_outcome)){
    
    # Arrange your 2x3 grid
    exp_rep <- (plot_list$exp_rep_t2m[[single_outcome]] | plot_list$exp_rep_tp[[single_outcome]] | plot_list$exp_rep_fg10[[single_outcome]]) + plot_annotation(title = paste0('Outcome: ', single_outcome), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep <- (plot_list$lag_rep_t2m[[single_outcome]] | plot_list$lag_rep_tp[[single_outcome]] | plot_list$lag_rep_fg10[[single_outcome]]) + plot_annotation(title = paste0('Outcome: ', single_outcome), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    all <- (exp_rep / lag_rep) + plot_annotation(title = paste0('Outcome: ', single_outcome), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    
    # Add a main title
    grid_plots <- list(all=all, exp_rep=exp_rep, lag_rep=lag_rep)
    
  } else {
    
    exp_rep_t2m <- multi_plot_format(plot_list$exp_rep_t2m) #+ plot_annotation(title = 'RR of t2m exposure-response (lag=5)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep_t2m <- multi_plot_format(plot_list$lag_rep_t2m) #+ plot_annotation(title = 'Lag-response of t2m (at 95th percentile)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    exp_rep_tp <- multi_plot_format(plot_list$exp_rep_tp) #+ plot_annotation(title = 'RR of tp exposure-response (lag=0)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep_tp <- multi_plot_format(plot_list$lag_rep_tp) #+ plot_annotation(title = 'Lag-response of tp (at 95th percentile)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    exp_rep_fg10 <- multi_plot_format(plot_list$exp_rep_fg10) #+ plot_annotation(title = 'RR of fg10 exposure-response (lag=0)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    lag_rep_fg10 <- multi_plot_format(plot_list$lag_rep_fg10) #+ plot_annotation(title = 'Lag-response of fg10 (at 95th percentile)', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
    
    grid_plots <- list(exp_rep_t2m=exp_rep_t2m, lag_rep_t2m=lag_rep_t2m, exp_rep_tp=exp_rep_tp, 
                       lag_rep_tp=lag_rep_tp, exp_rep_fg10=exp_rep_fg10, lag_rep_fg10=lag_rep_fg10)
    
  }
  
  return(grid_plots)
  
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/', file_import, '.RData'))

# Add lags to df
df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))

# Set primary outcome
outcome <- 'total_patients'


#### Get cumulative exposure-response (over all lags) ####

# Get results
results <- get_results(df, outcome)

# Get model plots
exp_rep_t2m <- results$exp_rep$t2m
exp_rep_tp <- results$exp_rep$tp
exp_rep_fg10 <- results$exp_rep$fg10

# Get model table
overall <- results$overall$overall
exp_rep_all <- results$overall$er

# Arrange plots
p_exp_rep_all <- (exp_rep_t2m | exp_rep_tp | exp_rep_fg10) #+ plot_annotation(title = 'Cumulative relative risk of weather exposure-response', theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))


#### Get the lag response (5th, 25th, 75th and 95th percentiles) ####

# Initialize lists
lag_rep_t2m <- list()
lag_rep_tp <- list()
lag_rep_fg10 <- list()
lag_rep_table <- data.frame()

quantiles <- c(0.05, 0.25, 0.75, 0.95)

for(q in quantiles){
  
  name <- paste0(as.character(q*100), 'th')
  # Get results
  results <- get_results(df, outcome, quant_v=q)
  
  # Get model plots
  lag_rep_t2m[[name]]<- results$lag_rep$t2m
  lag_rep_tp[[name]] <- results$lag_rep$tp
  lag_rep_fg10[[name]] <- results$lag_rep$fg10
  
  # Get model table
  lag_rep_all <- results$lag_rep$all
  lag_rep_all$q <- name
  lag_rep_table <- rbind(lag_rep_table, lag_rep_all)
  
}

# Arrange plots
p_lag_rep_t2m <- multi_plot_format(lag_rep_t2m, 1) + plot_annotation(title = 'Temperature', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_lag_rep_tp <- multi_plot_format(lag_rep_tp, 1) + plot_annotation(title = 'Rainfall', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_lag_rep_fg10 <- multi_plot_format(lag_rep_fg10, 1) + plot_annotation(title = 'Windspeed', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_lag_rep <- (p_lag_rep_t2m / p_lag_rep_tp / p_lag_rep_fg10) #+ plot_annotation(title = 'Lag-response of weather by percentile', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))


#### Get the exposure response (0, 5, 10, 14 lags) ####

# Initialize lists
exp_rep_t2m <- list()
exp_rep_tp <- list()
exp_rep_fg10 <- list()
exp_rep_table <- data.frame()

lags <- c(0, 5, 10, 14)

for(l in lags){
  
  name <- paste0('Lag ', as.character(l))
  # Get results
  results <- get_results(df, outcome, lag_v=c(l, l, l))
  
  # Get model plots
  exp_rep_t2m[[name]]<- results$exp_rep$t2m
  exp_rep_tp[[name]] <- results$exp_rep$tp
  exp_rep_fg10[[name]] <- results$exp_rep$fg10
  
  # Get model table
  exp_rep_all <- results$exp_rep$all
  exp_rep_all$l <- name
  exp_rep_table <- rbind(exp_rep_table, exp_rep_all)
  
}

# Arrange plots
p_exp_rep_t2m <- multi_plot_format(exp_rep_t2m, 1) + plot_annotation(title = 'Temperature', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_exp_rep_tp <- multi_plot_format(exp_rep_tp, 1) + plot_annotation(title = 'Rainfall', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_exp_rep_fg10 <- multi_plot_format(exp_rep_fg10, 1) + plot_annotation(title = 'Windspeed', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_exp_rep <- (p_exp_rep_t2m / p_exp_rep_tp / p_exp_rep_fg10) #+ plot_annotation(title = 'Exposure-response of weather by lag', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))


#### Check differences across responses ####

# Reinitialize dfs and lists
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
full_model <- list()
exp_rep_t2m <- list()
lag_rep_t2m <- list()
exp_rep_tp <- list()
lag_rep_tp <- list()
exp_rep_fg10 <- list()
lag_rep_fg10 <- list()
exp_rep_table <- data.frame()

# Extract data for all outcomes
outcomes <- c('total_patients', 'F', 'M', 'under_5', 'age5_14', 'age15_24', 'age25_64', 'over_65', 
              'general', 'ophthalmology', 'pediatrics', 'orthopedics', 'ent', 'gynecology', 'dermatology')

for (i in outcomes) {
  
  results <- get_results(df, i, quant_v=0.75, lag_v=c(5, 0, 0))
  
  # Get model fit
  test_result <- get_model_fit(results$model, i)
  modelfitstats <- rbind(modelfitstats, test_result)
  
  # Get model
  results$model
  
  # Get model plots
  exp_rep_t2m[[i]] <- results$exp_rep$t2m
  lag_rep_t2m[[i]] <- results$lag_rep$t2m
  exp_rep_tp[[i]] <- results$exp_rep$tp
  lag_rep_tp[[i]] <- results$lag_rep$tp
  exp_rep_fg10[[i]] <- results$exp_rep$fg10
  lag_rep_fg10[[i]] <- results$lag_rep$fg10
  
  # Get model table
  exp_rep_all <- results$exp_rep$all
  exp_rep_all$l <- i
  exp_rep_table <- rbind(exp_rep_table, exp_rep_all)
  
}

plots <- list(exp_rep_t2m=exp_rep_t2m, lag_rep_t2m=lag_rep_t2m, exp_rep_tp=exp_rep_tp, 
              lag_rep_tp=lag_rep_tp, exp_rep_fg10=exp_rep_fg10, lag_rep_fg10=lag_rep_fg10)

all_results <- get_all_plots(plots)

# Check specific results
unique(exp_rep_table$l)
filter <- exp_rep_table %>% filter(Variable=='tp',
                                   l== "dermatology")

# Get single outcome from this
# total_patients <- get_all_plots(plots, single_outcome = 'total_patients')

#### Export outcome plots #####
width <- 9
height <- 6

# Save cumulative exposure-response (total_patient)
ggsave(filename = paste0(path3, '/IV. Figures/MF_total_exp_rep_all.png'), plot = p_exp_rep_all,
       width=7, height=4, dpi = 300)

# Save lag-response by percentile (total_patient)
ggsave(filename = paste0(path3, '/IV. Figures/MF_total_lag_rep_perc.png'), plot = p_lag_rep,
       width=width, height=height, dpi = 300)

# Save exposure-response by lag (total_patient)
ggsave(filename = paste0(path3, '/IV. Figures/MF_total_exp_rep_lags.png'), plot = p_exp_rep,
       width=width, height=height, dpi = 300)

# Population plots
names(all_results) <- paste0("disag_", names(all_results))
for (nm in names(all_results)) {
  ggsave(filename = paste0(path3, '/IV. Figures/MF_', nm, '.png'), plot = all_results[[nm]],
         width=10, height=8, dpi = 300)
}


# # Single plots
# for (nm in names(total_patients)) {
#   ggsave(filename = paste0(path3, '/IV. Figures/total_patients_', nm, '.png'), plot = total_patients[[nm]],
#          width=width, height=height, dpi = 300)
# }


