### ------------------------------------------------------------------------ ### 
### HISP DLNM - Dependencies                                                 ###
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/

# Turn off scientific notation
options(scipen=999)

#### PATHS ####
user <- 'claer14'
### Project
path_data <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
path_ch2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 2 - Causal Inference/')
### Git
path_git <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path_project <- paste0(path_git, '2_DLNM/')
path_report <- paste0(path_project, '- Reporting/')

#### PACKAGES ####

# List packages
packages <- c('dplyr', 'tidyr', 'lubridate', 'zoo', 'dlnm', 'splines', 
              'ggplot2', 'patchwork')

# library(MASS); library(ggfortify); library(gridExtra)

# Install if needed
installed <- packages %in% rownames(installed.packages())
if(any(!installed)) {
  install.packages(packages[!installed])
}

# Load packages
lapply(packages, library, character.only = TRUE)

#### FUNCTIONS ####

# Function to replace exact matches from a list
replace_exact_matches <- function(df, dict) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      col <- ifelse(col %in% names(dict), dict[col], col)
    }
    col
  })
  df
}

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
get_results <- function(df, outcome, quant_v=0.95, lag_v=NULL){
  
  # Fit cross-basis
  cb.t2m <- crossbasis(df$t2m, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=4), group=df$year)
  cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=4), group=df$year)
  # cb.t2m <- crossbasis(df$mx2t, lag=c(0,14), argvar=list(fun='ns', df=3), arglag=list(fun='ns', df=8))
  # cb.tp <- crossbasis(df$tp, lag=c(0,14),argvar=list(fun='ns', df=4), arglag=list(fun='ns', df=4))
  # cb.fg10 <- crossbasis(df$fg10, lag=c(0,14),argvar=list(fun='ns', df=5), arglag=list(fun='ns', df=5))
  
  # Fit model
  model <- glm(formula=as.formula(paste0(outcome, ' ~  season + weekday + holiday + ns(time, df=9*9) + cb.t2m + cb.tp + cb.fg10 + patients_lag')), data = df, family = quasipoisson(link="log"))
  
  # Select quantile
  q_t2m <- quantile(df$t2m, probs = quant_v, na.rm = TRUE)
  q_tp <- quantile(df$tp, probs = quant_v, na.rm = TRUE)
  q_fg10 <- quantile(df$fg10, probs = quant_v, na.rm = TRUE)
  
  # DISPLAY ESTIMATE OF RR FOR BENCHMARK TP INCREMENT OF 23.4 (99th perc) 
  cp.t2m <- crosspred(cb.t2m, model, at=q_t2m)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.tp <- crosspred(cb.tp, model, at=q_tp)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp.fg10 <- crosspred(cb.fg10, model, at=q_fg10)[c('predvar', 'allRRfit','allRRlow','allRRhigh')]
  cp <- data.frame(rbind(cp.t2m, cp.tp, cp.fg10))
  df_cp <- data.frame(variables = rownames(cp), cp, row.names = NULL)
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK T2M value
  cp <- crosspred(cb.t2m, model, at=q_t2m)
  p_lr_t2m <- get_visuals(cp, model, lag_rep=TRUE, col="#D74B3C", title='Lag-response of t2m (at 95th percentile)')
  # GRAPH CUMULATIVE T2M-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & FG10) 
  cp <- crosspred(cb.t2m, model, at=seq(-10, 35, by=1))
  p_er_t2m <- get_visuals(cp, model, lag_rep=FALSE, lag=lag_v, col="#D74B3C", xlab="Temperature (C)", title="RR of t2m exposure-response (lag=5)")
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK TP value
  cp <- crosspred(cb.tp, model, at=q_tp)
  p_lr_tp <- get_visuals(cp, model, lag_rep=TRUE, col="#1D98EB", title='Lag-response of tp (at 95th percentile)')
  # GRAPH CUMULATIVE TP-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR T2M & FG10) 
  cp <- crosspred(cb.tp, model)
  p_er_tp <- get_visuals(cp, model, lag_rep=FALSE, lag=lag_v, col="#1D98EB", xlab="Rainfall (mm)", title="RR of tp exposure-response (lag=6)")
  
  # PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK FG10 value
  cp <- crosspred(cb.fg10, model, at=q_fg10)
  p_lr_fg10 <- get_visuals(cp, model, lag_rep=TRUE, col="#A6B727", title='Lag-response of fg10 (at 95th percentile)')
  # GRAPH CUMULATIVE FG10-ADMISSIONS CUMULTIVE ASSOCIATION (AFTER ADJUSTING FOR TP & T2M) 
  cp <- crosspred(cb.fg10, model)
  p_er_fg10 <- get_visuals(cp, model, lag_rep=FALSE, lag=lag_v, col="#A6B727", xlab="Wind speed (km/hr)", title="Overall RR of fg10 exposure-response")
  
  result <- list(model=model, overall=df_cp, 
                 exp_rep = list(t2m=p_er_t2m, tp=p_er_tp, fg10=p_er_fg10),
                 lag_rep = list(t2m=p_lr_t2m, tp=p_lr_tp, fg10=p_lr_fg10))
  
  return(result)
  
}

# Function to get cumulative exposure-response and lag-response effects
get_visuals <- function(cp, model, xlab, title, col, lag_rep=FALSE, lag=NULL){
  
  
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
      # scale_y_continuous(limits = c(-0.5, 2.5)) +
      labs(x = xlab, y = "RR",
           title = title) +
      theme_classic() +
      theme(
        axis.line = element_line(color = "gray20"),      # pale axis lines
        axis.text = element_text(color = "gray20"),      # pale axis labels
        axis.ticks = element_line(color = "gray20"),     # pale ticks
        plot.title = element_text(face = "bold", hjust = 0.5, size=12)  # center title
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
        plot.title = element_text(face = "bold", hjust = 0.5, size=12)  # center title
      )
    
  }
  
  return(p)
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

#### DEPENDENCIES ####

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD','#52C5B0FF', '#910030FF', '#127804FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9','#52C5B080', '#91003080', '#12780480') 

w_cols <- c('t2m', 'mn2t','mx2t','tp','fg10', 'rh_mean', 'sh_mean', 'pm25')
h_cols <- c('total_patients',
            'F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
            'gynaecology', 'orthopedics', 'dermatology', 'surgery',
            'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
            'ncd', 'antiretroviral', 'substance_counselling', 'physiotherapy')

location_levels <- c('Bilaspur','Chamba', 'Dharamshala', 'Hamirpur', 'Kullu', 'Kyelang', 
                     'Mandi', 'Nahan', 'Reckong Peo', 'Shimla', 'Solan', 'Una')
season_levels   <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
month_levels    <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
weekday_levels  <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
age_levels      <- c('under_5', 'age5_14', 'age15_24', 'age25_64',  'over_65','Unknown')
gender_levels   <- c('F', 'M', 'O')
ward_levels     <- c('general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
                     'gynaecology', 'orthopedics', 'dermatology', 'surgery',
                     'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
                     'ncd', 'mch', 'antiretroviral', 'substance_counselling', 'physiotherapy')

# Dictionary of names for figures and tables
nm <- c(season='Season', closed='Closed', holiday='Holiday', total_patients='Total patients', 
        # POPULATION GROUPS
        M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        # WARDS
        general='General', emergency='Emergency', ophthalmology='Ophthalmology', 
        paediatrics='Paediatric', orthopedics='Orthopedics', ent='ENT', 
        gynaecology='Gynaecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        pulmonary='Pulmonary', disability='Disability', ncd='NCD', mch='MCH',
        antiretroviral='HIV Treatment', substance_counselling='Substance counselling',
        physiotherapy='Physiotherapy', 
        # CONTINUOUS WEATHER
        t2m='Mean temp. (\u00B0C)', mx2t='Max temp. (\u00B0C)', mn2t='Min temp. (\u00B0C)',
        rh_mean='Relative humidity (%)', tp='Rainfall (mm)', fg10='Windspeed (km/hr)', pm25='Air pollution (µg/m³)',
        # IMD DEFINITIONS
        mx2t_IMD='IMD: Mx temp. >30\u00B0C', mn2t_IMD='IMD: Min temp. <0\u00B0C', 
        tp_IMD='IMD: Rainfall >15mm', fg10_IMD='IMD: Wind speed >40km/h', pm25_IMD='IMD: PM2.5 >60 µg/m³', 
        rh_low_IMD='IMD: Humidity <30%', rh_hig_IMD='IMD: Humidity >60%',
        # CLIM. DEFINITIONS
        mx2t_p90='Clim: Max temp. >90th', mn2t_p10='Clim: Min temp. <10th', 
        tp_p90='Clim: Rainfall >90th', fg10_p90='Clim: Wind speed >90th', pm25_p90='Clim: PM2.5 >90th', 
        rh_mean_p10 ='Clim: Humidity <10th', rh_mean_p90='Clim: Humidity >90th') 

