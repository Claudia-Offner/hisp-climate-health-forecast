### ------------------------------------------------------------------------ ### 
###  DNLM - Statistical Analysis                                             ###
### ------------------------------------------------------------------------ ### 

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# Set data path
setwd(path_git)

### FUNCTIONS ####
# Function to format plots for visualizing in grid
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


#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Add lags to df
df <- df %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1)) %>%
  filter(location=='Reckong Peo')

# Set primary outcome
outcome <- 'total_patients'

# 1. Get cumulative exposure-response (over all lags) ####

# Get results
results <- get_results(df, outcome)

# Get model plots
exp_rep_wb <- results$exp_rep$wb
exp_rep_tp <- results$exp_rep$tp
exp_rep_rh <- results$exp_rep$rh

# Get model table
overall <- results$overall$overall
exp_rep_all <- results$overall$er

# Arrange plots
(p_exp_rep_all <- (exp_rep_wb | exp_rep_tp | exp_rep_rh) + plot_annotation(title = 'Cumulative relative risk of weather exposure-response', theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))))


# 2. Get the lag response (5th, 25th, 75th and 95th percentiles) ####

# Initialize lists
lag_rep_wb <- list()
lag_rep_tp <- list()
lag_rep_rh <- list()
lag_rep_table <- data.frame()

quantiles <- c(0.05, 0.25, 0.75, 0.95)

for(q in quantiles){
  
  name <- paste0(as.character(q*100), 'th')
  # Get results
  results <- get_results(df, outcome, quant_v=q)
  
  # Get model plots
  lag_rep_wb[[name]]<- results$lag_rep$wb
  lag_rep_tp[[name]] <- results$lag_rep$tp
  lag_rep_rh[[name]] <- results$lag_rep$rh
  
  # Get model table
  lag_rep_all <- results$lag_rep$all
  lag_rep_all$q <- name
  lag_rep_table <- rbind(lag_rep_table, lag_rep_all)
  
}

# Arrange plots
p_lag_rep_wb <- multi_plot_format(lag_rep_wb, 1) + plot_annotation(title = 'Temperature', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_lag_rep_tp <- multi_plot_format(lag_rep_tp, 1) + plot_annotation(title = 'Rainfall', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_lag_rep_rh <- multi_plot_format(lag_rep_rh, 1) + plot_annotation(title = 'Windspeed', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
(p_lag_rep <- (p_lag_rep_wb / p_lag_rep_tp / p_lag_rep_rh) + plot_annotation(title = 'Lag-response of weather by percentile', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))


# 3. Get the exposure response (0, 5, 10, 14 lags) ####

# Initialize lists
exp_rep_wb <- list()
exp_rep_tp <- list()
exp_rep_rh <- list()
exp_rep_table <- data.frame()

lags <- c(0, 5, 10, 14)

for(l in lags){
  
  name <- paste0('Lag ', as.character(l))
  
  # Get results
  results <- get_results(df, outcome, lag_v=l)
  
  # Get model plots
  exp_rep_wb[[name]]<- results$exp_rep$wb
  exp_rep_tp[[name]] <- results$exp_rep$tp
  exp_rep_rh[[name]] <- results$exp_rep$rh
  
  # Get model table
  exp_rep_all <- results$exp_rep$all
  exp_rep_all$l <- name
  exp_rep_table <- rbind(exp_rep_table, exp_rep_all)
  
}

# Arrange plots
p_exp_rep_wb <- multi_plot_format(exp_rep_wb, 1) + plot_annotation(title = 'Temperature', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_exp_rep_tp <- multi_plot_format(exp_rep_tp, 1) + plot_annotation(title = 'Rainfall', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
p_exp_rep_rh <- multi_plot_format(exp_rep_rh, 1) + plot_annotation(title = 'Windspeed', theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
(p_exp_rep <- (p_exp_rep_wb / p_exp_rep_tp / p_exp_rep_rh) + plot_annotation(title = 'Exposure-response of weather by lag', theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))


# # 4. Check differences across responses ####
# 
# # Reinitialize dfs and lists
# setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c('model','obs','dfmodel','dispersion','deviance','AIC' ))
# full_model <- list()
# exp_rep_wb <- list()
# lag_rep_wb <- list()
# exp_rep_tp <- list()
# lag_rep_tp <- list()
# exp_rep_rh <- list()
# lag_rep_rh <- list()
# exp_rep_table <- data.frame()
# 
# # Extract data for all outcomes
# outcomes <- c('total_patients', 'F', 'M', 'under_5', 'age5_14', 'age15_24', 'age25_64', 'over_65',
#               'general', 'ophthalmology', 'paediatrics', 'orthopedics', 'ent', 'gynaecology', 'dermatology')
# # outcomes <- c('F', 'M', 'under_5', 'age5_14', 'age15_24', 'age25_64', 'over_65', 
# #               'pediatrics', 'orthopedics', 'ent', 'gynecology', 'dermatology') # FOR BAS POSTER
# 
# for (i in outcomes) {
#   
#   results <- get_results(df, i, quant_v=0.75, lag_v=c(5, 0, 0))
#   
#   # Get model fit
#   test_result <- get_model_fit(results$model, i)
#   modelfitstats <- rbind(modelfitstats, test_result)
#   
#   # Get model
#   results$model
#   
#   # Get model plots
#   exp_rep_wb[[i]] <- results$exp_rep$wb
#   lag_rep_wb[[i]] <- results$lag_rep$wb
#   exp_rep_tp[[i]] <- results$exp_rep$tp
#   lag_rep_tp[[i]] <- results$lag_rep$tp
#   exp_rep_rh[[i]] <- results$exp_rep$rh
#   lag_rep_rh[[i]] <- results$lag_rep$rh
#   
#   # Get model table
#   exp_rep_all <- results$exp_rep$all
#   exp_rep_all$l <- i
#   exp_rep_table <- rbind(exp_rep_table, exp_rep_all)
#   
# }
# 
# plots <- list(exp_rep_wb=exp_rep_wb, lag_rep_wb=lag_rep_wb, exp_rep_tp=exp_rep_tp, 
#               lag_rep_tp=lag_rep_tp, exp_rep_rh=exp_rep_rh, lag_rep_rh=lag_rep_rh)
# 
# all_results <- get_all_plots(plots)
# 
# # Check specific results
# unique(exp_rep_table$l)
# filter <- exp_rep_table %>% filter(Variable=='tp',
#                                    l== "dermatology")
# 
# # Get single outcome from this
# # total_patients <- get_all_plots(plots, single_outcome = 'total_patients')
# 
# #### EXPORT #####
# width <- 9
# height <- 6
# 
# # Save cumulative exposure-response (total_patient)
# ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_total_exp_rep_all.png'), plot = p_exp_rep_all, width=7, height=4, dpi = 300)
# ggsave(filename = paste0(path_report, '- figures/MF_total_exp_rep_all.png'), plot = p_exp_rep_all, width=7, height=4, dpi = 300)
# 
# # Save lag-response by percentile (total_patient)
# ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_total_lag_rep_perc.png'), plot = p_lag_rep, width=width, height=height, dpi = 300)
# ggsave(filename = paste0(path_report, '- figures/MF_total_lag_rep_perc.png'), plot = p_lag_rep, width=width, height=height, dpi = 300)
# 
# # Save exposure-response by lag (total_patient)
# ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_total_exp_rep_lags.png'), plot = p_exp_rep, width=width, height=height, dpi = 300)
# ggsave(filename = paste0(path_report, '- figures/MF_total_exp_rep_lags.png'), plot = p_exp_rep, width=width, height=height, dpi = 300)
# 
# # Population plots
# names(all_results) <- paste0("disag_", names(all_results))
# for (nm in names(all_results)) {
#   ggsave(filename = paste0(path_ch2, '/IV. Figures/MF_', nm, '.png'), plot = all_results[[nm]], width=10, height=8, dpi = 300)
#   ggsave(filename = paste0(path_report, '- figures/MF_', nm, '.png'), plot = all_results[[nm]], width=12, height=7, dpi = 300)
# }
# 
# 
# # # Single plots
# # for (nm in names(total_patients)) {
# #   ggsave(filename = paste0(path_ch2, '/IV. Figures/total_patients_', nm, '.png'), plot = total_patients[[nm]],
# #          width=width, height=height, dpi = 300)
# # }
# 
# 
