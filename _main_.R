### ------------------------------------------------------------------------ ### 
###  DLNM - Main script
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
options(scipen=999); set.seed(123)

# Set data path
user <- 'claer14'
path_data <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/', user, '/Documents/GitHub/hisp-climate-health-forecast/')
path_report <- paste0(path_git, 'results/')
setwd(path_git)
start <- Sys.time()

#### DEPENDENCIES ####

# Install and load packages
pkgs <- c('xgboost', 'dplyr', 'data.table', 'caret', 'purrr', 'gridExtra', 
          'scales','ggplotify', 'dlnm', 'splines', 'xtable', 'mixmeta', 'MASS', 
          'ggplot2', 'patchwork', 'TTR', 'openxlsx', 'colorspace', 'grid',
          'reshape2', 'gtable', 'cowplot', 'lmtest', 'flextable', 'officer', 
          'kableExtra', 'moments', 'stringr', 'tidyr', 'lubridate', 'zoo', 
          'psychrolib', 'ggpattern', 'ggnewscale', 'conflicted')
install.packages(pkgs[!pkgs %in% installed.packages()[, 'Package']])
invisible(lapply(pkgs, library, character.only = TRUE))
source("https://raw.githubusercontent.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/master/attrdl.R")
conflicted::conflicts_prefer(
  kableExtra::footnote,
  dplyr::filter,
  dplyr::select
)

# Function to identify appropriate reference for tp after log transformation
get_tp_ref <- function(df, locations) {
  cen_tp <- sapply(split(df, df$location)[locations], function(x) quantile(x$tp, (0:100)/100, na.rm = TRUE))
  first_valid <- which(apply(cen_tp, 1, function(row) all(abs(row) > 1e-10) && !anyDuplicated(round(row, 10))))[1]
  paste0(ceiling(as.numeric(sub('%', '', rownames(cen_tp)[first_valid])) / 5) * 5, '%')
}

####  MAIN CODE  #####

# 0. Load data and set parameters ####

# Load data
load(paste0('ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData'))
data <- df_dlnm; rm(df_dlnm)

# Set location info
locations <- data |> distinct(location, elev) |> arrange(elev) |> pull(location) |> as.character()
m <- length(locations)

# Set primary outcome
outcome <- 'total_patients'

# Keep variables
keep_vars <- ls()

# # !!! Percentile checks !!!
# df <- data
# df$var <- data$tp_roll7
# # mean(df$var, na.rm = TRUE) - 6.4
# cen <- sapply(split(df, df$location)[locations], function(x) quantile(x$var, (0:100)/100, na.rm = TRUE))
# cen <- cbind(cen, Pooled = rowMeans(cen))

# 1. Parameter selection ####

analyses <- list(
  t2m       = list(data       = data,
                   immediate  = list(tmp = 't2m',         tp = 'tp_log'),
                   cumulative = list(tmp = 't2m_roll7',   tp = 'tp_roll7_log'))
  # mn2t      = list(data       = data,
  #                  immediate  = list(tmp = 'mn2t',        tp = 'tp_log'),
  #                  cumulative = list(tmp = 'mn2t_roll7',  tp = 'tp_roll7_log')),
  # mx2t      = list(data       = data,
  #                  immediate  = list(tmp = 'mx2t',        tp = 'tp_log'),
  #                  cumulative = list(tmp = 'mx2t_roll7',  tp = 'tp_roll7_log'))
)

s <- '1_parameter_selection.R'
for (a in names(analyses)) {
  # Load model specifications and data
  model_specs <- analyses[[a]]
  df <- analyses[[a]]$data
  for (model in names(model_specs[-1])) {
    # Load parameters
    spec <- model_specs[[model]]
    df$tmp <- df[[spec$tmp]]
    df$tp   <- df[[spec$tp]]
    tp_ref  <- get_tp_ref(df, locations)
    source('_utilities_.R')
    source(s, echo = TRUE)
    rm(list = setdiff(ls(), c(keep_vars, 'keep_vars', 's', 'analyses', 'a', 'df', 'model_specs', 'model')))
  }
}

# 2. Run models ####

analyses <- list(
  t2m       = list(data       = data,
                   immediate  = list(tmp = 't2m',         tp = 'tp_log',        file = 't2m_immediate_final_params.RData'),
                   cumulative = list(tmp = 't2m_roll7',   tp = 'tp_roll7_log',  file = 't2m_cumulative_final_params.RData')),
  roll2  = list(data=data, cumulative = list(tmp = 't2m_roll2',  tp = 'tp_roll2_log',  file = 't2m_cumulative_final_params.RData')),
  roll3  = list(data=data, cumulative = list(tmp = 't2m_roll3',  tp = 'tp_roll3_log',  file = 't2m_cumulative_final_params.RData')),
  roll4  = list(data=data, cumulative = list(tmp = 't2m_roll4',  tp = 'tp_roll4_log',  file = 't2m_cumulative_final_params.RData')),
  roll5  = list(data=data, cumulative = list(tmp = 't2m_roll5',  tp = 'tp_roll5_log',  file = 't2m_cumulative_final_params.RData')),
  roll6  = list(data=data, cumulative = list(tmp = 't2m_roll6',  tp = 'tp_roll6_log',  file = 't2m_cumulative_final_params.RData')),
  mx2t      = list(data       = data,
                   immediate  = list(tmp = 'mx2t',        tp = 'tp_log',        file = 't2m_immediate_final_params.RData'),
                   cumulative = list(tmp = 'mx2t_roll7',  tp = 'tp_roll7_log',  file = 't2m_cumulative_final_params.RData')),
  mn2t      = list(data       = data,
                   immediate  = list(tmp = 'mn2t',        tp = 'tp_log',        file = 't2m_immediate_final_params.RData'),
                   cumulative = list(tmp = 'mn2t_roll7',  tp = 'tp_roll7_log',  file = 't2m_cumulative_final_params.RData'))
)

s <- '2_statistical_analysis'
for (a in names(analyses)) {
  # Load model specifications and data
  model_specs <- analyses[[a]]
  df <- analyses[[a]]$data
  for (model in names(model_specs[-1])) {
    # Load parameters
    spec <- model_specs[[model]]
    if (is.character(spec$file)) {
      load(paste0(path_project, 'Results/', spec$file))
      final_params$pooled$shape <- 'binom'
    } else{
      final_params <- list()
      final_params$pooled <- spec$file
    }
    df$tmp <- df[[spec$tmp]]
    df$tp   <- df[[spec$tp]]
    tp_ref  <- get_tp_ref(df, locations)
    source('_utilities_.R')
    source(s, echo = TRUE)
  }
}

rm(list = setdiff(ls(), c(keep_vars, 'keep_vars', 's', 'analyses', 'a', 'df', 'model_specs', 'model')))


# 3. Export figures and tables ####

s <- '3_figures_and_tables'
for (a in names(analyses)) {
  model_specs <- analyses[[a]]
  model1 <- names(model_specs)[2]; model2 <- names(model_specs)[3]
  if (any(grepl('immediate', names(model_specs)))) {
    folder1 <- model1; folder2 <- model2
  } else {
    folder1 <- a; folder2 <- a
  }
  df <- analyses[[a]]$data
  tp_ref  <- get_tp_ref(df, locations)
  source('_utilities_.R')
  source(s, echo = TRUE)
}

rm(list = setdiff(ls(), c(keep_vars, 'keep_vars', 's', 'analyses', 'a', 'model_specs', 'model')))


#### End timer: 2 hours 35 minutes 40 seconds (31/03/2026) ####
end <- Sys.time()
dt <- round(as.numeric(end - start, units = 'secs'))
hours   <- dt %/% 3600
minutes <- (dt %% 3600) %/% 60
seconds <- dt %% 60
sprintf('%d hours %d minutes %.d seconds', hours, minutes, seconds)

