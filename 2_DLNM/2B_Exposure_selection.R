### ------------------------------------------------------------------------ ### 
###  DNLM: TEST XGBoost algorithm for variable selection
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls())
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')
source('2_DLNM/DLNM - Dependencies.R')

# Load packages
library(xgboost); library(dplyr); library(data.table); library(caret); library(purrr)

#### FUNCTIONS ####

# Function to run XG Boost test on subset
xgboost <- function(data, title, var_list, params=NULL) {
  
  # 1. Prepare data ----
  data <- data %>% filter(!is.na(total_patients)) %>% arrange(date) %>%
    mutate(location = as.numeric(as.factor(location)))
  
  # Choose features
  features <- names(data)[grepl(var_list, names(data))]
  X <- as.matrix(data[, features] %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))))
  y <- data$total_patients
  
  # 2. Split into training and validation sets ----
  train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_valid <- X[-train_idx, ]
  y_valid <- y[-train_idx]
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dvalid <- xgb.DMatrix(data = X_valid, label = y_valid)
  watchlist <- list(train = dtrain, val = dvalid)
  
  # 3. Hyper parameter optimization ----
  
  if(is.null(params)){
    param_grid <- expand.grid(
      eta = c(0.01, 0.05, 0.1),
      max_depth = c(4, 6, 8),
      subsample = c(0.7, 0.8, 1),
      colsample_bytree = c(0.7, 0.8, 1)
    )
    
    results <- data.frame()
    
    for(i in 1:nrow(param_grid)){
      params <- list(
        objective = 'count:poisson',
        eval_metric = 'rmse',
        eta = param_grid$eta[i],
        max_depth = param_grid$max_depth[i],
        subsample = param_grid$subsample[i],
        colsample_bytree = param_grid$colsample_bytree[i]
      )
      
      model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = 500,               # high number, early stopping decides optimal
        watchlist = watchlist,
        early_stopping_rounds = 50,
        verbose = 0
      )
      
      results <- rbind(
        results,
        data.frame(
          eta = params$eta,
          max_depth = params$max_depth,
          subsample = params$subsample,
          colsample_bytree = params$colsample_bytree,
          best_iteration = as.numeric(xgb.attributes(model)$best_iteration),
          val_rmse = as.numeric(xgb.attributes(model)$best_score)
        )
      )
    }
    
    best_params <- results %>% arrange(val_rmse) %>% slice(1)
    print(best_params)
    
  } else {
    
    best_params <- params
    
  }
  print(best_params)
  
  # 4. Train final model on full training set ----
  final_params <- list(
    objective = 'count:poisson',
    eval_metric = 'rmse',
    eta = best_params$eta,
    max_depth = best_params$max_depth,
    subsample = best_params$subsample,
    colsample_bytree = best_params$colsample_bytree
  )
  best_iter <- best_params$best_iteration
  
  xgb_model <- xgb.train(
    params = final_params,
    data = dtrain,
    nrounds = best_iter,
    watchlist = watchlist,
    verbose = 1
  )
  
  # 5. Feature importance ----
  importance <- xgb.importance(model = xgb_model)
  importance <- replace_exact_matches(importance, nm)
  plot <- ggplot(importance, aes(x = Gain, y = reorder(Feature, Gain))) +
      geom_col(fill = "#1D98EB") +
      labs(x = "Gain", y = 'Exposure', title = title) +
      theme_classic() +
      theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(plot)
  
  return(list(best_params=best_params, importance=importance, plot=plot))
  
}

#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# 1. Non parameter correlation test ####

sub <- df %>% ungroup() %>% dplyr::select(total_patients, all_of(w_cols))
cor_matrix <- cor(sub, use = 'complete.obs', method = 'spearman')
melted_corr <- reshape2::melt(cor_matrix)
(cor_matrix <- ggplot(data = melted_corr, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low = '#1D98EB', high = '#D74B3C', mid = 'white', midpoint = 0, limit = c(-1,1)) +
    theme_minimal() +
    scale_x_discrete(labels = function(x) {
      # If x is in nm, replace with nm value, else keep original
      sapply(x, function(lbl) {
        if (lbl %in% names(nm)) {
          nm[[lbl]]
        } else {
          lbl
        }})}) +
    scale_y_discrete(labels = function(x) {
      # If x is in nm, replace with nm value, else keep original
      sapply(x, function(lbl) {
        if (lbl %in% names(nm)) {
          nm[[lbl]]
        } else {
          lbl
        }})}) +
    geom_text(aes(label=round(value,2)), color='black', size=4) +
    labs(title = 'Spearman Correlation Matrix', x='', y='') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # rotate x-axis labels
      axis.text.y = element_text(angle = 0),
      plot.title = ggplot2::element_text(hjust = 0.5)
      ))

# 2A. Lags: XG Boost feature importance test ####

# # Select parameters
# control_cols <- c('location', 'time', 'season', 'weekday', 'elev', 'popu')
# lagged_cols <- grep("_lag[1-7]$", names(df), value = TRUE)
# var_list <- paste(c(control_cols, w_cols, lagged_cols), collapse = "|")
# 
# # Get best parameters & results
# start_time <- Sys.time()
# xg_res <- list()
# xg_res[['all']] <- xgboost(df, "XGBoost Feature Importance", var_list, params=NULL)
# for (l in locations){
#   sub <- df %>% filter(location==l)
#   xg_res[[l]] <- xgboost(sub, l, var_list, params=NULL)
# }
# end_time <- Sys.time()
# # Print time & save (~80min)
# print(difftime(end_time, start_time, units = "mins"))
# save(xg_res, file = paste0(path_project, 'ind_xgboost_params_lags.RData'))
# 

# 2B. No Lags: XG Boost feature importance test ####

# # Select parameters
# var_list <- paste(c(control_cols, w_cols), collapse = "|")
# 
# # Get best parameters & results
# start_time <- Sys.time()
# xg_res <- list()
# xg_res[['all']] <- xgboost(df, "XGBoost Feature Importance", var_list, params=NULL)
# for (l in locations){
#   sub <- df %>% filter(location==l)
#   xg_res[[l]] <- xgboost(sub, l, var_list, params=NULL)
# }
# end_time <- Sys.time()
# # Print time & save (~35min)
# print(difftime(end_time, start_time, units = "mins"))
# save(xg_res, file = paste0(path_project, 'ind_xgboost_params.RData'))


# 3A. Lags: Format results ####
# Get results
load('2_DLNM/ind_xgboost_params_lags.RData')

# Assessment table of feature importance consistency 
imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
# Assessment table of feature importance consistency 
imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
(imp_rank_lag <- imp %>%
    group_by(subpop) %>%
    arrange(desc(Gain)) %>%
    mutate(
      top    = row_number() <= 20,
      bottom = row_number() > (n() - 120)
    ) %>%
    ungroup() %>%
    group_by(Feature) %>%
    summarise(
      `Top 20`    = sum(top),
      `Bottom 120` = sum(bottom),
      `Mean Gain`     = mean(Gain),
      .groups = "drop"
    ) %>% 
    mutate(Important= `Bottom 120` < `Top 20`) %>%
    filter(Important==TRUE) %>%
    arrange(desc(Important), desc(`Mean Gain`))) 
imp_rank_lag <- replace_exact_matches(imp_rank_lag, nm)

# 3B. No Lags: Format results ####
# Get results
load('2_DLNM/ind_xgboost_params.RData')

# Plot all feature importance
imp_plot <- xg_res$all$plot

# Plot location feature importance
imp_loc <- Map(function(p, nm) {
  # Base theme: remove y-axis title, tight margins
  base_theme <- theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )
  
  # Conditionally remove x-axis title for everything except Solan and Una
  x_theme <- ifelse(nm %in% c("Solan", "Una"),
                    list(theme()),  # keep x-axis title
                    list(theme(axis.title.x = element_blank())))  # remove it
  
  # Combine base + conditional theme
  p + base_theme + x_theme[[1]]
  
}, map(xg_res[-1], 'plot'), names(xg_res[-1]))
imp_loc_plot <- wrap_plots(imp_loc, nrow=2)

# Assessment table of feature importance consistency 
imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
# Assessment table of feature importance consistency 
imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
(imp_rank <- imp %>%
    group_by(subpop) %>%
    arrange(desc(Gain)) %>%
    mutate(
      top    = row_number() <= 7,
      bottom = row_number() > (n() - 7)
    ) %>%
    ungroup() %>%
    group_by(Feature) %>%
    summarise(
      `Top Half`    = sum(top),
      `Bottom Half` = sum(bottom),
      `Mean Gain`     = mean(Gain),
      .groups = "drop"
    ) %>% 
    mutate(Important= `Bottom Half` < `Top Half`) %>%
    arrange(desc(Important), desc(`Mean Gain`))) 
imp_rank <- replace_exact_matches(imp_rank, nm)

#### EXPORT ####

# # Figures
# export <- list('SF_corr_plot'=list(p=cor_matrix, w=7, h=5),
#                'SF_xgb_importance'=list(p=imp_plot, w=7, h=5),
#                'SF_xgb_importance_loc'=list(p=imp_loc_plot, w=7, h=5))
# for (name in names(export)){
#   e <- export[[name]]
#   ggsave(paste0(path_report, 'figures/sensitivity', name, '.png'), e$p, width = e$w, height = e$h, dpi = 300, bg='white')
# }
# 
# # Table
# write.csv(imp_rank, file = paste0(path_report, 'tables/ST_imp_rank.csv'), row.names = FALSE)
# write.csv(imp_rank_lag, file = paste0(path_report, 'tables/ST_imp_rank_lags.csv'), row.names = FALSE)

# # Interaction strength ----
# 
# shap_interactions <- predict(xgb_model, X, predinteraction = TRUE)
# shap_importance <- colMeans(abs(shap_interactions))
# 
# # Remove intercept col
# nfeat <- dim(shap_interactions)[2]  # total features including BIAS
# shap_interactions <- shap_interactions[, 1:(nfeat-1), 1:(nfeat-1)]
# 
# # Format for interpretation
# interaction_strength <- apply(abs(shap_interactions), c(2,3), mean)
# # diag(interaction_strength) <- NA # set diagonal (self-interactions, i.e. main effects) to NA
# 
# # Sort interactions
# interaction_list <- reshape2::melt(interaction_strength, na.rm = TRUE)
# interaction_list <- interaction_list[order(-interaction_list$value), ]
# # head(interaction_list, 14)   # Top 10 interactions
# 
# # Plot
# interaction_df <- reshape2::melt(interaction_strength, na.rm = TRUE)
# ggplot(interaction_df, aes(Var1, Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   theme_minimal() +
#   labs(title='XGBoost SHAP Interaction Heatmap', x='', y='')

# Overall: T2m (fg10, tp)
# Winter: T2m, fg10 and pm25
# Summer: T2m, rh, fg10 and tp

