### ------------------------------------------------------------------------ ### 
### DLNM - Exposure selection (HISP)
### ------------------------------------------------------------------------ ### 

####  DEPENDENCIES  ####

# Function to run XG Boost test on subset
xgboost <- function(data, title, var_list, params=NULL) {
  
  # 1. Prepare data ----
  data <- data %>% filter(!is.na(total_patients)) %>% arrange(date) %>%
    mutate(location = as.numeric(as.factor(location)))
  
  # Choose features
  # features <- names(data)[grepl(var_list, names(data))]
  # features <- data[, var_list, drop = FALSE]
  X <- as.matrix(data[, var_list] %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))))
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
  
  # Rename variables
  feat <- importance$Feature
  is_lag <- grepl("lag", feat)
  parts  <- strsplit(feat[is_lag], "_lag")
  feat[is_lag]  <- paste0(replace_exact_matches(sapply(parts, `[`, 1), nm), " (Lag ", sapply(parts, `[`, 2), ")")
  feat[!is_lag] <- replace_exact_matches(feat[!is_lag], nm)
  importance$Feature <- unlist(feat)
  
  plot <- ggplot(importance, aes(x = Gain, y = reorder(Feature, Gain))) +
    geom_col(fill = 'gray50') +
    labs(x = "Gain", y = 'Exposure', title = title) +
    theme_classic() +
    theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(plot)
  
  return(list(best_params=best_params, importance=importance, plot=plot))
  
}

#### MAIN CODE ####

# Set parameters
control_cols <- c('time', 'season', 'weekday', 'location', 'lat', 'lon', p_cols)
cumm_cols <- grep(paste0("(", paste(w_cols, collapse="|"), ")_roll"), names(data), value=TRUE)
cumm_cols <- cumm_cols[!grepl("^tp_roll[0-9]+$", cumm_cols)] # Remove tp not log
lagged_cols <- grep(paste0("(", paste(w_cols, collapse="|"), ")_lag([1-9]|1[0-4])$"), names(data), value=TRUE)
autocorr_cols <- grep(paste0("(", paste('total_patients', collapse="|"), ")_lag([1-9]|1[0-4])$"), names(data), value=TRUE)

# Create cumulative names
first_half <- sub("_.*", "", cumm_cols)
window_num <- as.integer(sub(".*roll([0-9]+).*", "\\1", cumm_cols))
new_entries <- setNames(
  paste0(window_num, "-day ", nm[first_half]),
  cumm_cols
)
nm <- c(nm, new_entries)

# 1. Non parameter correlation test ####

sub <- data %>% dplyr::ungroup() %>% dplyr::select(total_patients, all_of(w_cols))
# sub <- data %>% ungroup() %>% dplyr::select(total_patients, all_of(c(w_cols, cumm_cols)))
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

# 2A. Cumulative: XG Boost feature importance test ####

# Select parameters
var_list <- c(w_cols, control_cols, cumm_cols, 'total_patients_lag14')

# Get best parameters & results
start_time <- Sys.time()
xg_res <- list()
xg_res[['all']] <- xgboost(data, "Overall", var_list, params=NULL)
for (l in locations){
  sub <- data %>% filter(location==l)
  xg_res[[l]] <- xgboost(sub, l, var_list, params=NULL)
}
end_time <- Sys.time()
# Print time & save (~90min)
print(difftime(end_time, start_time, units = "mins"))
save(xg_res, file = paste0(path_project, 'Results/ind_xgboost_params.RData'))

# # Test
# imp <- map_dfr( xg_res, ~ as_tibble(.x$importance), .id = "subpop")
# imp %>%
#   filter(str_detect(Feature, "Precip")) %>%
#   filter(str_detect(subpop, "all")) %>%
#   arrange(desc(Gain))

# 2B. Cumulative: Format results ####

# Get results
load('2_DLNM/Results/ind_xgboost_params.RData')

# Feature importance overall
imp <- map_dfr( xg_res[1], ~ as_tibble(.x$importance), .id = "subpop")
imp_lag_all <- ggplot(imp[1:20,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_col(fill = 'gray50') +
  labs(x = "Gain", y = 'Exposure', title = NULL) +
  theme_classic() +
  theme(plot.title = ggplot2::element_text(hjust = 0.5))

# Feature importance consistency across locations
imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
(imp_lag_loc <- imp %>%
    # Group by location
    group_by(subpop) %>%
    arrange(desc(Gain)) %>%
    mutate(
      top    = row_number() <= (n()/2),
      bottom = row_number() > (n()/2)
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
    filter(Important==TRUE) %>%
    arrange(desc(`Top Half`), desc(`Mean Gain`)))
imp_lag_loc$`Mean Gain` <- round(imp_lag_loc$`Mean Gain`, 3)
imp_lag_loc <- imp_lag_loc[1:20,-5]

# Combine figure and table
tab <- as.ggplot(tableGrob(imp_lag_loc, rows = NULL))
imp_lag <- wrap_plots(imp_lag_all, tab, nrow = 1) + plot_annotation(tag_levels = "A")

# EXPORT
ggsave(paste0(path_report, 'II. dlnm/sensitivity/SM_xgboost_imp.png'), imp_lag, width = 14, height = 7, dpi = 300, bg='white')



# # Lags: XG Boost feature importance test ####
# 
# # Select parameters
# var_list <- c(w_cols, control_cols, autocorr_cols, lagged_cols)
# 
# # # Get best parameters & results
# # start_time <- Sys.time()
# # xg_res <- list()
# # xg_res[['all']] <- xgboost(data, "Overall", var_list, params=NULL)
# # for (l in locations){
# #   sub <- data %>% filter(location==l)
# #   xg_res[[l]] <- xgboost(sub, l, var_list, params=NULL)
# # }
# # end_time <- Sys.time()
# # # Print time & save (~80min)
# # print(difftime(end_time, start_time, units = "mins"))
# # save(xg_res, file = paste0(path_project, 'Results/ind_xgboost_params_lags.RData'))
# 
# 
# # Lags: Format results ####
# 
# # # Get results
# # load('2_DLNM/Results/ind_xgboost_params_lags.RData')
# # 
# # # Feature importance overall
# # imp <- map_dfr( xg_res[1], ~ as_tibble(.x$importance), .id = "subpop")
# # imp_lag_all <- ggplot(imp[1:20,], aes(x = Gain, y = reorder(Feature, Gain))) +
# #   geom_col(fill = 'gray50') +
# #   labs(x = "Gain", y = 'Exposure', title = NULL) +
# #   theme_classic() +
# #   theme(plot.title = ggplot2::element_text(hjust = 0.5))
# # 
# # # Feature importance consistency across locations
# # imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
# # (imp_lag_loc <- imp %>%
# #     # Group by location
# #     group_by(subpop) %>%
# #     arrange(desc(Gain)) %>%
# #     mutate(
# #       top    = row_number() <= (n()/2),
# #       bottom = row_number() > (n()/2)
# #     ) %>%
# #     ungroup() %>%
# #     group_by(Feature) %>%
# #     summarise(
# #       `Top Half`    = sum(top),
# #       `Bottom Half` = sum(bottom),
# #       `Mean Gain`     = mean(Gain),
# #       .groups = "drop"
# #     ) %>%
# #     mutate(Important= `Bottom Half` < `Top Half`) %>%
# #     filter(Important==TRUE) %>%
# #     arrange(desc(`Top Half`), desc(`Mean Gain`)))
# # imp_lag_loc$`Mean Gain` <- round(imp_lag_loc$`Mean Gain`, 3)
# # imp_lag_loc <- imp_lag_loc[1:20,-5]
# # 
# # # Combine figure and table
# # tab <- as.ggplot(tableGrob(imp_lag_loc, rows = NULL))
# # imp_lag <- wrap_plots(imp_lag_all, tab, nrow = 1) + plot_annotation(tag_levels = "A")
# # 
# # # EXPORT
# # ggsave(paste0(path_report, 'II. dlnm/sensitivity/SM_imp_lag.png'), imp_lag, width = 14, height = 7, dpi = 300, bg='white')

# # Final: XG Boost feature importance test ####
# 
# # Select parameters
# var_list <- c(w_cols, control_cols, 't2m_roll7', 'tp_roll7', 'total_patients_lag14')
# 
# # Get best parameters & results
# start_time <- Sys.time()
# xg_res <- list()
# xg_res[['all']] <- xgboost(data, "Overall", var_list, params=NULL)
# for (l in locations){
#   sub <- data %>% filter(location==l)
#   xg_res[[l]] <- xgboost(sub, l, var_list, params=NULL)
# }
# end_time <- Sys.time()
# # Print time & save (~35min)
# print(difftime(end_time, start_time, units = "mins"))
# save(xg_res, file = paste0(path_project, 'Results/ind_xgboost_params_final.RData'))
# 
# 
# # Final: Format results ####
# 
# # Get results
# load('2_DLNM/Results/ind_xgboost_params_final.RData')
# 
# # Feature importance overall
# imp <- map_dfr( xg_res[1], ~ as_tibble(.x$importance), .id = "subpop")
# imp_all <- ggplot(imp[1:20,], aes(x = Gain, y = reorder(Feature, Gain))) +
#   geom_col(fill = 'gray50') +
#   labs(x = "Gain", y = 'Exposure', title = NULL) +
#   theme_classic() +
#   theme(plot.title = ggplot2::element_text(hjust = 0.5))
# 
# # Feature importance consistency across locations
# imp <- map_dfr( xg_res[-1], ~ as_tibble(.x$importance), .id = "subpop")
# (imp_rank <- imp %>%
#     group_by(subpop) %>%
#     arrange(desc(Gain)) %>%
#     mutate(
#       top    = row_number() <= (n()/2),
#       bottom = row_number() > (n()/2)
#     ) %>%
#     ungroup() %>%
#     group_by(Feature) %>%
#     summarise(
#       `Top Half`    = sum(top),
#       `Bottom Half` = sum(bottom),
#       `Mean Gain`     = mean(Gain),
#       .groups = "drop"
#     ) %>% 
#     mutate(Important= `Bottom Half` < `Top Half`) %>%
#     arrange(desc(Important), desc(`Top Half`))) 
# imp_rank$`Mean Gain` <- round(imp_rank$`Mean Gain`, 3)
# imp_rank <- imp_rank[1:20,-5]
# 
# # Combine figure and table
# imp_rank[] <- lapply(imp_rank, function(x) ifelse(is.na(x), "", as.character(x)))
# tab <- as.ggplot(tableGrob(imp_rank, rows = NULL))
# imp_lag <- wrap_plots(imp_all, tab, nrow = 1) + plot_annotation(tag_levels = "A")
# 
# # EXPORT
# ggsave(paste0(path_report, 'II. dlnm/sensitivity/SM_imp_nolag.png'), imp_lag, width = 14, height = 7, dpi = 300, bg='white')
# 
