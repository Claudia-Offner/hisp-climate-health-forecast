### ------------------------------------------------------------------------ ### 
###  DNLM - TEST XGBoost algorithm for variable selection                    ###
### ------------------------------------------------------------------------ ### 

set.seed(123)

# LOAD THE PACKAGE
# install.packages(c('xgboost', 'SHAPforxgboost', 'data.table', 'dplyr', 'caret'))
library(xgboost); library(dplyr); library(data.table); library(caret)

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# MAIN CODE ####

# Load data
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

all_data <- df

# all_data <- all_data %>%
#   # filter(month %in% c('Jan', 'Feb', 'Mar', 'Oct', 'Nov', 'Dec'))
#   filter(month %in% c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'))
#   # filter(season %in% c('Monsoon (Jun-Sep)'))

# 1. Prepare data ----
all_data <- all_data %>% filter(!is.na(total_patients)) %>% arrange(date) %>%
  mutate(location = as.numeric(as.factor(location)))


# Choose features
features <- names(all_data)[grepl('season|weekday|time|patients_lag|elev|popu|t2m|mn2t|mx2t|wb|dt|tp|rh|pm25|fg10', names(all_data))]
X <- as.matrix(all_data[, features] %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))))
y <- all_data$total_patients

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
print(importance)
xgb.plot.importance(importance)

# # 6. Interaction strength ----
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

