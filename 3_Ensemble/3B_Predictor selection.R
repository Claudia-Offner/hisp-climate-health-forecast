
rolling_cv <- function(predictors, data, smoothed_obs, model_fn) {
  
  training_dates <- unique(data$date)
  preds_t7  <- c()
  preds_t14 <- c()
  obs_t7    <- c()
  obs_t14   <- c()
  
  for (t in training_dates) {
    
    train <- data %>% filter(date <= t - 2)
    if (nrow(train) < 10) next
    
    # Fit model using the provided function
    formula <- as.formula(paste("outcome ~", paste(predictors, collapse = " + ")))
    fit     <- model_fn(formula, data = train)
    
    forecast_dates <- seq(t - 1, t + 14, by = 1)
    forecast_data  <- data %>% filter(date %in% forecast_dates)
    forecast_data$pred <- predict(fit, newdata = forecast_data)
    
    obs <- smoothed_obs %>% filter(date %in% forecast_dates)
    
    pred_t7  <- forecast_data %>% filter(date == t + 7)  %>% pull(pred)
    pred_t14 <- forecast_data %>% filter(date == t + 14) %>% pull(pred)
    obs_t7_  <- obs           %>% filter(date == t + 7)  %>% pull(smoothed_value)
    obs_t14_ <- obs           %>% filter(date == t + 14) %>% pull(smoothed_value)
    
    if (length(pred_t7)  & length(obs_t7_))  { preds_t7  <- c(preds_t7,  pred_t7);  obs_t7  <- c(obs_t7,  obs_t7_)  }
    if (length(pred_t14) & length(obs_t14_)) { preds_t14 <- c(preds_t14, pred_t14); obs_t14 <- c(obs_t14, obs_t14_) }
  }
  
  rmse_t7  <- sqrt(mean((preds_t7  - obs_t7)^2,  na.rm = TRUE))
  rmse_t14 <- sqrt(mean((preds_t14 - obs_t14)^2, na.rm = TRUE))
  
  return(list(rmse_t7 = rmse_t7, rmse_t14 = rmse_t14))
}


forward_selection <- function(covariates, training_data, smoothed_obs, model_fn) {
  
  # --- Step 1: Univariate screening ---
  
  univariate_rmse <- tibble(covariate = covariates, rmse_t7 = NA, rmse_t14 = NA)
  
  for (cov in covariates) {
    rmse_results <- rolling_cv(cov, training_data, smoothed_obs, model_fn)
    univariate_rmse[univariate_rmse$covariate == cov, "rmse_t7"]  <- rmse_results$rmse_t7
    univariate_rmse[univariate_rmse$covariate == cov, "rmse_t14"] <- rmse_results$rmse_t14
  }
  
  univariate_rmse <- univariate_rmse %>%
    mutate(rmse_combined = (rmse_t7 + rmse_t14) / 2) %>%
    arrange(rmse_combined)
  
  top3_starts <- univariate_rmse$covariate[1:3]
  
  
  # --- Step 2: Forward selection from each starting point ---
  
  run_forward_selection <- function(starting_covariate) {
    
    selected     <- starting_covariate
    remaining    <- setdiff(covariates, starting_covariate)
    current_rmse <- univariate_rmse %>%
      filter(covariate == starting_covariate) %>%
      pull(rmse_combined)
    
    repeat {
      candidate_rmse <- tibble(covariate = remaining, rmse_combined = NA)
      
      for (cov in remaining) {
        rmse_results <- rolling_cv(c(selected, cov), training_data, smoothed_obs, model_fn)
        candidate_rmse[candidate_rmse$covariate == cov, "rmse_combined"] <-
          (rmse_results$rmse_t7 + rmse_results$rmse_t14) / 2
      }
      
      best_candidate <- candidate_rmse %>% arrange(rmse_combined) %>% slice(1)
      
      if (current_rmse - best_candidate$rmse_combined <= 1) break
      
      selected     <- c(selected, best_candidate$covariate)
      remaining    <- setdiff(remaining, best_candidate$covariate)
      current_rmse <- best_candidate$rmse_combined
      
      if (length(remaining) == 0) break
    }
    
    return(list(predictors = selected, rmse = current_rmse))
  }
  
  
  # --- Step 3: Run from each of the 3 starting points ---
  
  results <- lapply(top3_starts, run_forward_selection)
  
  
  # --- Step 4: Return the best of the 3 models ---
  
  best_model <- results[[which.min(sapply(results, function(x) x$rmse))]]
  return(best_model)
}

# Load data
load(paste0(path_project, file_import, '.RData'))
data_rt <- df; rm(df)

# Set parameters
training_data <- data_rt %>% filter(between(date, as.Date("2017-01-01"), as.Date("2017-06-01")))
covariates <- c('t2m_smooth', 'tp_smooth')
outcome <- 'total_patients_smooth'
smoothed_obs <- 'total_patients_smooth'

# Linear model
best <- forward_selection(covariates, training_data, smoothed_obs, model_fn = lm)

# Negative binomial (e.g. for count outcomes)
best <- forward_selection(covariates, training_data, smoothed_obs,
                          model_fn = function(formula, data) MASS::glm.nb(formula, data = data))

# Poisson
best <- forward_selection(covariates, training_data, smoothed_obs,
                          model_fn = function(formula, data) glm(formula, data = data, family = poisson()))