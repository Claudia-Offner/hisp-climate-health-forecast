### ------------------------------------------------------------------------ ### 
### GAM2: Multivariate Generalized additive model
### ------------------------------------------------------------------------ ### 

# Set Parameters 
covar_all <- c("t2m", "tp") 
alpha     <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)
res       <- tibble()

# For each date in analysis ...
for (iDate in seq_along(dates_analysis)) {
  
  dat_train_ext <- data_rt |> filter(date_analysis == dates_analysis[iDate])
  dat_train <- data_rt |> filter(date_analysis == dates_analysis[iDate],
                                 date <= dates_analysis[iDate]-2)
  
  # Identify best lag per covariate (max absolute Pearson correlation, lags 1–12)
  best_lag <- lapply(sort(covar_all), function(col) {
    lapply(1:12, function(lag) {
      sub <- left_join(
        dat_train[, c("date", "location", "outcome_r")],
        transform(dat_train[, c("date", "location", col)], date = date + lag),
        by = c("date", "location")
      ) |> drop_na()
      data.frame(pearson = cor(sub$outcome_r, sub[[col]]), lag = lag, var = col)
    }) |> bind_rows()
  }) |>
    bind_rows() |>
    group_by(var) |>
    slice(which.max(abs(pearson)))
  
  list_lag <- best_lag$lag
  list_col <- best_lag$var
  
  # Apply best lags and fit GAM across all locations
  sub <- dat_train[, c("date", "location", "outcome_r")]
  for (i in seq_along(list_lag)) {
    lagged <- transform(dat_train[, c("date", "location", list_col[i])],
                        date = date + list_lag[i])
    sub <- left_join(sub, lagged[, c("date", "location", list_col[i])],
                     by = c("date", "location"))
  }
  sub <- sub |> tidyr::drop_na(outcome_r, all_of(list_col)) |> dplyr::select(-date, -location)
  # sub$time<-1:nrow(sub)
  
  mod <- gam(as.formula(paste0("outcome_r ~ ", paste0("s(", list_col, ")", collapse = "+"))), data = sub) #s(time,bs='ps',k=round(nrow(sub)/10)) + 
  summary(mod)

  # For each location ...
  for (myloc in locations) {
    
    print(paste(model, myloc, dates_analysis[iDate]))
    dat_train_ext_reg <- filter(dat_train_ext, location == myloc)
    dat_train_reg <- filter(dat_train, location == myloc)
    
    # Build forecast-window covariates with last-value carry-forward for NAs
    newdata <- data.frame(
      date    = seq(dates_analysis[iDate] - remove_last_n + 1, dates_analysis[iDate] + 14, by = "day"),
      outcome_r = NA
    )
    for (i in seq_along(list_lag)) {
      lagged  <- transform(dat_train_ext_reg[, c("date", list_col[i])], date = date + list_lag[i])
      newdata <- left_join(newdata, lagged, by = "date")
      col     <- ncol(newdata)
      newdata[is.na(newdata[, col]), col] <- last(newdata[!is.na(newdata[, col]), col])
    }
    
    # GAM predictions: point estimate + symmetric CIs for each alpha level
    proj <- predict(mod, newdata, se.fit = TRUE)
    z    <- qnorm(1 - (alpha / 100) / 2)
    ci   <- mapply(\(z_i) cbind(proj$fit - z_i * proj$se.fit,
                                proj$fit + z_i * proj$se.fit),
                   z, SIMPLIFY = FALSE)
    
    # Convert log-growth-rate predictions to hospitalisation counts
    pred        <- matrix(NA, nrow = max_prediction_horizon, ncol = length(alpha) * 2 + 1)
    cum_gr      <- cumsum(proj$fit)
    pred[, 1]   <- dat_train_reg$outcome[nrow(dat_train_reg)]    * exp(cum_gr)
    for (i in seq_along(alpha)) {
      pred[, 2*i]   <- dat_train_reg$outcome[nrow(dat_train_reg)] * exp(cumsum(ci[[i]][, 1]))
      pred[, 2*i+1] <- dat_train_reg$outcome[nrow(dat_train_reg)] * exp(cumsum(ci[[i]][, 2]))
    }
    pred <- apply(pred, 2, \(x) pmax(0, round(x)))
    colnames(pred) <- c("point", paste0(c("lower_", "upper_"), rep(alpha, each = 2)))
    
    res <- rbind(res, data.frame(
      var                = "outcome",
      location             = myloc,
      date               = dates_analysis[iDate],
      prediction_horizon = seq(1, max_prediction_horizon) - remove_last_n,
      pred
    ))
  }
}

res$model <- model
write.csv(res, file = paste0(path_results, model, "_period", period, ".csv"), row.names = FALSE)



