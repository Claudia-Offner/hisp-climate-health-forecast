### ------------------------------------------------------------------------ ### 
### TSGC: Time series growth curves
### ------------------------------------------------------------------------ ### 

library(quantmod); library(KFAS); library(xts); library(abind); library(tsgc)

# Function to update model noise
update <- function(pars, model, q, sea.type) {
  # Update Kalman filter noise parameters via MLE
  #
  # Sets seasonal, observation, and slope noise in a KFAS SSModel. Called by
  # fitSSM() during optimisation; NAs in H/Q indicate parameters to estimate.
  #
  # @param pars  Numeric vector of log-scale noise parameters (from BFGS).
  # @param model SSModel object with NA placeholders for free parameters.
  # @param q     Signal-to-noise ratio (slope/irregular variance). If non-NULL,
  #              slope noise is derived from observation noise rather than estimated.
  # @param sea.type "trigonometric" (adds 1 seasonal param) or "none".
  # @return Updated SSModel with H and Q matrices filled.
  
  if (!any(is.na(model$H)) && !any(is.na(model$Q))) return(model)
  
  nparQ  <- as.integer(sea.type == "trigonometric")
  Q      <- as.matrix(model$Q[,,1])
  naQd   <- which(is.na(diag(Q)))
  
  # Seasonal noise
  Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
  diag(Q)[naQd] <- exp(0.5 * pars[1])
  if (any(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd])))
    stop("NotImplementedError: Unexpected off-diagonal element when updating")
  model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
  
  # Observation noise
  H    <- as.matrix(model$H[,,1])
  naHd <- which(is.na(diag(H)))
  H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
  diag(H)[naHd] <- exp(0.5 * pars[nparQ + 1])
  model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
  
  # Slope noise
  model$Q[2, 2, 1] <- if (!is.null(q)) crossprod(H[naHd, naHd]) * q else
    exp(0.5 * pars[nparQ + 2])
  
  model
}

# Function to forecast cumulative-level values
predict_tsgc <- function(output, y.level.est, n.ahead, new.data, confidence.level, sea.on) {
  
  # Forecast cumulative hospitalisations from a fitted TSGC model
  #
  # Extends the fitted SSModel n.ahead steps, re-fits noise parameters on the
  # forecast horizon, then recursively builds the cumulative level by applying
  # the predicted daily growth rate to the last observed value.
  #
  # Supports both univariate (ncol(att) == 8) and multivariate models; the
  # latter appends covariate observations via the Z array.
  #
  # @param output          KFS output from the fitted model.
  # @param y.level.est     xts of observed cumulative values (training window).
  # @param n.ahead         Integer forecast horizon (days).
  # @param new.data        xts of covariate values over the forecast window.
  # @param confidence.level Coverage probability for prediction intervals (e.g. 0.68).
  # @param sea.on          If TRUE, include seasonal states when predicting.
  # @return xts with columns fit / lower / upper (cumulative level scale).

  multivariate  <- ncol(output$att) > 8
  extend.model  <- output$model
  extend.model$y <- rbind(extend.model$y,
                          matrix(NA, nrow = n.ahead, ncol = ncol(extend.model$y)) |> as.ts())
  attr(extend.model, "n") <- nrow(extend.model$y)
  
  if (multivariate) {
    new.Z <- array(NA, dim = c(1, ncol(extend.model$Z), n.ahead))
    new.Z[, 2:ncol(extend.model$Z), ] <- extend.model$Z[, 2:ncol(extend.model$Z), 1]
    extend.model$Z <- abind(extend.model$Z, new.Z, along = 3)
  }
  
  model_output <- KFS(extend.model)
  
  new <- if (multivariate) {
    SSModel(rep(NA, n.ahead) ~
              SSMregression(~ new.data, Q = diag(0, length(covar_all))) +
              SSMtrend(degree = 2, Q = list(matrix(0), matrix(NA))) +
              SSMseasonal(period = sea.period, Q = NA, sea.type = sea.type), H = matrix(NA))
  } else {
    SSModel(rep(NA, n.ahead) ~
              SSMtrend(degree = 2, Q = list(matrix(0), matrix(NA))) +
              SSMseasonal(period = sea.period, Q = NA, sea.type = sea.type), H = matrix(NA))
  }
  
  # updatefn <- purrr::partial(update, ..., q = q, sea.type = sea.type)
  updatefn <- function(pars, model) update(pars, model, q = q, sea.type = sea.type)
  new_fit  <- fitSSM(new, inits = c(0, 0), updatefn = updatefn, method = "BFGS")
  
  # Growth-rate trajectory
  dates  <- seq(index(y.level.est)[1], by = "day", length.out = nrow(y.level.est) + n.ahead)
  y.t.t  <- output$att %*% output$model$Z[1,,1]
  y.hat.kfas <- predict(output$model, interval = "confidence",
                        newdata = new_fit$model, level = confidence.level,
                        states = if (sea.on) "all" else "level")
  g.t <- exp(xts(rbind(y.t.t, as.matrix(y.hat.kfas[, 1])), order.by = dates))
  
  # Confidence intervals from Ptt
  v_dates_end <- seq(last(index(y.level.est)), last(index(g.t)), by = "day")
  i.level     <- grep("level", colnames(output$model$T))
  timespan    <- nrow(y.level.est) + 0:n.ahead
  ci <- xts(qnorm((1 - confidence.level) / 2) *
              sqrt(model_output$Ptt[i.level, i.level, timespan]) %o% c(1, -1),
            order.by = v_dates_end)
  
  # Recursive cumulative forecast
  y.hat <- xts(matrix(NA, nrow = n.ahead + 1, ncol = 3), order.by = v_dates_end)
  y.hat[1, ] <- as.numeric(dat_train_reg$outcome[nrow(dat_train_reg)])
  
  for (i in seq_len(n.ahead)) {
    d  <- v_dates_end[i + 1]
    d0 <- v_dates_end[i]
    y.hat[d, 1]   <- as.numeric(y.hat[d0, 1]) * as.numeric(1 + g.t[d, ])
    y.hat[d, 2:3] <- as.numeric(y.hat[d0, 1]) * as.numeric(1 + g.t[d, ] %*% ci[d, ])
  }
  
  y.hat <- y.hat[-1, ]
  names(y.hat) <- c("fit", "lower", "upper")
  y.hat
}

# Set Parameters 
covar_all   <- c("t2m", "tp")
q           <- 0.005
n.ahead     <- max_prediction_horizon
sea.type    <- "trigonometric"
sea.period  <- 7
sea.on      <- FALSE
alpha_levels <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)

# Create cumulative hospitalisation variable
data_rt <- data_rt |>
  dplyr::group_by(date_analysis, location) |>
  dplyr::mutate(outcome_c = cumsum(replace_na(outcome, 1)))

#### MAIN CODE #### 

res <- tibble()

# For each date in analysis...
for (iDate in seq_along(dates_analysis)) {
  
  dat_train_ext <- data_rt |> filter(date_analysis == dates_analysis[iDate])
  dat_train <- data_rt |> filter(date_analysis == dates_analysis[iDate],
                                 date <= dates_analysis[iDate]-2)
  
  # Identify best lag per covariate
  best_lag <- lapply(sort(covar_all), function(col) {
    lapply(1:12, function(lag) {
      sub <- left_join(
        dat_train[, c("date", "location", "outcome_c")],
        transform(dat_train[, c("date", "location", col)], date = date + lag),
        by = c("date", "location")
      ) |> drop_na()
      data.frame(pearson = cor(sub$outcome_c, sub[[col]]), lag = lag, var = col)
    }) |> bind_rows()
  }) |>
    bind_rows() |>
    group_by(var) |>
    slice(which.max(abs(pearson)))
  
  list_lag <- best_lag$lag
  list_col <- best_lag$var
  
  # For each location...
  for (myloc in locations) {

    print(paste(model, myloc, dates_analysis[iDate]))
    dat_train_ext_reg <- filter(dat_train_ext, location == myloc)
    dat_train_reg <- filter(dat_train, location == myloc)
    
    # Apply covariate lags
    sub <- dat_train_reg[, c("date", "location", "outcome_c")]
    for (i in seq_along(list_lag)) {
      lagged <- transform(dat_train_reg[, c("date", "location", list_col[i])], date = date + list_lag[i])
      sub <- left_join(sub, lagged[, c("date", "location", list_col[i])], by = c("date", "location"))
    }
    sub <- drop_na(sub, outcome_c, all_of(list_col))
    
    # Fit state-space model
    Y    <- xts(subset(sub, select = c("outcome_c", covar_all)), order.by = sub$date)
    y.lg <- tsgc::df2ldl(Y[, "outcome_c"])
    y    <- as.numeric(y.lg)
    lapply(covar_all, \(v) assign(v, as.numeric(Y[, v]), envir = .GlobalEnv))
    
    mod <- SSModel(y ~
                     SSMregression(as.formula(paste("~", paste(covar_all, collapse = "+"))),
                                   Q = diag(0, length(covar_all))) +
                     SSMtrend(degree = 2, Q = list(matrix(0), matrix(NA))) +
                     SSMseasonal(period = sea.period, Q = NA, sea.type = sea.type),
                   H = matrix(NA))
    
    updatefn <- function(pars, model) update(pars, model, q = q, sea.type = sea.type)
    output <- KFS(fitSSM(mod, inits = c(0, 0), updatefn = updatefn, method = "BFGS")$model)
    
    # Get new data for prediction (last-value carry-forward for NAs)
    newdata <- data.frame(
      date = seq(dates_analysis[iDate] - remove_last_n + 1, dates_analysis[iDate] + 14, by = "day"),
      outcome_c = NA
    )
    for (i in seq_along(list_lag)) {
      lagged <- transform(dat_train_ext_reg[, c("date", list_col[i])], date = date + list_lag[i])
      newdata <- left_join(newdata, lagged, by = "date")
      col <- ncol(newdata)
      newdata[is.na(newdata[, col]), col] <- last(newdata[!is.na(newdata[, col]), col])
    }
    new.data    <- xts(subset(newdata, select = covar_all), order.by = newdata$date)
    y.level.est <- Y[as.Date(if (!is.null(mod$index)) mod$index else index(y.lg))]
    
    # Predict across alpha levels
    pred <- array(NA, dim = c(n.ahead, length(alpha_levels) * 2 + 1))
    for (i in seq_along(alpha_levels)) {
      yh <- predict_tsgc(output, y.level.est, n.ahead, new.data,
                         confidence.level = 1 - alpha_levels[i] / 100, sea.on)
      if (i == 1) pred[, 1] <- round(as.numeric(yh[, "fit"]))
      pred[, 2*i]   <- round(as.numeric(yh[, "lower"]))
      pred[, 2*i+1] <- round(as.numeric(yh[, "upper"]))
    }
    colnames(pred) <- c("point", paste0(c("lower_", "upper_"), rep(alpha_levels, each = 2)))
    
    res <- rbind(res, data.frame(
      var                = "outcome",
      location             = myloc,
      date               = dates_analysis[iDate],
      prediction_horizon = seq(1, n.ahead) - remove_last_n,
      pred
    ))
  }
}

res$model <- model
write.csv(res, file = paste0(path_results, model, "_period", period, ".csv"), row.names = FALSE)



