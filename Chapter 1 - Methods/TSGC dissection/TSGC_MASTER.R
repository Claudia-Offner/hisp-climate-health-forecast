
# Set working directory
setwd('C:/Users/offne/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')

# ------------------------------------------------------------------------------ 
# -------------------------------- DEPENDENCIES --------------------------------
# ------------------------------------------------------------------------------

# LOAD PACKAGES

library(KFAS)
library(dplyr)
library(xts)

# LOAD MODEL FUNCTIONS

### Function specifying dynamic Gompertz model
get_dynamic_gompertz_model <- function(y, q = NULL, sea.type = 'trigonometric', sea.period = 7, 
                                       a1 = NULL, P1 = NULL, Q = NULL, H = NULL) {
  
  # Description:
  # Function to generate a dynamic Gompertz curve model - The model incorporates
  # a trend component with a random walk (IRW) and seasonal adjustment. The 
  # function allows various parameters to be specified for flexibility and 
  # customization.
  
  # Parameters:
  #   - y: Cumulated variable
  #   - q: Signal-to-noise ratio (optional, estimated if NULL)
  #   - sea.type: Seasonal type ('trigonometric' or 'none')
  #   - sea.period: Period of seasonality (e.g., 7 for weekly seasonality)
  #   - a1, P1: Prior means and covariance (optional, diffuse prior if NULL)
  #   - Q, H: State error variances and irregular variance (optional)
  
  # Return: State-space model, formatted as a KFS model object.
  
  # OPTIONAL: Identify state error and irregular variance if applicable
  Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
  Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
  Ht <- if (is.null(H)) { NA } else { H }
  
  # OPTIONAL: Set prior on state as ~ N(a1, P1) if a1 supplied.
  use.prior <- if (!is.null(a1)) { TRUE } else { FALSE }
  
  # 1. If we are using prior information, fit the following models depending on seasonality...
  if (use.prior) {
    if (sea.type == 'trigonometric') {
      ss_model <- SSModel(
        y ~
          SSMtrend(
            degree = 2,
            Q = list(matrix(0), matrix(Qt.slope)),
            a1 = a1[1:2],
            P1 = P1[1:2, 1:2]
          ) +
          SSMseasonal(
            period = sea.period,
            Q = Qt.seas,
            sea.type = sea.type,
            a1 = a1[3:dim(a1)[1]],
            P1 = P1[3:dim(a1)[1], 3:dim(a1)[1]]
          ),
        H = Ht
      )
      n.pars <- 0
    } else if (sea.type == 'none') {
      ss_model <- SSModel(
        y ~
          SSMtrend(
            degree = 2,
            Q = list(matrix(0), matrix(Qt.slope)),
            a1 = a1[1:2],
            P1 = P1[1:2, 1:2]
          ),
        H = Ht
      )
      n.pars <- 0
    } else {
      stop(sprintf("sea.type= '%s' not implemented", sea.type))
    }
  } 
  
  # 2. If we are NOT using prior information, fit the following models depending on seasonality...
  else {
    if (sea.type == 'trigonometric') {
      ss_model <- SSModel(
        y ~
          SSMtrend(
            degree = 2,
            Q = list(matrix(0), matrix(Qt.slope))
          ) +
          SSMseasonal(
            period = sea.period,
            Q = Qt.seas,
            sea.type = sea.type),
        H = matrix(Ht)
      )
      n.pars <- if (is.null(q)) { 3 } else { 2 }
    } else if (sea.type == 'none') {
      ss_model <- SSModel(
        y ~
          SSMtrend(
            degree = 2,
            Q = list(matrix(0), matrix(Qt.slope))
          ),
        H = matrix(Ht)
      )
      n.pars <- if (is.null(q)) { 2 } else { 1 }
    } else {
      stop(sprintf("sea.type= '%s' not implemented", sea.type))
    }
  }
  
  out <- list(model = ss_model, inits = rep(0, n.pars))
  return(out)
}

### Function estimating model noise parameters prior to kalman filter
update <- function(pars, model, q, sea.type) {
  
  # Description: 
  # Function to update the Kalman filter with new parameters - This method updates 
  # the Kalman filter to implement 'get_dynamic_gompertz_model'. The method allows 
  # for the adjustment of the following parameters:
  # - Observation noise (1 parameter)
  # - Transition equation slope
  # - Seasonal noise
  
  # Parameters:
  # - pars: A vector containing the parameters used in the model.
  # - model: The KFS model object representing the Kalman filter state.
  # - q: The signal-to-noise ratio, which is the ratio of slope variance to irregular variance.
  #      If provided, q is used to adjust the slope noise.
  # - sea.type: Specifies the type of seasonality in the model:
  #     - "trigonometric": Includes a trigonometric seasonal component.
  #     - "none": No seasonal component.
  
  # Return: KFS model object, reflecting the updated model.
  
  estH <- any(is.na(model$H))
  estQ <- any(is.na(model$Q))
  
  if ((!estH) & (!estQ)) {
    # If nothing to update then return model - this only applies to NA values though?
    return(model)
  } else {
    nparQ <- if (sea.type == 'trigonometric') { 1 } else { 0 }
    # 1. Set seasonal noise
    if (estQ) {
      Q <- as.matrix(model$Q[, , 1])
      # Update diagonal elements
      naQd <- which(is.na(diag(Q)))
      Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
      diag(Q)[naQd] <- exp(0.5 * pars[1])
      # Check for off-diagonal elements and raise error if found.
      naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
      if (length(naQnd) > 0) {
        stop("NotImplmentedError: Unexpected off-diaganol element updating")
      }
    }
    
    # 2. Set observation noise
    H <- as.matrix(model$H[, , 1])
    if (estH) {
      # Update diagonal elements
      naHd <- which(is.na(diag(H)))
      H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
      diag(H)[naHd] <- exp(0.5 * pars[(nparQ + 1)])
      model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
    }
    
    # 3. Set slope noise
    # Get index of slope, 1 before the seasonal component.
    model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
    i.slope <- 2
    # Estimate slope if no signal to noise ratio specified.
    if (is.null(q)) {
      Q.slope <- exp(0.5 * pars[(nparQ + 2)])
    } else {
      Q.slope <- crossprod(H[naHd, naHd]) * q
    }
    model$Q[i.slope, i.slope, 1] <- Q.slope
  }
  return(model)
}

# LOAD FORMATTING FUNCTIONS

### Function predicting y cases, state variables & conditional state covariance matrices
predict_all <- function(index, output, n.ahead, sea.on = FALSE, return.all = FALSE) {
  
  # Description:
  # Function returning the incidence variable y, the state variables 
  # and the conditional covariance matrix for the states.
  
  # Parameters:
  # - `n.ahead`: Number of forecasts to generate from the end of the sample period.
  # - `sea.on`: Logical; whether to include seasonal components in the state-space model (default: `TRUE`).
  # - `return.all`: Logical; whether to return all filtered estimates and forecasts (`TRUE`) 
  #   or only forecasts (`FALSE`, default).
  
  # Return Value:
  # An 'xts' object that contains:
  # - The forecasted (and filtered, if applicable) level of y (`y.hat`)
  # - The level component (`level.t.t`)
  # - The slope component (`slope.t.t`)
  # - A vector of states, including seasonals if applicable (`a.t.t`)
  # - The covariance matrix of all states, including seasonals if applicable (`P.t.t`)
  
  idx <- index
  model <- output
  
  new.model <- output$model
  new.model$y <- rbind(new.model$y, 
                       matrix(NA, ncol = ncol(new.model$y), nrow = n.ahead) %>% 
                         as.ts())
  attr(new.model, 'n') <- length(output$model$y) + n.ahead %>%
    as.integer()
  model_output <- KFS(new.model)
  
  if (sea.on == TRUE) {
    y.hat.kfas <- predict(
      output$model, interval = 'prediction',
      n.ahead = n.ahead, level = 0.68, states = 'all')
  } else {
    y.hat.kfas <- predict(
      output$model, interval = 'prediction',
      n.ahead = n.ahead, level = 0.68, states = 'level')
  }
  
  n <- attr(model$model, "n")
  dates <- seq(idx[1], by = 'day', length.out = (n + n.ahead))
  
  # Assumes time invariant Z.t
  y.t.t <- output$att %*% drop(output$model$Z)
  
  y.hat <- xts::xts(
    rbind(y.t.t, y.hat.kfas[, 1] %>% as.matrix()),
    order.by = dates)
  
  i.level <- grep("level", colnames(model_output$att))
  level.t.t <- xts::xts(model_output$att[, i.level], order.by = dates) %>%
    as.xts()
  i.slope <- grep("slope", colnames(model_output$att))
  slope.t.t <- xts::xts(model_output$att[, i.slope], order.by = dates) %>%
    as.xts()
  
  if (!return.all) {
    y.hat <- y.hat %>%
      subset(index(.) > tail(idx, 1))
    level.t.t <- level.t.t %>%
      subset(index(.) > tail(idx, 1))
    slope.t.t <- slope.t.t %>%
      subset(index(.) > tail(idx, 1))
  }
  
  out <- list(
    y.hat = y.hat,
    level.t.t = level.t.t,
    slope.t.t = slope.t.t,
    a.t.t = model_output$att,
    P.t.t = model_output$Ptt
  )
  
  return(out)
}

### Function predicting level of cumulative cases 
predict_level <- function(index, output, y.cum, n.ahead, confidence.level, return.diff = FALSE) {
  
  # Description:
  # Function to predict the level of a cumulative variable
  
  # Parameters:
  # y.cum: The cumulative variable (e.g., total cases in an epidemic)
  # n.ahead: Number of periods ahead to forecast
  # confidence.level: Confidence level for the forecast interval
  # return.diff: If TRUE, returns the difference (incidence) instead of cumulative values
  
  # Return: xts object containing the point forecasts and upper and lower bounds of the forecast interval.

  model <- output$model
  n <- attr(model, "n")
  p <- attr(model, "p")
  
  freq <- unclass(periodicity(y.cum))$label
  endtime <- end(model$y) + c(0, n.ahead)
  filtered.out <- predict_all(index, output, n.ahead, sea.on = TRUE, return.all = TRUE)
  
  # Extract parameters.
  timespan <- n + 0:n.ahead
  
  # Calculate g.t as exponent of y.t
  g.t <- exp(filtered.out$y.hat)
  
  # Forecast dates
  v_dates_end <- seq(last(index(y.cum)), last(
    index(filtered.out$y.hat)), by = freq
  )
  
  # Construct CI
  Ptt <- filtered.out$P.t.t
  i.level <- grep("level", colnames(model$T))
  ci <- qnorm((1 - confidence.level) / 2) * sqrt(Ptt[i.level, i.level, timespan]) %o% c(1, -1)
  ci <- xts(exp(ci), order.by = v_dates_end)
  
  y.hat <- xts(matrix(NA, nrow = n.ahead + 1, ncol = 3),
               order.by = v_dates_end)
  y.hat[v_dates_end[1],] <- y.cum[v_dates_end[1]]
  for (i in seq_len(length(v_dates_end[-1]))) {
    date.forecast <- v_dates_end[i + 1]
    date.lag <- date.forecast - 1
    # Update level
    y.hat[date.forecast, 1] <- as.numeric(y.hat[date.lag, 1]) *
      as.numeric(1 + g.t[date.forecast,])
    
    # Make prediction intervals
    y.hat[date.forecast, 2:3] <- as.numeric(y.hat[date.lag, 1]) *
      as.numeric(1 + g.t[date.forecast,] %*% ci[date.forecast,])
  }
  
  # Difference output if requested
  d <- if (return.diff) { diff(y.hat[, 1])[-1] } else { (y.hat[, 1])[-1] }
  
  ci_bounds <- if (return.diff) {
    (y.hat[, 2:3] - as.vector(y.hat[, 1]))[-1] + as.vector(d)
  } else { y.hat[2:dim(y.hat)[1], 2:3] }
  
  pred <- vector("list", length = p)
  pred[[p]] <- cbind(fit = d, prediction = ci_bounds)
  pred <- lapply(pred, ts, end = endtime, frequency = 1)
  
  y.hat <- xts(pred[[p]], order.by = v_dates_end[-1])
  names(y.hat)[2:3] <- list('lower', 'upper')
  
  return(as.xts(y.hat))
  
  
  
}

### Function returning table of Kalman estimated variance parameters
print_estimation_results <- function(index, output) {
  
  # Description:
  # Prints a formatted table of kalman estimated variance parameters
  # (Optional - ready to paste into LaTeX)
  
  # Parameters:
  # - None
  
  # Returns: Table
  
  H <- output$model$H[, , 1]
  Q_gamma <- output$model$Q[2, 2, 1]
  Q_seasonal <- output$model$Q[3, 3, 1]
  
  tbl <- data.frame(
    a = format(H, digits = 3),
    b = format(Q_gamma, digits = 3),
    c = format(Q_seasonal, digits = 3),
    d = format(Q_gamma / H, digits = 4))
  header.names <- c('$\\sigma_\\varepsilon^2$',
                    '$\\sigma_\\gamma^2$',
                    '$\\sigma_{seas}^2$',
                    'q')
  
  out <- tbl %>%
    kableExtra::kbl(
      caption = "Estimated parameters",
      col.names = header.names,
      # format = 'latex',
      booktabs = TRUE,
      escape = FALSE
    ) %>%
    kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria") #%>%
    # kableExtra::footnote(general = " ")
  
  return(out)
}

### Function extracting smoothed/filtered growth rate of cases and associated components (δ and γ)
get_growth_y <- function(index, output, smoothed = FALSE, return.components = FALSE) {
  
  # Description:
  # Returns the growth rate of the incidence (`y`) of the cumulated variable (`Y`).
  # Computed as: g_t = exp(δ_t) + γ_t
  
  # Parameters:
  # - `smoothed`: Logical; whether to use smoothed estimates of δ and γ (`TRUE`) 
  #   or contemporaneous filtered estimates (`FALSE`, default).
  # - `return.components`: Logical; whether to return δ, γ, and growth rate estimates (`TRUE`), 
  #   or just the growth rate (`FALSE`, default).
  
  # Return Value: An `xts` object containing smoothed/filtered growth rates and 
  # components (δ and γ), if applicable.
  
  kfs_out <- output
  idx <- index
  
  if (smoothed) {
    att <- kfs_out$alphahat
  } else {
    att <- kfs_out$att
  }
  
  filtered.slope <- xts(att[, "slope"], order.by = idx)
  filtered.level <- xts(att[, "level"], order.by = idx)
  g.t <- exp(filtered.level)
  gy.t <- g.t + filtered.slope
  names(gy.t) <- if (smoothed) { "smoothed gy.t" } else { "filtered gy.t" }
  names(g.t) <- if (smoothed) { "smoothed g.t" } else { "filtered g.t" }
  names(filtered.slope) <- if (smoothed) { "smoothed gamma.t" } else {"filtered gamma.t" }
  if (return.components) {
    return(list(gy.t, g.t, filtered.slope))
  } else {
    return(gy.t)
  }
}

### Function getting CIs of smoothed/filtered growth rate of cases
get_gy_ci <- function(index, output, smoothed = FALSE, confidence_level = 0.68) {
  
  # Description:
  # Returns the growth rate of the incidence (`y`) of the cumulated variable (`Y`).
  # Computed as: g_t = exp(δ_t) + γ_t
  
  # Parameters:
  # - `smoothed`: Logical; whether to use smoothed estimates of δ and γ (`TRUE`) 
  #   or contemporaneous filtered estimates (`FALSE`, default).
  # - `confidence_level`: Confidence level for the confidence interval. Default is 0.68 
  #   (one standard deviation for a normally distributed random variable).
  
  # Return Value:
  # - An `xts` object containing smoothed/filtered growth rates and the upper and 
  #   lower bounds for the confidence intervals.
  
  kfs_out <- output
  idx <- index
  
  if (smoothed) {
    att <- kfs_out$alphahat
  } else {
    att <- kfs_out$att
  }
  
  filtered_slope <- xts(att[, "slope"], order.by = idx)
  filtered.level <- xts(att[, "level"], order.by = idx)
  g.t <- exp(filtered.level)
  gy.t <- g.t + filtered_slope
  
  idx.slope <- grep("slope", colnames(kfs_out$att))
  ci <- qnorm((1 - confidence_level) / 2) * sqrt(kfs_out$Ptt[idx.slope, idx.slope,]) %o% c(1, -1)
  ci_bounds <- as.vector(gy.t) + ci
  
  pred <- xts(cbind(gy.t, ci_bounds), order.by = idx)
  colnames(pred) <- c("fit","lower","upper")
  
  return(pred)
}


# ------------------------------------------------------------------------------ 
# --------------------------------- PART 1 ------------------------------------- 
#     Univariate model of A&E attendances in Scotland, of a single hospital
# ------------------------------------------------------------------------------

# 0. Import & clean data ####

# Load data
df <- read.csv('201_heal/agy_adm4_NHS-Scotland_2015-2025_daily-hosptial-ae.csv')

# Select hospital
df <- df %>% 
  filter(sitecode=='A111H') %>% # relevant sites
  select(-sitecode)

# Set date type, create XTS object & select outcome
df <- df[ , -2] # Drop irrelevant data
df$date <- as.Date(df$date, format = "%Y-%m-%d")  # Adjust format as needed
Y <- xts(df[ , -1], order.by = df$date)

# Select outcome and subset desired dates
Y <- Y["2017-02-01/2020-02-01"]
Y <- Y$c_NumberOfAttendancesEpisode

# 1. Fit base model ####

# Set parameters 
estimation.date.start <- as.Date("2018-01-01")
estimation.date.end <- as.Date("2018-05-25")
n.ahead <- 14
q <- 0.005
confidence.level <- 0.68
plt.length <- 30
sea.type <- 'trigonometric'

# Set sample period & get log growth rate of cumulative series
y.cum <- Y[estimation.date.start <= index(Y) & index(Y) <= estimation.date.end]
y.lg <- tsgc::df2ldl(y.cum) 

# Fit initial state-space model
model <- get_dynamic_gompertz_model(y.lg, q)


# 2. Update  model with known state-space param estimates ####

# Partial update of known params (seasonal, observational & slope noise)
updatefn <- purrr::partial(update, ... =, q = q, sea.type = sea.type)

# Full update, using MLE for unknown model params
model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                    method = 'BFGS')

# 3. Run the Kalman Filter ####
output <- KFS(model_fit$model)

# 4. Filter results for reporting tables and figures ####

# Get truncated index from model, if using a reinitialisation in self$get_model
index <- if (!is.null(model$index)) { model$index } else { index(y.lg) }

# Check kalman noise variance parameters
print_estimation_results(index, output)

# Get the smoothed growth rate of cases, with delta & gamma parameters
get_growth_y(index, output, smoothed = TRUE, return.components = TRUE)

# Get the growth rate of cases, with CIs
get_gy_ci(index, output)

# Get predicted n.ahead cumulative cases, with CIs
predict_level(index, output, y.cum, n.ahead, confidence.level, return.diff = FALSE)

# Get full KFS model with predicted n.ahead values (no CI)
predict_all(index, output, n.ahead, sea.on = FALSE, return.all = FALSE)






# ------------------------------------------------------------------------------ 
# --------------------------------- PART 2 ------------------------------------- 
#    Multivariate model of mortality & temperature in England, nationally
# ------------------------------------------------------------------------------

# 0. Import & clean data ####

# Load data
df <- read.csv('eng_adm0_ONS_2010-2019_daily-national-mortality&temp.csv')

# Set date type, create XTS object & select outcome
df <- df[ , -2] # Drop irrelevant data
df$date <- as.Date(df$date, format = "%Y-%m-%d")  # Adjust format as needed
Y <- xts(df[ , -1], order.by = df$date)

# Select outcome and subset desired dates
Y <- Y["2017-02-01/2019-12-27"]

# 1. Fit base model ####

# Set parameters 
estimation.date.start <- as.Date("2018-01-01")
estimation.date.end <- as.Date("2018-05-25")
n.ahead <- 14
q <- 0.005
confidence.level <- 0.68
plt.length <- 30
sea.type <- 'trigonometric'
sea.period <- 7

# Set sample period & get log growth rate of cumulative series
y.cum <- Y[estimation.date.start <= index(Y) & index(Y) <= estimation.date.end]
y.lg <- tsgc::df2ldl(y.cum$c_total_deaths) 

################################################################################

y <- as.numeric(tsgc::df2ldl(y.cum$c_total_deaths))
x1 <- as.numeric(y.cum$temp)

model <- SSModel(y ~ SSMregression(~ x1, Q = diag(NA,1)) + 
                   SSMtrend(degree = 2, Q = list(matrix(0), matrix(NA))) +
                   SSMseasonal(period = sea.period, Q = NA, sea.type = sea.type), H=matrix(NA))

updatefn <- purrr::partial(update, ... =, q = q, sea.type = sea.type)
model_fit <- fitSSM(model, inits = c(0,0), updatefn = updatefn, method = 'BFGS')

output <- KFS(model_fit$model)
output$model

ts.plot(output$alphahat[,1])

index <- if (!is.null(model$index)) { model$index } else { index(y.lg) }

predict_level(index, output, y.cum$c_total_deaths, n.ahead, confidence.level, return.diff = FALSE)
