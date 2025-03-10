# Created by: Craig Thamotheram
# Created on: 27/07/2022

# ------------------------------------------------------------------------------
# -------------------------------- DESCRIPTION ---------------------------------
# ------------------------------------------------------------------------------

#' @title Base Class for Estimating Time-Series Growth Curve Models
#'
#' @description This is the **base class** for modeling and estimating 
#' time-series growth curves. It serves as the foundation for the classes 
#' `SSModelDynamicGompertz` and `SSModelDynGompertzReinit`, which extend its 
#' functionality.
#'
#' @section Dependencies:
#' This class relies on the following packages and functions:
#' - `xts::periodicity`, `xts::last`, `xts::xts` – Functions for handling 
#'   time-series data.
#' - `zoo::index` – Retrieves time-series index values.
#' - `methods::new` – Used for creating new objects in R.
#' - `KFAS::SSModel`, `KFAS::fitSSM`, `KFAS::KFS` – Functions for state-space 
#'   modeling and Kalman filtering.

# ------------------------------------------------------------------------------
# --------------------------------- MAIN CODE ---------------------------------- 
# ------------------------------------------------------------------------------
# NB: this is a deconstruction of the base class to improve understanding

# Set OOP
setOldClass("xts")

# Create reference class
SSModelBase <- setRefClass(
  "SSModelBase",
  fields = list(Y = "xts", q = "ANY"),
  methods = list(,,,)
)

# THE FOLLOWING FUNCTIONS ARE LISTED IN METHODS LIST OF SSModelBase CLASS

# FUNCTION: `initialize` ` ####

# Initialize the class with time-series data Y and optional q parameter
initialize = function(Y, q = NULL)
{
  Y <<- Y
  q <<- q
}

# FUNCTION: `get_model` ####

# Placeholder function for getting a model, meant to be implemented in subclasses
get_model = function(y, q = NULL)
{
  return(stop('Not Implemented Error: Implement in subclass with your model'))
}

# FUNCTION: `get_dynamic_gompertz_model` ** ####

# Description:
# Function to generate a dynamic Gompertz curve model - The model incorporates
# a trend component with a random walk (IRW) and seasonal adjustment. The 
# function allows various parameters to be specified for flexibility and customization.

# Parameters:
#   - y: Cumulated variable
#   - q: Signal-to-noise ratio (optional, estimated if NULL)
#   - sea.type: Seasonal type ('trigonometric' or 'none')
#   - sea.period: Period of seasonality (e.g., 7 for weekly seasonality)
#   - a1, P1: Prior means and covariance (optional, diffuse prior if NULL)
#   - Q, H: State error variances and irregular variance (optional)

# Return: State-space model, formatted as a KFS model object.

get_dynamic_gompertz_model = function(y, q = NULL, sea.type = 'trigonometric',
                                      sea.period = 7, a1 = NULL, P1 = NULL,
                                      Q = NULL, H = NULL) {
  
  # OPTIONAL: Identify state error and irregular variance if applicable
  Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
  Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
  Ht <- if (is.null(H)) { NA } else { H }
  
  # OPTOINAL: Set prior on state as ~ N(a1, P1) if a1 supplied.
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

# FUNCTION: `update` ** ####

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

update = function(pars, model, q, sea.type) {
  
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


# FUNCTION: `estimate` ** ####

# Description:
# This function estimates the dynamic Gompertz curve model when applied to an 
# object of class SSModelDynamicGompertz or SSModelDynGompertzReinit.
# 
# Parameters:
# - sea.type: The type of seasonality.
#   - "trigonometric": Includes a trigonometric seasonal component.
#   - "none": No seasonal component.
# - sea.period: The period of seasonality (e.g., 7 for a day-of-the-week effect with daily data).
#   Not required if sea.type is "none".
# 
# Return Value: object of class FilterResults, containing estimated dynamic Gompertz curve model.

estimate = function(sea.type = 'trigonometric', sea.period = 7) {
  
  # 1. Get LDL of cumulative series Y.
  y <- tsgc::df2ldl(Y)
  
  # 2. Add update / model methods via estimated MLE of unknown params
  updatefn <- purrr::partial(
    .self$update, ... =, q = q, sea.type = sea.type
  )
  model <- .self$get_model(y, q = q, sea.type, sea.period)
  model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                      method = 'BFGS')
  
  # 3. Run smoother/filter
  model_output <- KFS(model_fit$model)
  
  # 4. Get truncated index from model, if using a reinitialisation in self$get_model
  date.index <- if (!is.null(model$index)) { model$index } else { index(y) }
  
  results <- FilterResults$new(
    index = date.index,
    output = model_output
  )
  return(results)
}

# ------------------------------------------------------------------------------
# --------------------------------- PRACTICE  ---------------------------------- 
# ------------------------------------------------------------------------------

library(KFAS)
library(tsgc)
library(zoo)

# Set parameters
estimation.date.start <- as.Date("2021-08-01")
estimation.date.end <- as.Date("2021-11-20")
X <- england$cum_cases
Y <- england$cum_deaths
n.forecasts <- 14
q <- 0.005
confidence.level <- 0.68
plt.length <- 30

# Defaults
sea.type <- 'trigonometric'
sea.period <- 7
a1 <- NULL
P1 <- NULL
Q <- NULL
H <- NULL
# Optional
Qt.slope <- NA
Qt.seas <- NA
Ht <- NA
use.prior <- FALSE

# Subset sample period 
x <- X[estimation.date.start <= index(X) & index(X) <= estimation.date.end]
y <- Y[estimation.date.start <= index(Y) & index(Y) <= estimation.date.end]

# Function SSModel creates a state space object object of class SSModel which can be used as an input object for various functions of KFAS package.
# https://rdrr.io/cran/KFAS/man/SSModel.html 

# NB: SSModel(formula, data, H, u, distribution, tol = .Machine$double.eps^0.5)
# NB: SSMseasonal(period, Q, sea.type = c("dummy", "trigonometric"), type, index, a1, P1, P1inf, n = 1, state_names = NULL, ynames, harmonics )
# NB: SSMtrend(degree = 1, Q, type, index, a1, P1, P1inf, n = 1, state_names = NULL, ynames )


model <- SSModel(y ~ x + SSMtrend(degree = 2, Q = list(matrix(0), matrix(Qt.slope))) +
                        SSMseasonal(period = sea.period, Q = Qt.seas, sea.type = sea.type), 
                    H = matrix(Ht))

# Not sure what this one is supposed to represent
pars <- if (is.null(q)) { 3 } else { 2 }

out <- list(model = ss_model, inits = rep(0, n.pars))

# NB: Get Model = fit model
# NB: Estimate = update to addresses noise in Q & H and the seasonal component 
# (smoothsness), run via MLE (fitSSM), run smoother/filter (KFS), run tsgc filterResults()





