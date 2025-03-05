# Created by: Craig Thamotheram
# Created on: 27/07/2022


# ------------------------------------------------------------------------------ 
# -------------------------------- DESCRIPTION ---------------------------------
# ------------------------------------------------------------------------------

#' @title Class for Re-initialised Dynamic Gompertz Curve Model
#'
#' @description This class allows for the implementation of the **re-initialisation** procedure 
#' described in the vignette. In this model, let:
#' - `t=r` be the re-initialisation date.
#' - `r_0` is the date where the cumulative series is set to 0.
#'
#' The growth rate of cumulative cases is defined as:
#' \eqn{g_t \equiv \frac{y_t}{Y_{t-1}}}
#'
#' Re-initialisation occurs at the date `r`, and the formula for cumulative cases after re-initialisation is:
#' \deqn{Y_{t}^{r} = Y_{t-1}^{r} + y_t \quad t = r, \dots, T}
#' The cumulative cases are reset to 0 at `r_0 = r-1`, and we re-initialise the model by specifying the 
#' appropriate prior distribution for the initial states.
#'
#' @section Methods:
#' - **`new(Y, q = NULL, reinit.date = NULL, original.results = NULL, use.presample.info = TRUE)`**:
#'   Creates an instance of the `SSModelDynGompertzReinit` class.
#'   - **Parameters**:
#'     - `Y`: The cumulative variable.
#'     - `q`: The signal-to-noise ratio (ratio of slope to irregular variance). Defaults to `NULL`, 
#'       in which case it will be estimated.
#'     - `reinit.date`: The re-initialisation date `r`. Should be an object of class `"Date"`. 
#'       This parameter must be specified.
#'     - `original.results`: Instead of re-estimating the model up to the `reinit.date`, a 
#'       `FilterResults` object can be provided here to take the reinitialisation parameters from this object. 
#'       Default is `NULL`.
#'     - `use.presample.info`: Logical value indicating whether to use information from before the reinitialisation 
#'       date. Default is `TRUE`. If `FALSE`, the model is estimated from scratch from the reinitialisation date.
#' - **`get_model(y, q = NULL, sea.type = NULL, sea.period)`**:
#'   Retrieves the model object, which is a dynamic Gompertz curve model re-initialised at the 
#'   `self$reinit.date`.
#'   - **Parameters**:
#'     - `y`: The cumulative variable.
#'     - `q`: The signal-to-noise ratio. Defaults to `NULL`, and will be estimated if not provided.
#'     - `sea.type`: Seasonal type. Options are:
#'       - `"trigonometric"`: Includes a trigonometric seasonal component.
#'       - `"none"`: No seasonal component.
#'     - `sea.period`: The period of seasonality. For example, for a day-of-the-week effect with daily data, this would be 7. 
#'       Not required if `sea.type = 'none'`.
#'   - **Return Value**: Returns a `KFS` model object.

# ------------------------------------------------------------------------------ 
# --------------------------------- MAIN CODE ---------------------------------- 
# ------------------------------------------------------------------------------
# NB: this is a deconstruction of the base class to improve understanding

# Set OOP
setOldClass("KFS")

# Create reference class
SSModelDynGompertzReinit <- setRefClass(
  "SSModelDynGompertzReinit",
  contains="SSModelBase",
  fields = list(
    reinit.date = "ANY",
    original.results = "ANY",
    use.presample.info = "ANY"
  ),
  methods = list()
)

# THE FOLLOWING FUNCTIONS ARE LISTED IN METHODS LIST OF SSModelDynGompertzReinit CLASS


# FUNCTION: `initialize` ####
#' @title Create an Instance of the SSModelDynGompertzReinit Class
#' 
#' @description This function creates an instance of the `SSModelDynGompertzReinit` class, which
#' allows for the re-initialisation of the dynamic Gompertz curve model.
#'
#' @section Parameters:
#' - **`Y`**: The cumulative variable (the dependent variable).
#' - **`q`**: The signal-to-noise ratio, which is the ratio of the slope to irregular variance. 
#'   Defaults to `NULL`, in which case the ratio will be estimated instead of being imposed.
#' - **`reinit.date`**: The reinitialisation date `r`. This should be specified as an object of class 
#'   `"Date"`. It must be provided if `original.results = NULL` and `use.presample.info = TRUE`.
#' - **`original.results`**: A `KFS` results object that can be used in place of a reinitialisation date. 
#'   If provided, the reinitialisation parameters will be taken from this object. 
#'   This parameter must be specified if `reinit.date = NULL` and `use.presample.info = TRUE`.
#' - **`use.presample.info`**: A logical value indicating whether information from before the reinitialisation 
#'   date should be used in the reinitialisation procedure. The default is `TRUE`. If set to `FALSE`, the model is 
#'   estimated from scratch from the reinitialisation date and no information prior to this date is used.
#'
#' @section Usage:
#' To create an instance of the `SSModelDynGompertzReinit` class:
#' \code{SSModelDynGompertzReinit$new(y, q = 0.005, reinit.date = as.Date("2021-05-12", format = date.format))}


initialize = function(Y, q = NULL, reinit.date=NULL, original.results=NULL, use.presample.info=TRUE)
{
  reinit.date <<- reinit.date
  original.results <<- original.results
  use.presample.info <<- use.presample.info
  callSuper(Y, q)
}

# FUNCTION: `get_model` ####

#' @title Retrieve the Reinitialised Dynamic Gompertz Curve Model
#'
#' @description This method retrieves the model object, which is a dynamic
#' Gompertz curve model reinitialised at the `self$reinit.date`.
#'
#' @section Parameters:
#' - **`y`**: The cumulative variable (the dependent variable).
#' - **`q`**: The signal-to-noise ratio, which is the ratio of the slope to 
#'   irregular variance. Defaults to `NULL`, in which case the ratio is estimated 
#'   instead of being imposed.
#' - **`sea.type`**: The type of seasonality. Options include:
#'   - `'trigonometric'`: This will yield a model with a trigonometric seasonal component.
#'   - `'none'`: This will yield a model with no seasonal component.
#' - **`sea.period`**: The period of seasonality. For example, for a day-of-the-week 
#'   effect with daily data, this would be `7`. This parameter is not required 
#'   if `sea.type = 'none'`.
#'
#' @section Return Value:
#' The function returns a `KFS` model object, which is the fitted model based 
#' on the dynamic Gompertz curve reinitialised at `self$reinit.date`.

get_model = function(y, q=NULL, sea.type = NULL, sea.period)
{
  # 4.1. Index for reinitialisation, t_0
  stopifnot(length(Y[reinit.date]) == 1)
  Y.t.r_0 <- as.numeric(Y[reinit.date - 1])
  
  # 4.2 Reinitialisation:
  #   ln g_t^r = ln g_t + ln (Y_{t-1}/Y_{t-1}^r), where Y_t^r=Y_t-Y_{r_0}.
  idx.dates <- (index(y) > reinit.date)
  lag.Y <- stats::lag(Y)[idx.dates]
  y.reinit <- y[index(y) > reinit.date] + log(lag.Y / (lag.Y - Y.t.r_0))
  
  # 4.3 Run Kalman filter/smoother on new series with non-diffuse prior
  if (use.presample.info) {
    # Either estimate full model here or take results from previous model.
    if (is.null(original.results)) {
      # NB. Restrict sample to t<=r - date of reinitialisation.
      idx.est <- zoo::index(Y) <= reinit.date
      model <- SSModelDynamicGompertz$new(Y = Y[idx.est], q = q)
      res.original <- model$estimate(sea.type = sea.type,
                                     sea.period = sea.period)
      model_output <- res.original$output
    } else {
      model_output <- original.results$output
      model_seasonal <- attr(
        model_output$model$terms, "specials")$SSMseasonal
      sea.type <- if (is.null(model_seasonal)) {'none'} else {
        'trigonometric'}
      sea.period <- if (!is.null(model_seasonal)) {
        ncol(model_output$att) - 1}
    }
    
    # 4.3 Reset slope to 0 and add constant to initial value for level.
    # where reinit.date is t=r
    idx <- which(reinit.date == index(y))
    stopifnot(length(idx) == 1)
    att <- model_output$att[idx,]
    Ptt <- model_output$Ptt[, , idx]
    Tt <- drop(model_output$model$T)
    Rt <- drop(model_output$model$R)
    Qt <- drop(model_output$model$Q)
    Ht <- drop(model_output$model$H)
    
    # a. Take a_{r|r} and P_{r|r} through prediction step to get a_{r+1}
    # and P_{r+1}
    a1 <- Tt %*% att
    P1 <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
    
    # b. Set slope to 0 and add correction (\ln(Y_r/y_r) to level.
    a1["slope",] <- 0
    a1["level",] <- a1["level",] + log(Y[idx] / (Y[idx] - Y.t.r_0))
  } else {
    # Don't use presample info
    a1 <- NULL; P1 <- NULL; Qt <- NULL; Ht <- NULL
  }
  out <- .self$get_dynamic_gompertz_model(
    y = y.reinit, q = q, sea.type = sea.type, sea.period = sea.period,
    a1 = a1, P1 = P1, Q = Qt, H = Ht
  )
  out[['index']] <- index(y.reinit)
  return(out)
}

# EXAMPLES ####
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2021-05-20")
#'
#' # Specify a model
#' model.reinit <- SSModelDynGompertzReinit$new(Y = gauteng[idx.est], q = 0.005,
#'   reinit.date = as.Date("2021-04-29"))
#' # Estimate the specified model
#' res.reinit <- model.reinit$estimate()
#'
#' ## Alternatively, use prior results instead of specifying a reinitialisation date.
#'
#' # Specify the initial model
#' idx.orig <- zoo::index(gauteng) <= as.Date("2021-04-29")
#' model.orig <- SSModelDynamicGompertz$new(Y = gauteng[idx.orig], q = 0.005)
#' res.orig <- model.orig$estimate()
#' # Estimate the initial model
#' model.reinit2 <- SSModelDynGompertzReinit$new(Y = gauteng[idx.est],
#' q = 0.005, reinit.date = as.Date("2021-04-29"), original.results = res.orig)
#' res.reinit2 <- model.reinit2$estimate()
#'