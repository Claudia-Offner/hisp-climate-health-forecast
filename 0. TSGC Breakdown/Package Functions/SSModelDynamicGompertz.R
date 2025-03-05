# Created by: Craig Thamotheram
# Created on: 27/07/2022

# ------------------------------------------------------------------------------ 
# -------------------------------- DESCRIPTION ---------------------------------
# ------------------------------------------------------------------------------
#' @title Class for Dynamic Gompertz Curve State Space Model
#'
#' @description This class represents a **dynamic Gompertz curve state-space model**.
#' It is used for modeling and forecasting time-series data where the Gompertz 
#' curve is a suitable fit.
#'
#' @section Methods:
#' - **`get_model(y, q = NULL, sea.type = "trigonometric", sea.period = 7)`**:  
#'   Retrieves the model object.
#'
#' @section Parameters:
#' - **`y`**: The **cumulated variable** (i.e., the observed data being modeled).
#' - **`q`**: The **signal-to-noise ratio** (ratio of slope to irregular variance). 
#'            Defaults to `NULL`, meaning it will be estimated automatically.
#' - **`sea.type`**: Specifies the **type of seasonality**:
#'   - `"trigonometric"` (default) → Uses a trigonometric seasonal component.
#'   - `"none"` → No seasonal component.
#' - **`sea.period`**: The **seasonality period** (e.g., `7` for weekly patterns). 
#'                     This is **not required** if `sea.type = "none"`.
#'
#' @section Return Value:
#' - Returns a **`KFS` model object**, which contains the estimated parameters  
#'   and state-space model structure.

# ------------------------------------------------------------------------------ 
# --------------------------------- MAIN CODE ---------------------------------- 
# ------------------------------------------------------------------------------
setOldClass("KFS") # Only for OOP: Register KFS as an old-style S3 class

# Load required packages
library(KFAS)
library(magrittr)
library(xts)
library(methods)

# Define the dynamic Gompertz state space model function
SSModelDynamicGompertz <- function(y, q = NULL, sea.type = "trigonometric", sea.period = 7) {
  
  # PARAMETERS:
  # - y: The cumulated variable (i.e., the observed data being modeled).
  # - q: The signal-to-noise ratio (ratio of slope to irregular variance).
  #      Defaults to NULL, meaning it will be estimated automatically.
  # - sea.type: Specifies the type of seasonality:
  #      - "trigonometric" (default) → Uses a trigonometric seasonal component.
  #      - "none" → No seasonal component.
  # - sea.period: The period of seasonality (e.g., 7 for weekly seasonality).
  #      Not required if sea.type = "none".
  
  # RETURNS: A "KFS" model object, representing the estimated state-space model.
  
  # Input validation
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.null(q) && !is.numeric(q)) stop("q must be numeric or NULL")
  if (!sea.type %in% c("trigonometric", "none")) stop("Invalid seasonal type")
  if (sea.type == "trigonometric" && (!is.numeric(sea.period) || sea.period <= 0)) {
    stop("Invalid seasonal period")
  }
  
  # Create the model object
  model <- get_dynamic_gompertz_model(y, q, sea.type, sea.period)
  
  return(model)
}

# ------------------------------------------------------------------------------ 
# ---------------------------------- EXAMPLE ----------------------------------- 
# ------------------------------------------------------------------------------

# Example usage
library(tsgc)
data(gauteng, package = "tsgc")
idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")

# Specify a model
model <- SSModelDynamicGompertz(y = gauteng[idx.est], q = 0.005)

# Output the model
print(model)
