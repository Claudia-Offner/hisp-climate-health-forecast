### ------------------------------------------------------------------------ ### 
### TSGC: Time series growth curves
### ------------------------------------------------------------------------ ### 

library(quantmod); library(KFAS); library(xts); library(abind); library(tsgc);

### Function estimating model noise parameters prior to Kalman filter
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
  #     - 'trigonometric': Includes a trigonometric seasonal component.
  #     - 'none': No seasonal component.
  
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
        stop('NotImplmentedError: Unexpected off-diaganol element updating')
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

### Function predicting y cases, state variables & conditional state covariance matrices
predict_all_MULTI <- function(index, output, n.ahead, new.data=NULL, sea.on=FALSE, return.all=FALSE) {
  
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
  
  # Extend matrices for n.ahead KFS & SSM
  extend.model <- output$model

  if(ncol(output$att)==8){

    # UNIVARIATE - Extend y & n
    extend.model$y <- rbind(extend.model$y, matrix(NA, ncol=ncol(extend.model$y), nrow=n.ahead) %>% as.ts())
    attr(extend.model, 'n') <- length(output$model$y) + n.ahead %>% as.integer()
    # Fit KFS & SSM
    model_output <- KFS(extend.model)
    new <- SSModel(rep(NA, n.ahead) ~
                     SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                     SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))

  } else if (ncol(output$att) > 8) {

    # MULTIVARIATE - Extend y, z & n
    extend.model$y <- rbind(extend.model$y, matrix(NA, ncol=ncol(extend.model$y), nrow=n.ahead) %>% as.ts())
    new.Z <- array(NA, dim=c(1, ncol(extend.model$Z), n.ahead))
    new.Z[, 2:ncol(extend.model$Z), ] <- extend.model$Z[, 2:ncol(extend.model$Z), 1]
    extend.model$Z <- abind(extend.model$Z, new.Z, along=3)
    attr(extend.model, 'n') <- length(output$model$y) + n.ahead %>% as.integer()
    # Fit KFS & SSM
    model_output <- KFS(extend.model)
    new <- SSModel(rep(NA, n.ahead) ~ SSMregression(~ new.data, Q=diag(0, length(covar_all))) + #Q=diag(0,length(c(exposures, controls)))
                     SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                     SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))

  }
  # Get n.ahead forecast with new fit
  updatefn <- purrr::partial(update, ... =, q=q, sea.type=sea.type)
  new_fit <- fitSSM(new, inits=c(0, 0), updatefn=updatefn, method='BFGS')
  states_arg <- if (sea.on) 'all' else 'level'
  y.hat.kfas <- predict(output$model, interval='confidence',newdata=new_fit$model, level=0.68, states=states_arg)
  
  # Get full dates (n + n.ahead)
  n <- attr(output$model, 'n')
  dates <- seq(index[1], by='day', length.out=(n + n.ahead))
  
  # Combine fitted values (y.t.t) with forecast (y.hat.kfas)
  y.t.t <- output$att %*% output$model$Z[1, , 1]
  y.hat <- xts::xts(rbind(y.t.t, y.hat.kfas[, 1] %>% as.matrix()), order.by=dates)
  
  # Extract level and slope components
  i.level <- grep('level', colnames(model_output$att))
  level.t.t <- xts::xts(model_output$att[, i.level], order.by=dates) %>% as.xts()
  i.slope <- grep('slope', colnames(model_output$att))
  slope.t.t <- xts::xts(model_output$att[, i.slope], order.by=dates) %>% as.xts()
  
  # Filter only to the forecast period (by removing sample indexed dates)
  if (!return.all) {
    y.hat <- y.hat %>% subset(index(.) > tail(index, 1))
    level.t.t <- level.t.t %>% subset(index(.) > tail(index, 1))
    slope.t.t <- slope.t.t %>% subset(index(.) > tail(index, 1))
  }
  
  # Final output
  out <- list(
    y.hat=y.hat,
    level.t.t=level.t.t,
    slope.t.t=slope.t.t,
    a.t.t=model_output$att,
    P.t.t=model_output$Ptt
  )
  
  return(out)
}

### Function predicting level of cumulative cases 
predict_level_MULTI <- function(index, output, y.level.est, n.ahead, new.data=NULL, confidence.level, sea.on=FALSE, return.diff=FALSE) {
  
  # Description:
  # Function to predict the level of a cumulative variable
  
  # Parameters:
  # y.level.est: The cumulative variable (e.g., total cases in an epidemic)
  # n.ahead: Number of periods ahead to forecast
  # confidence.level: Confidence level for the forecast interval
  # return.diff: If TRUE, returns the difference (incidence) instead of cumulative values
  
  # Return: xts object containing the point forecasts and upper and lower bounds of the forecast interval.
  
  model <- output$model
  n <- attr(model, 'n')
  p <- attr(model, 'p')
  
  # Make predication
  freq <- unclass(periodicity(y.level.est))$label
  endtime <- end(model$y) + c(0, n.ahead)
  filtered.out <- predict_all_MULTI(index, output, n.ahead, new.data=new.data, sea.on=sea.on, return.all=TRUE)
  
  # Extract parameters.
  timespan <- n + 0:n.ahead
  g.t <- exp(filtered.out$y.hat) # Calculate g.t as exponent of y.t
  v_dates_end <- seq(last(index(y.level.est)), last(index(filtered.out$y.hat)), by=freq) # Forecast dates
  
  # Construct CI
  Ptt <- filtered.out$P.t.t
  i.level <- grep('level', colnames(model$T))
  ci <- qnorm((1 - confidence.level) / 2) * sqrt(Ptt[i.level, i.level, timespan]) %o% c(1, -1)
  ci <- xts(exp(ci), order.by=v_dates_end)
  
  # Forecast loop (start by initializing the first cumulative value)
  y.hat <- xts(matrix(NA, nrow=n.ahead + 1, ncol=3), order.by=v_dates_end)
  y.hat[v_dates_end[1],] <- y.level.est[v_dates_end[1],1]
  
  # Extract growth rate and prediction intervals (recursive build cumulative curve)
  for (i in seq_len(length(v_dates_end[-1]))) {
    date.forecast <- v_dates_end[i + 1]
    date.lag <- date.forecast - 1
    # Update level
    y.hat[date.forecast, 1] <- as.numeric(y.hat[date.lag, 1]) * 
      as.numeric(1 + g.t[date.forecast,])
    # Get prediction intervals
    y.hat[date.forecast, 2:3] <- as.numeric(y.hat[date.lag, 1]) *
      as.numeric(1 + g.t[date.forecast,] %*% ci[date.forecast,])}
  
  # Difference y.level.est if requested
  d <- if (return.diff){diff(y.hat[, 1])[-1] } else {(y.hat[, 1])[-1]}
  ci_bounds <- if (return.diff) {
    (y.hat[, 2:3] - as.vector(y.hat[, 1]))[-1] + as.vector(d)
  } else { y.hat[2:dim(y.hat)[1], 2:3] }
  # Piece components together
  pred <- vector('list', length=p)
  pred[[p]] <- cbind(fit=d, prediction=ci_bounds)
  pred <- lapply(pred, ts, end=endtime, frequency=1)
  
  y.hat <- xts(pred[[p]], order.by=v_dates_end[-1])
  names(y.hat)[2:3] <- list('lower', 'upper')
  
  return(as.xts(y.hat))
  
}


if (period==1) {
  covar_all<-c('t2m_smooth', 'tp_smooth')
} else if (period==2){
  covar_all<-c('t2m_smooth', 'tp_smooth')
}


# Set identified parameters
q <- 0.005
n.ahead <- max_prediction_horizon
confidence.level <- 0.68
sea.type <- 'trigonometric'
sea.period <- 7
sea.on <- FALSE


data_rt <- data_rt %>%
  dplyr::group_by(date_analysis, region) %>%
  dplyr::mutate(across(iHosp_smooth, ~ cumsum(replace_na(., 1)), .names = '{.col}_c'))


# FILTER TO ANALYSIS DATE ----
res<-tibble()
for (iDate in seq_along(dates_analysis)){
  
  print(dates_analysis[iDate]) 
  dat_train<-data_rt %>%
    filter(date_analysis==dates_analysis[iDate]) %>%
    mutate(iHosp_smooth_c_covar=iHosp_smooth_c)
  
  # IDENTIFY COVARIATE LAGS -----
  res_lag<-tibble()
  list_col<-sort(covar_all)
  for (iCol in seq_along(list_col)){
    
    mycol<-list_col[iCol]
    
    for (lag in 1:12){
      
      sub1<-dat_train
      sub2<-sub1
      sub2$date<-sub2$date+lag
      
      sub<-left_join(sub1[,c("date","region","iHosp_smooth_c")],sub2[,c("date","region",mycol)], by = c("date", "region"))
      sub<-sub %>% drop_na(iHosp_smooth_c,mycol) %>%
        as.data.frame()
      
      cor<-cor.test(y=sub$iHosp_smooth_c,x=sub[,mycol],method="pearson")$estimate
      
      res_lag_tmp<-data.frame(pearson=cor,region="regions",lag=lag,var=mycol)
      res_lag<-bind_rows(res_lag,res_lag_tmp)
      
    }
  }
  
  best_lag<-res_lag %>% 
    group_by(var) %>%
    slice(which.max(abs(pearson)))
  
  if ("iHosp_smooth_c_covar" %in% list_col){
    best_lag$lag[best_lag$var=="iHosp_smooth_c_covar"]<-1
  }
  
  list_lag<-best_lag$lag

  # FIT MODEL FOR EACH REGION -----
  for (iRegion in seq_along(regions)) {
    
    myreg<-regions[iRegion]
    print(myreg) 
    
    dat_train_reg<-dat_train %>%
      filter(region==myreg)
    
    # Prepare lagged data
    sub<-dat_train_reg[,c("date","region","iHosp_smooth_c")]
    sub1<-dat_train_reg[,c("date","region","iHosp_smooth_c",list_col)]
    
    for (i in 1:length(list_lag)){
      
      sub2<-sub1
      sub2$date<-sub2$date+list_lag[i]
      
      sub<-left_join(sub,sub2[,c("date","region",list_col[i])],by=c("region"="region","date"="date"))
      
    }
    
    sub<-sub %>% drop_na(iHosp_smooth_c,all_of(list_col))
    
    #  ----- SET PARAMETERS
    
    # Create xts object and subset relevant variables as numeric
    Y <- xts(subset(sub, select=c('iHosp_smooth_c', covar_all)),  order.by=sub$date)
    # Get growth rate
    y.lg <- tsgc::df2ldl(Y[,'iHosp_smooth_c'])
    
    #  ----- FIT INITIAL MODEL & UPDATE PARAMETER ESTIMATES -----
    
    # Set state-space model
    y <- as.numeric(y.lg)
    lapply(covar_all, function(exp) assign(exp, as.numeric(Y[, exp]), envir=.GlobalEnv))
    formula <- paste('~', paste(covar_all, collapse='+'))
    if (is.null(covar_all)){
      model <- SSModel(y ~ 
                         SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                         SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))
    } else {
      model <- SSModel(y ~ SSMregression(as.formula(formula), Q=diag(0,length(covar_all))) + #Q=diag(0,length(c(exposures, controls)))
                         SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                         SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))
    }
    
    # Partial update of known params (seasonal, observational & slope noise)
    updatefn <- purrr::partial(update, ... =, q=q, sea.type=sea.type)
    
    # Full update, using MLE for unknown model params
    model_fit <- fitSSM(model, inits=c(0,0), updatefn=updatefn, method='BFGS')
    
    # Kalman filter
    output <- KFS(model_fit$model)
    
    # Get truncated index from model (if using a re-initialisation, self$get_model)
    index <- if (!is.null(model$index)) { model$index } else { index(y.lg) }
    
    #  ----- BULD NEW DATA FOR PREDICTIONS -----
    newdata<-data.frame(seq.Date(dates_analysis[iDate]-remove_last_n+1,dates_analysis[iDate]+14,by="day"),NA)
    names(newdata)<-c("date","iHosp_smooth_c")
    
    sub1<-dat_train_reg[,c("date","iHosp_smooth_c",list_col)]
    
    for (i in 1:length(list_lag)){
      
      sub2<-sub1
      sub2$date<-sub2$date+list_lag[i]
      
      newdata<-left_join(newdata,sub2[,c("date",list_col[i])],by=c("date"="date"))
      
      # Fill NA by last available value
      newdata[is.na( newdata[,ncol(newdata)]),ncol(newdata)]<- last(newdata[!is.na( newdata[,ncol(newdata)]),ncol(newdata)])
      
    }
    new.data <- xts(subset(newdata, select=covar_all),  order.by=newdata$date)
    
    
    # ----- PERCENTILE-BASED PREDICTIONS -----
    alpha <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)
    
    est.date.index <- index %>% as.Date()
    y.level.est    <- Y[est.date.index]
    
    # Allocate output array (point + lower/upper per alpha, same as GAM)
    pred <- array(data = NA, dim = c(n.ahead, length(alpha) * 2 + 1))
    
    for (iAlpha in seq_along(alpha)) {
      
      # alpha=10 -> 90% CI -> confidence.level=0.90
      conf_i <- 1 - alpha[iAlpha] / 100
      
      y.hat.diff <- predict_level_MULTI(
        index            = index,
        output           = output,
        y.level.est      = y.level.est,
        n.ahead          = n.ahead,
        new.data         = new.data,
        confidence.level = conf_i,
        return.diff      = TRUE,
        sea.on           = sea.on
      )
      
      # Point forecast is the same for all alpha — only store once
      if (iAlpha == 1) pred[, 1] <- pmax(0, round(as.numeric(y.hat.diff[, 'fit'])))
      
      pred[, 2 * iAlpha]     <- pmax(0, round(as.numeric(y.hat.diff[, 'lower'])))
      pred[, 2 * iAlpha + 1] <- pmax(0, round(as.numeric(y.hat.diff[, 'upper'])))
      
    }
    
    colnames(pred) <- c("point", paste0(c("lower_", "upper_"), rep(alpha, each = 2)))
    
    # ----- BIND RESULTS -----
    res_tmp <- data.frame(
      var                = "iHosp",
      region             = myreg,
      date               = dates_analysis[iDate],
      prediction_horizon = seq(1, n.ahead) - remove_last_n,
      pred
    )
    
    res <- rbind(res, res_tmp)
    
  }
  
}

res<-res %>% mutate(model=model) 


write.csv2(res,file=paste0(path_results,model,"_period",period,".csv"),row.names=F)