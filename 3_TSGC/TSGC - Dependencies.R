### ------------------------------------------------------------------------ ### 
### TSGC: Dependencies
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation and set seed
options(scipen=999)
set.seed(123)


#### PATHS ####
user <- 'claer14'
path_data <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/', user, '/Documents/GitHub/hisp-climate-health-forecast/')
path_project <- paste0(path_git, '3_TSGC/')
path_report <- paste0(path_git, '0_Reporting/TSGC/')

#### PACKAGES ####

# Packages
library(quantmod); library(KFAS); library(xts); library(abind);
library(ggplot2); library(ggthemes); library(zoo); library(utils); 
library(lubridate); library(tsgc); library(gridExtra); library(patchwork);
library(purrr); library(tidyr); library(dplyr); library(openxlsx); library(stringr)

# FUNCTIONS: MODEL ####

### Function to fit model with set parameters
get_results <- function(Y, y, exposures, controls, y.lg, y.cum, y.eval, new.data,
                        n.ahead, q, confidence.level, sea.period, sea.type, 
                        plt.start.date){
  
  #  !!!!!!!!!!!! FIT MODEL & UPDATE PARAMETER ESTIMATES  !!!!!!!!!!!! 
  # Format y & x variables for model (dynamic adds exposures)
  y <- as.numeric(y.lg)
  lapply(c(exposures, controls), function(exp) assign(exp, as.numeric(y.cum[, exp]), envir=.GlobalEnv))
  formula <- paste('~', paste(c(exposures, controls), collapse='+'))
  
  # Set state-space model
  if (is.null(c(exposures, controls))){
    model <- SSModel(y ~ 
                       SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                       SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))
  } else {
    model <- SSModel(y ~ SSMregression(as.formula(formula), Q=diag(0,length(c(exposures, controls)))) + #Q=diag(0,length(c(exposures, controls)))
                       SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                       SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))
  }
  
  # Partial update of known params (seasonal, observational & slope noise)
  updatefn <- purrr::partial(update, ... =, q=q, sea.type=sea.type)
  
  # Full update, using MLE for unknown model params
  model_fit <- fitSSM(model, inits=c(0,0), updatefn=updatefn, method='BFGS')
  
  # Kalman filter
  (output <- KFS(model_fit$model))
  # diff(output$model$y)
  # mean(output$alphahat[, 1] )
  
  # Get truncated index from model (if using a re-initialisation, self$get_model)
  index <- if (!is.null(model$index)) { model$index } else { index(y.lg) }
  
  # !!!!!!!!!!!! VISUALS !!!!!!!!!!!! 
  # Plot growth rate forecast
  p1 <- plot_forecast_MULTI(output=output, index=index, y.eval=y.eval, n.ahead=n.ahead,
                            plt.start.date=plt.start.date)
  
  # Plot new cases forecast
  p2 <- plot_new_cases_MULTI(output=output, index=index, Y=y.cum, n.ahead=n.ahead,
                             confidence.level=confidence.level, date_format='%Y-%m-%d',
                             plt.start.date=plt.start.date, og=FALSE)
  
  # Plot forecast accuracy (relative to forecast, not trend)
  p3 <- plot_holdout_MULTI(output=output, index=index, Y=y.cum, Y.eval=Y[(estimation.date.end+0:n.ahead)],
                           confidence.level=confidence.level, date_format='%Y-%m-%d',og=FALSE, MAPE=NULL)
  mape.sea <- p3$MAPE
  mape.trend <- p3$MAPE.TREND
  p3 <- p3$figure
  
  # !!!!!!!!!!!! COEFFICIENTS !!!!!!!!!!!! 
  n <- dim(output$alphahat)[1] # number of time points
  estimates <- output$alphahat[n, ] # smoothed state estimates at time n
  std_errors <- sqrt(diag(output$V[, ,n])) # get standard errors (square root of diag covariance matrix of n)
  # Lower and upper bounds
  z_star <- qnorm(0.95)
  lower <- estimates - z_star * std_errors
  upper <- estimates + z_star * std_errors
  res <- data.frame(State = names(estimates), 
                    Estimate = as.numeric(estimates), 
                    CI_Lower = lower,
                    CI_Upper = upper,
                    Std_Error = as.numeric(std_errors))
  rownames(res) <- NULL
  
  # !!!!!!!!!!!! CRITERION !!!!!!!!!!!! 
  # NB: length(model_fit$optim.out$par) is the number of unique NAs KFAS  treats as independent hyper-parameters.
  mod <- model_fit$model            
  LL  <- logLik(mod) # higher is better                
  AIC <- -2*LL + 2*length(model_fit$optim.out$par) # lower is better
  BIC <- -2*LL + log(nrow(y.lg))*length(model_fit$optim.out$par) # lower is better
  
  # Create the data frame
  comp <- data.frame(Model = plot_title, MAPE=paste0(mape.sea, '%'), LogLik = LL, AIC = AIC, BIC = BIC)
  
  
  combined_plot <- (p1 | p2 | p3) +
    plot_annotation(
      title = plot_title,
      subtitle = paste0('MAPE: ', mape.sea, '%', sep=''),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 10)),  # Adjust vertical margin to move title down
        plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = -10)) # Adjust subtitle margin
      )) & theme(plot.margin = margin(0, 0, 0, 0))
  
  p3 <- plot_holdout_MULTI(output=output, index=index, Y=y.cum, Y.eval=Y[(estimation.date.end+0:n.ahead)],
                           confidence.level=confidence.level, date_format='%Y-%m-%d',
                           title=NULL, og=FALSE, MAPE=TRUE)
  
  return(list(model=output, estimates=res, criteria=comp, plots=combined_plot,
              p1=p1, p2=p2, p3=p3))
  
}

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

# FUNCTIONS: PLOTTING ####

### Function for plotting data with moving average
plot_moving_avg <- function(variable, MA=7, y.lab='', colour='red', plt.start.date=NULL, est.date.index, estimation.date.end){
  
  # Determine plot length
  if (is.null(plt.start.date)) {plt.start.date <- head(est.date.index, 1)}
  tmp.date <- min(estimation.date.end - 4, plt.start.date)
  s <- sprintf('%s/', format(tmp.date, '%Y-%m-%d'))
  variable <- variable[s]
  
  # Calculate the 7 day centered moving average
  ma.cent <- zoo::rollmean(variable, MA, align='center') 
  
  # Create df for plot (Data data vs. MA data)
  d <- setNames(cbind(variable, ma.cent), c('New Cases', 'Centered 7-day MA'))
  df.f <- data.frame(Date=index(d), 
                     New.Cases=coredata(d[,1]),
                     Centered.7.day.MA=coredata(d[,2]))
  
  
  # Plot Figure 1
  p1 <- ggplot(df.f, aes(x = Date)) +
    geom_line(aes(y = New.Cases, color = 'Real data'), linewidth = 0.1) +
    geom_line(aes(y = Centered.7.day.MA, color = 'Centered 7 day MA'), linewidth = 1) +
    scale_color_manual(name = '', values = c('Real data' = 'grey', 'Centered 7 day MA' = colour)) +
    scale_y_continuous(n.breaks = 6) +
    scale_x_date(labels = scales::date_format('%d %b %y'), 
                 breaks = seq(min(df.f$Date), max(df.f$Date), by = '3 months')) +
    xlab('Day') + 
    ylab(y.lab) +
    theme_economist_white(gray_bg = FALSE, base_size = 12) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, vjust = 0.75),
      axis.title = element_text(size = rel(1), margin = margin(t = 10, r = 10)),
      axis.text = element_text(size = rel(1)),
      plot.title = element_text(face = 'bold', margin = margin(b = 5)),
      plot.caption = element_text(size = rel(1)),
      text = element_text(size = rel(1.1))
    ) +
    scale_linetype_manual(values = rep('solid', 3)) +
    scale_size_manual(values = rep(1, 3))
  
  return(p1) 
}

### Function plotting Data and filtered values of the log cumulative growth rate
plot_forecast_MULTI <- function(index, output, y.eval, n.ahead=14, plt.start.date=NULL, title='', caption='') {
  
  # Description: Function plots Data and filtered values of the log cumulative growth rate 
  # (ln(g_t)) within the estimation sample, along with forecasted and realised 
  # values outside the estimation sample.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - y.eval: Out-of-sample realisation of the log growth rate of the cumulated variable 
  #           (Data values for comparison with forecasts).
  # - n.ahead: Number of time periods ahead to forecast from the end of the sample 
  #            (default is 14).
  # - plt.start.date: Start date for the plot (default: NULL, meaning the start of 
  #                   the estimation sample).
  # - title: Title for the plot (entered as a text string).
  # - caption: Caption for the plot (entered as a text string).
  #
  # Returns: A ggplot2 plot visualizing the forecast and realised values of the 
  # log cumulative growth rate.
  
  Date <- NULL
  model <- output$model
  est.date.index <- index
  
  y <- xts::xts(output$model$y %>% as.numeric(), order.by=est.date.index)
  p <- attr(output$model, 'p')
  
  y.hat.all <- predict_all_MULTI(index, output, n.ahead, new.data, sea.on=FALSE, return.all=TRUE)
  
  y.pred <-  subset(y.hat.all$y.hat,index(y.hat.all$y.hat) > tail(index,1))
  filtered.level <- y.hat.all$level
  
  if (p == 1) {
    EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
    d <- cbind(y, filtered.level, y.pred,
               y.eval[index(y.eval)>tail(index(est.date.index),1),])
    if (!is.null(plt.start.date)) { d <- d[index(d) >= plt.start.date] }
    d <- d[index(d) <= tail(index(y.pred),1)]
    names(d) <- c(
      'EstimationSample', 'FilteredLevel', 'Forecast', 'RealisedData'
    )
    
    df_plot <- as.data.frame(d)
    df_plot$Date <- as.Date(rownames(df_plot))
    
    p1 <- ggplot2::ggplot(data=df_plot, aes(x=Date)) +
      ggplot2::geom_line(aes(y=EstimationSample, color='Data'), lwd=0.5) +
      ggplot2::geom_line(aes(y=FilteredLevel, color='Filtered\nTrend'), lwd=0.85) +
      ggplot2::geom_line(aes(y=Forecast, color='Forecast\nTrend'), lwd=0.85) +
      ggplot2::geom_line(aes(y=RealisedData, color='Forecast'), lwd=0.75) +
      ggplot2::scale_color_manual(
        values=c('Data'=1, 'Filtered\nTrend'=2, 'Forecast\nTrend'=3,'Forecast'='grey'),
        limits=c('Data', 'Filtered\nTrend','Forecast\nTrend', 'Forecast')) +
      scale_linetype_manual(values=c('solid', 'solid', 'solid', 'dashed')) +
      labs(x='Date', y='Log Growth Rate', caption=caption, title=title) +
      scale_x_date(labels=scales::date_format('%d %b %y'), breaks=seq(min(df_plot$Date), max(df_plot$Date), by='2 weeks')) +
      theme_economist_white(gray_bg = FALSE, base_size = 12) +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        text = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        axis.title.y = element_text(size = rel(1), margin = margin(r = 10)),
        axis.title.x = element_text(size = rel(1), margin = margin(t = 10)),
        plot.title = element_text(margin = margin(b = 5)),
        plot.caption = element_text(size = rel(1))
      )
  } else if (p == 2) {
    g_1 <- g_2 <- delta <- Forecast <- RealisedData <- NULL
    d <- cbind(y, filtered.level, y.pred[,2],
               y.eval[index(y.eval)>tail(index(est.date.index),1),2])
    d <- d[index(d) <= tail(index(y.pred),1)]
    names(d) <- c('g_1', 'g_2', 'delta', 'Forecast', 'RealisedData')
    
    df_plot <- as.data.frame(d)
    df_plot$Date <- as.Date(rownames(df_plot))
    
    p1 <- ggplot2::ggplot(data=df_plot, aes(x=Date)) +
      ggplot2::geom_line(aes(y=g_1, color='g_1')) +
      ggplot2::geom_line(aes(y=g_2, color='g_2')) +
      ggplot2::geom_line(aes(y=g_2, color='delta')) +
      ggplot2::geom_line(aes(y=Forecast, color='Forecast')) +
      ggplot2::geom_line(aes(y=RealisedData, color='Realised\nData')) +
      ggplot2::scale_color_manual(
        values=c(1, 2, 3, 4, 'grey')) +
      ggplot2::scale_linetype_manual(values=c('solid', 'solid', 'solid', 'solid', 'dashed')) +
      labs(x='Date', y='Log Growth Rate', caption=caption, title=title) +
      scale_x_date(labels=scales::date_format('%d %b %y'), breaks=seq(min(df_plot$Date), max(df_plot$Date), by='2 weeks')) +
      theme_economist_white(gray_bg = FALSE, base_size = 12) +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        text = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        axis.title.y = element_text(size = rel(1), margin = margin(r = 10)),
        axis.title.x = element_text(size = rel(1), margin = margin(t = 10)),
        plot.title = element_text(margin = margin(b = 5)),
        plot.caption = element_text(size = rel(1))
      )
  } else { stop('NotImplemented Error') }
  
  return(p1)
}

### Function plotting new cases
plot_new_cases_MULTI <- function(index, output, Y, n.ahead, ylim=NULL, confidence.level=0.68, date_format='%Y-%m-%d', title=NULL, plt.start.date=NULL, og=FALSE) {
  
  # Description: Function plots the Data and forecasted values of new cases 
  # (i.e., the difference in the cumulated variable). It includes forecasts with 
  # and without seasonal adjustments and shows forecast intervals.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - Y: Cumulated variable (time series data).
  # - n.ahead: Number of future periods to forecast.
  # - confidence.level: Confidence level for the prediction intervals (default is 0.68, approximately one standard deviation).
  # - date_format: Format for date labels on the plot (default: '%Y-%m-%d').
  # - title: Title for the plot (default: NULL, meaning no title).
  # - plt.start.date: Start date for the Data data shown in the plot 
  #
  # Returns: A ggplot2 plot visualizing the forecast.
  
  Date <- Data <- Forecast <- ForecastTrend <- lower <- upper <- NULL
  if (is.null(title)) {title <- ''}
  est.date.index <- index %>% as.Date()
  estimation.date.end <- tail(est.date.index, 1)
  y.level.est <- Y[est.date.index]
  if (is.null(plt.start.date)) {plt.start.date <- head(est.date.index, 1)}
  
  # Get the CIs for the trend (y.hat.diff.final.ci) AND forecast (y.hat.diff.final.ci)
  y.hat.diff.final.ci <- predict_level_MULTI(index, output, y.cum=y.level.est, 
                                             n.ahead, new.data=new.data, 
                                             confidence.level, sea.on=FALSE, 
                                             return.diff=TRUE)
  
  y.hat.diff.final <- predict_level_MULTI(index, output, y.cum=y.level.est, 
                                          n.ahead, new.data=new.data, 
                                          confidence.level, sea.on=TRUE, 
                                          return.diff=TRUE)
  
  tmp.date <- min(estimation.date.end - 4, plt.start.date)
  s <- sprintf('%s/', format(tmp.date, '%Y-%m-%d'))
  d.plot <- cbind(
    diff(y.level.est[,1])[s],
    # y.level.est[,2][s], # this adds other expalantory variables - but not necessary is it?
    y.hat.diff.final[, 1],
    y.hat.diff.final.ci[, 1]
  )
  names(d.plot) <- c('Data', 'Forecast', 'ForecastTrend') #'Expl',
  
  # !!! Original TSGC package looks at CI of TREND (y.hat.diff.final.ci) !!!
  # !!! We are interested in CI of the FORECAST (y.hat.diff.final) !!!
  if (og){
    ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                              (as.Date(index(y.hat.diff.final.ci),
                                       format=date_format))))
  } else{
    ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final[, 2:3]),
                              (as.Date(index(y.hat.diff.final),
                                       format=date_format))))
  }
  colnames(ci) <- c('lower', 'upper', 'date')
  ci[, 'date'] <- as.Date(ci[, 'date'], format=date_format, origin='1970-01-01')
  
  df_plot <- as.data.frame(d.plot)
  df_plot$Date <- as.Date(rownames(df_plot), format=date_format)
  
  p1 <- ggplot2::ggplot(data=df_plot, aes(x=Date)) +
    ggplot2::geom_line(aes(y=Data, color='Data'), lwd=0.5) +
    ggplot2::geom_line(aes(y=Forecast, color='Forecast'), lwd=0.85) +
    # ggplot2::geom_line(aes(y=ForecastTrend, color='Forecast\nTrend'),lwd=0.85) +
    ggplot2::scale_color_manual(values=c('black', '#AA2045')) + # 'lightgreen'
    ggplot2::geom_ribbon(data=ci, aes(x=date, ymin=lower, ymax=upper),
                         linetype=0, linewidth=0, fill='#AA2045', alpha=0.1) +
    labs(x='Date', y='Patients', title=title) +
    scale_x_date(labels=scales::date_format('%d %b %y'), breaks=seq(min(df_plot$Date), max(df_plot$Date), by='2 weeks')) +
    ggplot2::scale_linetype_manual(values=c('solid', 'solid', 'solid')) +
    ggplot2::scale_size_manual(values=c(1, 1, 1)) +
    theme_economist_white(gray_bg=FALSE, base_size=12) +
    theme(
      legend.title=element_blank(),
      legend.text = element_text(size = 10),
      text=element_text(size=rel(1.1)),
      axis.text=element_text(size=rel(1)),
      axis.text.x = element_text(angle = 45, vjust = 0.75),
      axis.title.y=element_text(size=rel(1),margin=margin(r=10)),
      axis.title.x=element_text(size=rel(1),margin=margin(t=10)),
      plot.title=element_text(margin=margin(b=5)),
      plot.caption=element_text(size=rel(1))
    )
  
  if(!is.null(ylim)){
    p1 <- p1 + coord_cartesian(ylim=c(0, ylim))
  }
  
  return(p1)
}

# Function plotting percentage prediction error
plot_holdout_MULTI <- function(index, output, Y, Y.eval, ylim=NULL, confidence.level=0.68, date_format='%Y-%m-%d', title= NULL, caption=NULL, og=FALSE, MAPE=NULL) {
  
  # Description: Function plots the forecast of new cases (difference in the cumulated variable) 
  # over a holdout sample. It includes Data values, forecasted values (with and 
  # without seasonal components), forecast intervals, and Data outcomes from 
  # the holdout sample. Additionally, it reports the mean absolute percentage 
  # prediction error over the holdout sample.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - Y: Values of the cumulated variable used in the estimation window.
  # - Y.eval: Values of the cumulated variable in the holdout sample 
  #           (Data values for comparison with forecasts).
  # - confidence.level: Width of the prediction interval for ln(g_t) used in 
  #                     forecasts of y_t=ΔY_t (default: 0.68, approximately 
  #                     one standard deviation).
  # - series.name: Name of the variable being forecasted, used for the y-axis 
  #                label (e.g., 'Cases' results in 'New Cases' on the y-axis).
  # - date_format: Format for date labels on the plot (default: '%Y-%m-%d').
  # - title: Title for the forecast plot (entered as a text string, default: NULL, 
  #          meaning no title).
  # - caption: Caption for the forecast plot (entered as a text string, default: NULL, 
  #            meaning no caption).
  #
  # Returns: A ggplot2 plot visualizing the forecast, Data values, and prediction 
  # intervals over the holdout sample.
  
  
  Date <- Data <- Forecast <- ForecastTrend <- lower <- upper <- NULL
  
  # Y.eval <- Y[(estimation.date.end+0:n.ahead)]
  # model <- output$model
  
  # Get sample period
  y.level.est <- Y[index]
  p <- attr(output$model, 'p')
  if(p!=1) { stop('NotImplementedError') }
  
  n.ahead <- tail(index(Y.eval),1)-tail(index(y.level.est),1) ## !!!!!!
  y.eval.diff <- diff(Y.eval) %>% na.omit
  
  est.date.index <- index %>% as.Date()
  estimation.date.end <- tail(index %>% as.Date(), 1)
  
  y.hat.diff.final.ci <- predict_level_MULTI(index, output, y.cum=y.level.est, n.ahead,
                                             new.data, confidence.level, return.diff=TRUE)
  
  y.hat.diff.final <- predict_level_MULTI(index, output, y.cum=y.level.est, n.ahead,
                                          new.data, confidence.level, return.diff=TRUE, sea.on=TRUE)
  
  d <- cbind(
    y.eval.diff[index(y.eval.diff)>estimation.date.end,],
    y.hat.diff.final[, 1],
    y.hat.diff.final.ci[, 1]
  )
  colnames(d)[1] <- 'Data'
  colnames(d)[(ncol(d)-1):ncol(d)] <- c('Forecast', 'ForecastTrend')
  
  df_plot <- as.data.frame(d)
  df_plot$Date <- as.Date(rownames(df_plot), format=date_format)
  
  d.eval <- na.omit(d)
  mape.trend <- 100*(abs(d.eval$Data - d.eval$`ForecastTrend`) / d.eval$Data) %>% mean %>% round(2)
  mape.sea <- 100*(abs(d.eval$Data - d.eval$Forecast) / d.eval$Data) %>% mean %>% round(2)
  
  # !!! Original TSGC package looks at CI of TREND (y.hat.diff.final.ci) !!!
  # !!! We are interested in CI of the FORECAST (y.hat.diff.final) !!!
  if (og){
    ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                              (as.Date(index(y.hat.diff.final.ci),
                                       format=date_format))))
  } else{
    ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final[, 2:3]),
                              (as.Date(index(y.hat.diff.final),
                                       format=date_format))))
  }
  colnames(ci) <- c('lower', 'upper', 'date')
  ci[, 'date'] <- as.Date(ci[, 'date'], format=date_format, origin='1970-01-01')
  
  p1 <- ggplot2::ggplot(data=df_plot, aes(x=Date)) +
    ggplot2::geom_line(aes(y=Data, color='Data'),lwd=0.75) +
    ggplot2::geom_line(aes(y=Forecast, color='Forecast'),lwd=0.85) +
    # ggplot2::geom_line(aes(y=ForecastTrend, color='Forecast\nTrend'),lwd=0.85) +
    ggplot2::scale_color_manual(values=c('black', '#AA2045')) + # 'lightgreen'
    ggplot2::geom_ribbon(data=ci, aes(x=date, ymin=lower, ymax=upper),
                         linetype=0, linewidth=0, fill='#AA2045',
                         alpha=0.1) +
    labs(x='Date', y='Patients', title=title) + # subtitle=paste('MAPE: ', mape.sea, '%. Trend MAPE: ', mape.trend, '%.',sep='')
    scale_x_date(labels=scales::date_format('%d %b %y'), breaks=seq(min(df_plot$Date), max(df_plot$Date), by='4 days')) +
    scale_linetype_manual(values=c('solid', 'solid', 'solid')) +
    scale_size_manual(values=c(1, 1.5, 1)) +
    theme_economist_white(gray_bg=FALSE, base_size=14) +
    theme(
      legend.title=element_blank(),
      legend.text = element_text(size = 10),
      text=element_text(size=rel(1)),
      axis.text=element_text(size=rel(1)),
      axis.text.x = element_text(angle = 45, vjust = 0.75),
      axis.title.y=element_text(size=rel(1), margin=margin(r=10)),
      axis.title.x=element_text(size=rel(1), margin=margin(t=10)),
      plot.title=element_text(margin=margin(b=5)),
      plot.subtitle=element_text(size=rel(1), hjust=0.5,  margin=margin(t=3)))
  
  if(!is.null(ylim)){
    p1 <- p1 + coord_cartesian(ylim=c(0, ylim))
  }
  
  
  if (!is.null(MAPE)){
    p1 <- p1 + labs(x='Date', y='Patients', title=title, subtitle=paste('MAPE: ', mape.sea, '%')
    ) 
    
  }
  
  return(list(figure=p1, MAPE=mape.sea, MAPE.TREND=mape.trend))
}


# FUNCTIONS: FORMATTING ####

# Function to replace variable names with set names
replace_exact_matches <- function(df, dict) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      col <- ifelse(col %in% names(dict), dict[col], col)
    }
    col
  })
  df
}

# Function that rounds all numeric data in df
round_numeric <- function(df, digit) {
  df[] <- lapply(df, function(x) if(is.numeric(x)) round(x, digit) else x)
  df
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
    new <- SSModel(rep(NA, n.ahead) ~ SSMregression(~ new.data, Q=diag(0, length(c(exposures, controls)))) + #Q=diag(0,length(c(exposures, controls)))
                     SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))) +
                     SSMseasonal(period=sea.period, Q=NA, sea.type=sea.type), H=matrix(NA))
    
  }
  # Get n.ahead forecast with new fit
  updatefn <- purrr::partial(update, ... =, q=q, sea.type=sea.type)
  new_fit <- fitSSM(new, inits=c(0, 0), updatefn=updatefn, method='BFGS')
  y.hat.kfas <- if (sea.on) {
    predict(output$model, interval='prediction', newdata=new_fit$model, level=0.68, states='all')
  } else {
    predict(output$model, interval='prediction', newdata=new_fit$model, level=0.68, states='level')
  }
  
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
predict_level_MULTI <- function(index, output, y.cum, n.ahead, new.data=NULL, confidence.level, sea.on=FALSE, return.diff=FALSE) {
  
  # Description:
  # Function to predict the level of a cumulative variable
  
  # Parameters:
  # y.cum: The cumulative variable (e.g., total cases in an epidemic)
  # n.ahead: Number of periods ahead to forecast
  # confidence.level: Confidence level for the forecast interval
  # return.diff: If TRUE, returns the difference (incidence) instead of cumulative values
  
  # Return: xts object containing the point forecasts and upper and lower bounds of the forecast interval.
  
  model <- output$model
  n <- attr(model, 'n')
  p <- attr(model, 'p')
  
  # Make predication
  freq <- unclass(periodicity(y.cum))$label
  endtime <- end(model$y) + c(0, n.ahead)
  filtered.out <- predict_all_MULTI(index, output, n.ahead, new.data=new.data, sea.on=sea.on, return.all=TRUE)
  
  # Extract parameters.
  timespan <- n + 0:n.ahead
  g.t <- exp(filtered.out$y.hat) # Calculate g.t as exponent of y.t
  v_dates_end <- seq(last(index(y.cum)), last(index(filtered.out$y.hat)), by=freq) # Forecast dates
  
  # Construct CI
  Ptt <- filtered.out$P.t.t
  i.level <- grep('level', colnames(model$T))
  ci <- qnorm((1 - confidence.level) / 2) * sqrt(Ptt[i.level, i.level, timespan]) %o% c(1, -1)
  ci <- xts(exp(ci), order.by=v_dates_end)
  
  # Forecast loop (start by initializing the first cumulative value)
  y.hat <- xts(matrix(NA, nrow=n.ahead + 1, ncol=3), order.by=v_dates_end)
  y.hat[v_dates_end[1],] <- y.cum[v_dates_end[1],1]
  
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
  
  # Difference output if requested
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


#### DEPENDENCIES ####
mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9')

nm <- c(total_patients='Total patients', M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        general='General medicine', emergency='Emergency', ophthalmology='Ophthalmology', 
        paediatrics='Paediatric', orthopedics='Orthopedics', ent='ENT', 
        gynecology='Gynecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        mx2t_IMD='IMD: Mx2t >30\u00B0C', mn2t_IMD='IMD: Mn2t <0\u00B0C', tp_IMD='IMD: Tp >15mm', fg10_IMD='IMD: Fg10 >40km/h',
        mx2t_p90='Clim: Mx2t >90th', mn2t_p10='Clim: Mn2t <10th', tp_p90='Clim: Tp >90th', fg10_p90='Clim: Fg10 >90th',
        mx2t='Mx2t (\u00B0C)', mn2t='Mn2t (\u00B0C)', tp='Tp (mm)', fg10='Fg10 (km/h)', 
        closed='Closed', holiday='Holiday') 

weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
h_cols <- c('total_patients', 'F', 'M', 'y1_l', 'y1_17', 'y18_60', 'y60_g', 
            'general_medicine', 'ophthalmology', 'pediatrics', 'orthopedics', 
            'ent', 'gynecology', 'dermatology', 'surgery', 'dentistry', 
            'emergency', 'psychiatry', 'oncology')



#############################################################################
# ### Function returning table of Kalman estimated variance parameters
# print_estimation_results <- function(index, output) {
#   
#   # Description:
#   # Prints a formatted table of kalman estimated variance parameters
#   # (Optional - ready to paste into LaTeX)
#   
#   # Parameters:
#   # - None
#   
#   # Returns: Table
#   
#   H <- output$model$H[, , 1]
#   Q_gamma <- output$model$Q[2, 2, 1]
#   Q_seasonal <- output$model$Q[3, 3, 1]
#   
#   tbl <- data.frame(
#     a=format(H, digits=3),
#     b=format(Q_gamma, digits=3),
#     c=format(Q_seasonal, digits=3),
#     d=format(Q_gamma / H, digits=4))
#   header.names <- c('$\\sigma_\\varepsilon^2$',
#                     '$\\sigma_\\gamma^2$',
#                     '$\\sigma_{seas}^2$',
#                     'q')
#   
#   out <- tbl %>%
#     kableExtra::kbl(
#       caption='Estimated parameters',
#       col.names=header.names,
#       # format='latex',
#       booktabs=TRUE,
#       escape=FALSE
#     ) %>%
#     kableExtra::kable_classic(full_width=FALSE, html_font='Cambria') #%>%
#   # kableExtra::footnote(general=' ')
#   
#   return(out)
# }
# 
# ### Function extracting smoothed/filtered growth rate of cases and associated components (δ and γ)
# get_growth_y <- function(index, output, smoothed=FALSE, return.components=FALSE) {
#   
#   # Description:
#   # Returns the growth rate of the incidence (`y`) of the cumulated variable (`Y`).
#   # Computed as: g_t=exp(δ_t) + γ_t
#   
#   # Parameters:
#   # - `smoothed`: Logical; whether to use smoothed estimates of δ and γ (`TRUE`) 
#   #   or contemporaneous filtered estimates (`FALSE`, default).
#   # - `return.components`: Logical; whether to return δ, γ, and growth rate estimates (`TRUE`), 
#   #   or just the growth rate (`FALSE`, default).
#   
#   # Return Value: An `xts` object containing smoothed/filtered growth rates and 
#   # components (δ and γ), if applicable.
#   
#   kfs_out <- output
#   idx <- index
#   
#   if (smoothed) {
#     att <- kfs_out$alphahat
#   } else {
#     att <- kfs_out$att
#   }
#   
#   filtered.slope <- xts(att[, 'slope'], order.by=idx)
#   filtered.level <- xts(att[, 'level'], order.by=idx)
#   g.t <- exp(filtered.level)
#   gy.t <- g.t + filtered.slope
#   names(gy.t) <- if (smoothed) { 'smoothed gy.t' } else { 'filtered gy.t' }
#   names(g.t) <- if (smoothed) { 'smoothed g.t' } else { 'filtered g.t' }
#   names(filtered.slope) <- if (smoothed) { 'smoothed gamma.t' } else {'filtered gamma.t' }
#   if (return.components) {
#     return(list(gy.t, g.t, filtered.slope))
#   } else {
#     return(gy.t)
#   }
# }
# 
# ### Function getting CIs of smoothed/filtered growth rate of cases
# get_gy_ci <- function(index, output, smoothed=FALSE, confidence_level=0.68) {
#   
#   # Description:
#   # Returns the growth rate of the incidence (`y`) of the cumulated variable (`Y`).
#   # Computed as: g_t=exp(δ_t) + γ_t
#   
#   # Parameters:
#   # - `smoothed`: Logical; whether to use smoothed estimates of δ and γ (`TRUE`) 
#   #   or contemporaneous filtered estimates (`FALSE`, default).
#   # - `confidence_level`: Confidence level for the confidence interval. Default is 0.68 
#   #   (one standard deviation for a normally distributed random variable).
#   
#   # Return Value:
#   # - An `xts` object containing smoothed/filtered growth rates and the upper and 
#   #   lower bounds for the confidence intervals.
#   
#   kfs_out <- output
#   idx <- index
#   
#   if (smoothed) {
#     att <- kfs_out$alphahat
#   } else {
#     att <- kfs_out$att
#   }
#   
#   filtered_slope <- xts(att[, 'slope'], order.by=idx)
#   filtered.level <- xts(att[, 'level'], order.by=idx)
#   g.t <- exp(filtered.level)
#   gy.t <- g.t + filtered_slope
#   
#   idx.slope <- grep('slope', colnames(kfs_out$att))
#   ci <- qnorm((1 - confidence_level) / 2) * sqrt(kfs_out$Ptt[idx.slope, idx.slope,]) %o% c(1, -1)
#   ci_bounds <- as.vector(gy.t) + ci
#   
#   pred <- xts(cbind(gy.t, ci_bounds), order.by=idx)
#   colnames(pred) <- c('fit','lower','upper')
#   
#   return(pred)
# }
