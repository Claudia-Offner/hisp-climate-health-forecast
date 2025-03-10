# Created by: Craig Thamotheram
# Created on: 11/02/2022

# ------------------------------------------------------------------------------
# -------------------------------- DESCRIPTION ---------------------------------
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggthemes)
library(zoo)
library(utils)

# ------------------------------------------------------------------------------
# --------------------------------- MAIN CODE ---------------------------------- 
# ------------------------------------------------------------------------------

# Function plotting new cases
plot_new_cases <- function(res, Y, n.ahead, confidence.level = 0.68, date_format = "%Y-%m-%d", title=NULL, plt.start.date=NULL) {
  
  # Description: Function plots the actual and forecasted values of new cases 
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
  # - plt.start.date: Start date for the actual data shown in the plot 
  #
  # Returns: A ggplot2 plot visualizing the forecast.
  
  Date <- Data <- Forecast <- ForecastTrend <- lower <- upper <- NULL
  if (is.null(title)) {title <- ""}
  est.date.index <- res$index %>% as.Date()
  estimation.date.end <- tail(est.date.index, 1)
  y.level.est <- Y[est.date.index]
  if (is.null(plt.start.date)) {plt.start.date <- head(est.date.index, 1)}
  
  y.hat.diff.final.ci <- res$predict_level(
    y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
    return.diff = TRUE
  )
  y.hat.diff.final <- res$predict_level(
    y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
    sea.on = TRUE, return.diff = TRUE
  )
  
  tmp.date <- min(estimation.date.end - 4, plt.start.date)
  s <- sprintf("%s/", format(tmp.date, "%Y-%m-%d"))
  d.plot <- cbind(
    diff(y.level.est)[s],
    y.hat.diff.final[, 1],
    y.hat.diff.final.ci[, 1]
  )
  names(d.plot) <- c('Data', 'Forecast', 'ForecastTrend')
  
  ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                            (as.Date(index(y.hat.diff.final.ci),
                                     format = date_format))))
  colnames(ci) <- c('lower', 'upper', 'date')
  ci[, 'date'] <- as.Date(
    ci[, 'date'], format = date_format, origin = "1970-01-01"
  )
  
  df_plot <- as.data.frame(d.plot)
  df_plot$Date <- as.Date(rownames(df_plot), format = date_format)
  
  ggplot2::ggplot(data = df_plot, aes(x = Date)) +
    ggplot2::geom_line(aes(y = Data, color = "Data"), lwd = 0.85) +
    ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
    ggplot2::geom_line(
      aes(y = ForecastTrend, color = "Forecast\nTrend"), lwd = 0.85
    ) +
    ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, aes(x = date, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    labs(x = "Date", y = "New Cases", title = title) +
    theme_economist_white(gray_bg = FALSE, base_size = 12) +
    theme(legend.title = element_blank()) +
    theme(
      text = element_text(size = rel(1.1)),
      axis.text = element_text(size = rel(1)),
      axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
      axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
      plot.title = element_text(margin=margin(b=5)),
      plot.caption = element_text(size = rel(1))
    ) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "solid", "solid")) +
    ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
    ggplot2::scale_size_manual(values = c(1, 1, 1))
  
}

# Function plotting actual and filtered values of the log cumulative growth rate
plot_forecast <- function(res, y.eval, n.ahead = 14, plt.start.date=NULL, title="", caption = "") {
  
  # Description: Function plots actual and filtered values of the log cumulative growth rate 
  # (ln(g_t)) within the estimation sample, along with forecasted and realised 
  # values outside the estimation sample.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - y.eval: Out-of-sample realisation of the log growth rate of the cumulated variable 
  #           (actual values for comparison with forecasts).
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
  model <- res$output$model
  est.date.index <- res$index
  
  y <- xts::xts(res$output$model$y %>% as.numeric(), order.by = est.date.index)
  p <- attr(res$output$model, 'p')
  
  y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
  y.pred <-  subset(y.hat.all$y.hat,index(y.hat.all$y.hat) > tail(res$index,1))
  filtered.level <- y.hat.all$level
  
  if (p == 1) {
    EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
    d <- cbind(y, filtered.level, y.pred,
               y.eval[index(y.eval)>tail(index(est.date.index),1),])
    if (!is.null(plt.start.date)) { d <- d[index(d) > plt.start.date] }
    d <- d[index(d) <= tail(index(y.pred),1)]
    names(d) <- c(
      'EstimationSample', 'FilteredLevel', 'Forecast', 'RealisedData'
    )
    
    df_plot <- as.data.frame(d)
    df_plot$Date <- as.Date(rownames(df_plot))
    
    p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
      ggplot2::geom_line(aes(
        y = EstimationSample, color = "Estimation\nSample"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = FilteredLevel, color = "Filtered\nLevel"),
                         lwd = 0.85) +
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData"),
                         lwd = 0.85) +
      ggplot2::scale_color_manual(values = c(1, 2, 3, 'grey')) +
      scale_linetype_manual(
        values = c("solid", "solid", "solid", "dashed")) +
      scale_x_date(labels = scales::date_format("%d %b %y")) +
      labs(x = "Date", y = "Log Growth Rate", caption = caption,
           title = title
      ) +
      theme_economist_white(gray_bg = FALSE) +
      scale_fill_economist() +
      theme(legend.title = element_blank()) +
      theme(
        text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
        axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
        plot.title = element_text(margin=margin(b=5)),
        plot.caption = element_text(size = rel(1)),
      )
  } else if (p == 2) {
    g_1 <- g_2 <- delta <- Forecast <- RealisedData <- NULL
    d <- cbind(y, filtered.level, y.pred[,2],
               y.eval[index(y.eval)>tail(index(est.date.index),1),2])
    d <- d[index(d) <= tail(index(y.pred),1)]
    names(d) <- c('g_1', 'g_2', 'delta', 'Forecast', 'RealisedData')
    
    df_plot <- as.data.frame(d)
    df_plot$Date <- as.Date(rownames(df_plot))
    
    p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
      ggplot2::geom_line(aes(y = g_1, color = "g_1")) +
      ggplot2::geom_line(aes(y = g_2, color = "g_2")) +
      ggplot2::geom_line(aes(y = g_2, color = "delta")) +
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast")) +
      ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData")) +
      ggplot2::scale_color_manual(
        values = c(1, 2, 3, 4, 'grey')) +
      ggplot2::scale_linetype_manual(
        values = c("solid", "solid", "solid", "solid", "dashed")
      ) +
      ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
      labs(x = "Date", y = "Log Growth Rate", caption = caption,
           title = title
      ) +
      theme_economist_white(gray_bg = FALSE) +
      scale_fill_economist() +
      theme(legend.title = element_blank()) +
      theme(
        text = element_text(size = rel(1.)),
        axis.text = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
        axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
        plot.title = element_text(margin=margin(b=5)),
        plot.caption = element_text(size = rel(1))
      )
  } else { stop('NotImplemented Error') }
  
  return(p1)
}

# Function plotting smoothed/filtered growth rates of g_y and g differences
plot_gy_components <- function(res, plt.start.date = NULL, smoothed = FALSE, title = NULL){
  
  # Description: Function plots the smoothed/filtered growth rates of the difference 
  # in the cumulated variable (g_y), the cumulated variable (g), and the slope of 
  # ln(g) (γ). The growth rate g_y is computed as per Harvey and Kattuman (2021):
  # g_y,t = exp(δ_t) + γ_t.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - plt.start.date: Start date for the plot (default: NULL, meaning the start 
  #                   of the estimation sample).
  # - smoothed: Logical value indicating whether to use smoothed estimates of δ 
  #             and γ (default: FALSE, meaning filtered estimates are used).
  # - title: Title for the plot (entered as a text string, default: NULL, meaning 
  #          no title).
  #
  # Returns: A ggplot2 plot visualizing the growth rates and slope of the log  
  # cumulative growth rate.
  
  Date <- Value <- Variable <- NULL
  # Determine plot start date
  if(is.null(plt.start.date)) plt.start.date <- res$index[1]
  
  # Get gy.t, g.t and gamma
  gy.components <- res$get_growth_y(return.components = TRUE, smoothed =
                                      smoothed)
  gy.t <- gy.components[[1]]
  g.t <- gy.components[[2]]
  gamma.t <- gy.components[[3]]
  
  d <- cbind(gy.t,g.t,gamma.t)
  names(d) <- c('gy.t','g.t','gamma.t')
  
  df_plot <- as.data.frame(d)
  df_plot$Date <- as.Date(rownames(df_plot))
  
  df_long <- df_plot %>%
    dplyr::filter(Date >= plt.start.date) %>%
    pivot_longer(cols = c(gy.t, g.t, gamma.t), names_to = "Variable",
                 values_to = "Value")
  
  p1 <- ggplot(df_long, aes(x = Date, y = Value, color = Variable)) +
    geom_line(lwd=0.85) +
    ggplot2::facet_wrap(~ factor(
      Variable, c("gy.t", "g.t", "gamma.t")), ncol = 1, scales = "free_y") +
    labs(title = title, y=ggplot2::element_blank()) +
    scale_color_manual(values = c("#AA2045","darkgrey","black")) +
    scale_x_date(labels = scales::date_format("%d %b %y")) +
    scale_y_continuous(breaks = waiver(), n.breaks = 4) +
    theme_economist_white(gray_bg = FALSE, base_size = 14) +
    theme(text = element_text(size= rel(1), margin=ggplot2::margin(b=5)),
          axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
          legend.position = "none")
  
  return(p1)
}

# Function plotting rate of difference in g_y and g with CIs
plot_gy_ci <- function(res, plt.start.date = NULL, smoothed = FALSE, title = NULL, series.name = NULL, pad.right = NULL){
  
  # Description: Function plots the smoothed/filtered growth rate of the 
  # difference in the cumulated variable (g_y) along with the associated 
  # confidence intervals.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - plt.start.date: Start date for the plot (default: NULL, meaning the start 
  #                   of the estimation sample).
  # - smoothed: Logical value indicating whether to use smoothed estimates of δ 
  #             and γ (default: FALSE, meaning filtered estimates are used).
  # - title: Title for the plot (entered as a text string, default: NULL, meaning 
  #          no title).
  # - series.name: Name of the series for which the growth rate is computed 
  #                (e.g., 'New cases').
  # - pad.right: Numerical value indicating the number of time periods of blank 
  #              space to leave on the right side of the graph (extends the 
  #              horizontal axis).
  #
  # Returns: A ggplot2 plot visualizing the growth rate and confidence intervals.
  
  
  Date <- fit <- upper <- lower <- NULL
  
  # Determine plot start date
  if(is.null(plt.start.date)) plt.start.date <- res$index[1]
  
  # Get confidence intervals to plot
  gy.ci<- res$get_gy_ci(smoothed = smoothed)
  
  y.lab <- if(is.null(series.name)) { c("Growth rate") } else {
    paste("Growth rate of"," ",series.name,sep="")
  }
  
  df_plot <- as.data.frame(gy.ci)
  df_plot$Date <- as.Date(rownames(df_plot))
  
  p1 <- ggplot2::ggplot(df_plot[df_plot$Date>=plt.start.date,], aes(x=Date)) +
    ggplot2::geom_line(aes(y = fit), lwd = 0.85) +
    ggplot2::geom_hline(yintercept=0, linetype="solid",
                        color = "green", linewidth=1)+
    ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045",
                         alpha = 0.3) +
    ggplot2::scale_color_manual(values = c("black")) +
    geom_hline(
      aes(yintercept = 0.0), linetype = "solid", color = "green", lwd = 1.
    ) +
    labs(title=title, x="Date", y=y.lab) +
    theme_economist_white(gray_bg = FALSE, base_size = 14) +
    theme(
      legend.title = element_blank(),
      text = element_text(size = rel(1.)),
      axis.text = element_text(size = rel(1.)),
      axis.title.y = element_text(
        size = rel(1.),margin = ggplot2::margin(r=10)),
      axis.title.x = element_text(
        size = rel(1.),margin = ggplot2::margin(t=10)),
      plot.caption = element_text(size = rel(1))
    ) +
    theme(panel.grid.major.x = ggplot2::element_line(
      color = "gray50", linewidth = 0.5)) +
    scale_linetype_manual(
      values = c("solid")) +
    scale_x_date(labels = scales::date_format("%d %b %y"))
  
  if (!is.null(pad.right)) {
    end.date <- tail(index(gy.ci),1)
    p1 <- p1 +
      ggplot2::scale_x_date(
        limits = c(as.Date(plt.start.date), end.date + pad.right))
  }
  
  return(p1)
}

# Function plotting percentage prediction error
plot_holdout <- function(res, Y, Y.eval, confidence.level = 0.68, date_format = "%Y-%m-%d", series.name = NULL, title= NULL, caption = NULL) {
  
  # Description: Function plots the forecast of new cases (difference in the cumulated variable) 
  # over a holdout sample. It includes actual values, forecasted values (with and 
  # without seasonal components), forecast intervals, and actual outcomes from 
  # the holdout sample. Additionally, it reports the mean absolute percentage 
  # prediction error over the holdout sample.
  #
  # Parameters:
  # - res: Results object estimated using the estimate() method.
  # - Y: Values of the cumulated variable used in the estimation window.
  # - Y.eval: Values of the cumulated variable in the holdout sample 
  #           (actual values for comparison with forecasts).
  # - confidence.level: Width of the prediction interval for ln(g_t) used in 
  #                     forecasts of y_t = ΔY_t (default: 0.68, approximately 
  #                     one standard deviation).
  # - series.name: Name of the variable being forecasted, used for the y-axis 
  #                label (e.g., "Cases" results in "New Cases" on the y-axis).
  # - date_format: Format for date labels on the plot (default: '%Y-%m-%d').
  # - title: Title for the forecast plot (entered as a text string, default: NULL, 
  #          meaning no title).
  # - caption: Caption for the forecast plot (entered as a text string, default: NULL, 
  #            meaning no caption).
  #
  # Returns: A ggplot2 plot visualizing the forecast, actual values, and prediction 
  # intervals over the holdout sample.
  
  
  Date <- Actual <- Forecast <- ForecastTrend <- lower <- upper <- NULL
  
  model <- res$output$model
  est.date.index <- res$index
  
  y.level.est <- Y[est.date.index]
  
  p <- attr(res$output$model, 'p')
  if(p!=1) { stop('NotImplementedError') }
  
  n.ahead <- tail(index(Y.eval),1)-tail(index(y.level.est),1)
  
  y.eval.diff <- diff(Y.eval) %>% na.omit
  
  est.date.index <- res$index %>% as.Date()
  estimation.date.end <- tail(est.date.index, 1)
  
  y.hat.diff.final.ci <- res$predict_level(
    y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
    return.diff = TRUE
  )
  y.hat.diff.final <- res$predict_level(
    y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
    sea.on = TRUE,
    return.diff = TRUE
  )
  
  d <- cbind(
    y.eval.diff[index(y.eval.diff)>estimation.date.end,],
    y.hat.diff.final[, 1],
    y.hat.diff.final.ci[, 1]
  )
  names(d) <- c('Actual', 'Forecast', 'ForecastTrend')
  
  df_plot <- as.data.frame(d)
  df_plot$Date <- as.Date(rownames(df_plot), format = date_format)
  
  d.eval <- na.omit(d)
  mape.trend <- 100*(abs(d.eval$Actual - d.eval$`ForecastTrend`)/
                       d.eval$Actual) %>% mean %>% round(2)
  mape.sea <- 100*(abs(d.eval$Actual - d.eval$Forecast)/d.eval$Actual) %>%
    mean %>% round(2)
  
  ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                            (as.Date(index(y.hat.diff.final.ci),
                                     format = date_format))))
  colnames(ci) <- c('lower', 'upper', 'date')
  ci[, 'date'] <- as.Date(ci[, 'date'], format = date_format,
                          origin = "1970-01-01")
  
  p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
    ggplot2::geom_line(aes(y = Actual, color = "Actual"),lwd = 0.85) +
    ggplot2::geom_line(aes(y = Forecast, color = "Forecast"),lwd = 0.85) +
    ggplot2::geom_line(
      aes(y = ForecastTrend, color = "Forecast\nTrend"),lwd = 0.85) +
    ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, aes(x = date, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045",
                         alpha = 0.1) +
    labs(x = "Date", y = paste("New",series.name), title = title,
         subtitle = paste("MAPE: ",mape.sea,"%. Trend MAPE: ",
                          mape.trend,"%.",sep="")) +
    theme_economist_white(gray_bg = FALSE, base_size = 14) +
    theme(legend.title = element_blank()) +
    theme(
      text = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      axis.title.y = element_text(size = rel(1), margin = margin(r=10)),
      axis.title.x = element_text(size = rel(1), margin = margin(t=10)),
      plot.title = element_text(margin=margin(b=5)),
      plot.subtitle = element_text(
        size = rel(1), hjust=0,  margin = margin(t=3))
    ) +
    scale_linetype_manual(
      values = c("solid", "solid", "solid")) +
    scale_x_date(labels = scales::date_format("%d %b %y")) +
    scale_size_manual(values = c(1, 1.5, 1))
  
  
  return(p1)
}