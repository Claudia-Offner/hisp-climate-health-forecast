### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city meta-analysis                                         ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork); 
library(colorspace); library(png); library(grid); library(reshape2)


setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

#### FUNCTIONS ####


#### MAIN CODE ####

load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# Run model ----

# Total patients
primary <- 'total_patients'
output <- run_model(df, locations, primary)
output$meta_comp$sum_aic

# Main results
main_res <- format_results(output, primary)
(total_er_overall <- main_res$p_er_overall + plot_annotation(title=replace_exact_matches(primary, nm)))
(total_lr <- main_res$p_lr + plot_annotation(title=paste0('Lag-response: ', replace_exact_matches(primary, nm))))
(total_er <- main_res$p_er + plot_annotation(title=paste0('Exposure-response: ', replace_exact_matches(primary, nm))))

# Main percentile look up
datalist <- lapply(locations, function(region) df[df$location==region,])
cen_wb <- sapply(datalist,function(x) {
  quantile(x$wb,c(0:100)/100,na.rm=T)
})
colnames(cen_wb) <- locations
cen_wb

# Winter/pre-monsoon results (months=5)
df_win <- df %>%
  filter(season %in% c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)'))
win <- run_model(df_win, locations, primary)
win$meta_comp$sum_aic
win_res <- format_results(win, primary)
(win_er_overall <- win_res$p_er_overall + plot_annotation(title='Winter'))
(win_lr <- win_res$p_lr + plot_annotation(title='Lag-response: Winter'))
(win_er <- win_res$p_er + plot_annotation(title='Exposure-response: Winter'))
# Winter percentile look up
datalist <- lapply(locations, function(region) df_win[df_win$location==region,])
cen_wb <- sapply(datalist,function(x) {
  quantile(x$wb,c(0:100)/100,na.rm=T)
})
colnames(cen_wb) <- locations
cen_wb

# Monsoon/Post-monsoon, summer filter (months=7)
df_sum <- df %>%
  filter(season %in% c('Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)'))
sum <- run_model(df_sum, locations, primary)
sum$meta_comp$sum_aic
sum_res <- format_results(sum, primary)
(sum_er_overall <- sum_res$p_er_overall + plot_annotation(title='Summer'))
(sum_lr <- sum_res$p_lr + plot_annotation(title='Lag-response: Summer'))
(sum_er <- sum_res$p_er + plot_annotation(title='Exposure-response: Summer'))
# Summer percentile look up
datalist <- lapply(locations, function(region) df_sum[df_sum$location==region,])
cen_wb <- sapply(datalist,function(x) {
  quantile(x$wb,c(0:100)/100,na.rm=T)
})
colnames(cen_wb) <- locations
cen_wb






# Example contour plots
city_con_res(meta=output$meta_wb, xlab='Temperature (%)', title='Mean temp')
city_con_res(meta=output$meta_tp, xlab='Rainfall (%)', title='5-day cummulative Rainfall')

# Main exposure-response, by city
(city_wb_er <- main_res$p_er_wb_city + plot_annotation(title=replace_exact_matches('wb', nm)))
(city_tp_er <- main_res$p_er_tp_city + plot_annotation(title=replace_exact_matches('tp', nm)))
(city_rh_er <- main_res$p_er_rh_city + plot_annotation(title=replace_exact_matches('rh', nm)))
# (city_fg10_er <- main_res$p_er_fg10_city + plot_annotation(title=replace_exact_matches('fg10', nm)))
# (city_pm25_er <- main_res$p_er_pm25_city + plot_annotation(title=replace_exact_matches('pm25', nm)))

# Main lag-response, by city
# (city_wb_lr <- main_res$p_lr_wb_city + plot_annotation(title=replace_exact_matches('wb', nm)))
# (city_tp_lr <- main_res$p_lr_tp_city + plot_annotation(title=replace_exact_matches('tp', nm)))
# (city_rh_lr <- main_res$p_lr_rh_city + plot_annotation(title=replace_exact_matches('rh', nm)))
# (city_fg10_lr <- main_res$p_lr_fg10_city + plot_annotation(title=replace_exact_matches('fg10', nm)))
# (city_pm25_lr <- main_res$p_lr_pm25_city + plot_annotation(title=replace_exact_matches('pm25', nm)))

# City case study: Nahan & Recking Peo
nahan_res <- format_results(output, primary, city='Nahan')
(nahan_lr <- wrap_plots(nahan_res$p_lr_wb_city, nahan_res$p_lr_tp_city, nahan_res$p_lr_rh_city,
                        nrow=3)+plot_annotation(title='Lag-response: Nahan', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
(nahan_er <- wrap_plots(nahan_res$p_er_wb_city, nahan_res$p_er_tp_city, nahan_res$p_er_rh_city,
                        nrow=3)+plot_annotation(title='Exposure-response: Nahan', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
reck_res <- format_results(output, primary, city='Reckong Peo')
(reck_lr <- wrap_plots(reck_res$p_lr_wb_city, reck_res$p_lr_tp_city, reck_res$p_lr_rh_city,
                        nrow=3)+plot_annotation(title='Lag-response: Reckong Peo', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
(reck_er <- wrap_plots(reck_res$p_er_wb_city, reck_res$p_er_tp_city, reck_res$p_er_rh_city,
                        nrow=3)+plot_annotation(title='Exposure-response: Reckong Peo', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))

# Check map of elevation
path <- paste0(path_git, '1_Data/HISP_Himachal_Pradesh/- Reporting/- figures/desc_study_area.png')
img <- readPNG(path)
grid.newpage()
grid.raster(img)

# Check population
sort(sapply(locations, function(loc) {mean(unique(df[['popu']][df$location==loc]))}))

# META REGRESSION: Elevation and population density
meta_reg <- function(meta, xlab='Exposure (%)', predictor='elev', perc=c(0.1, 0.9), 
                     col='darkred'){
  
  # Get unique predictors
  pred <- sapply(locations, function(loc) {
    if(predictor=='popu'){
      mean(unique(df[[predictor]][df$location==loc]))
    } else{
      unique(df[[predictor]][df$location==loc])
    }
  })
  
  # Remove duplicates before predictions
  rows_with_dup <- apply(meta$cen, 1, function(row) any(duplicated(row)))
  cen <- meta$cen[!rows_with_dup, ]
  
  # Get pooled regression estimates
  (mv_reg <- mvmeta(meta$ymat~pred,meta$Slist,method="reml"))
  pred_perc <- round(quantile(pred,perc,1))
  predict <- predict(mv_reg,data.frame(pred=pred_perc),vcov=T)
  cp_p1 <- crosspred(meta$b_mean,coef=predict[[1]]$fit, vcov=predict[[1]]$vcov,model.link="log", 
                     at=rowMeans(cen), cen=rowMeans(cen)["50%"])
  cp_p2 <- crosspred(meta$b_mean,coef=predict[[2]]$fit, vcov=predict[[2]]$vcov,model.link="log", 
                     at=rowMeans(cen), cen=rowMeans(cen)["50%"])
  
  # Format tables for plotting
  dfp1 <- data.frame(exp=as.numeric(sub("%", "", names(rowMeans(cen)))),
                     perc=paste0(perc[1]*100, "%"), RR =cp_p1$allRRfit, 
                     low=cp_p1$allRRlow, high=cp_p1$allRRhigh)
  dfp2 <- data.frame(exp=as.numeric(sub("%", "", names(rowMeans(cen)))),
                     perc=paste0(perc[2]*100, "%"), RR =cp_p2$allRRfit, 
                     low=cp_p2$allRRlow, high=cp_p2$allRRhigh)
  df_plot <- rbind(dfp1, dfp2)
  df_plot$cross_one <- ifelse(df_plot$low <= 1 & df_plot$high >= 1, "cross", "no_cross")
  
  # Plot
  ggplot(df_plot, aes(x = exp, y = RR, color = perc)) + #linetype = perc
    geom_hline(yintercept = 1, color = 'gray20') +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = perc), alpha = 0.2, color = NA) +
    labs(
      x = xlab,
      y = "RR",
      title = "",
      color = "",
      fill  = "",
      # linetype = ""
    ) +
    scale_color_manual(values = c('gray20', col)) +
    scale_fill_manual(values = c('gray20', col)) +
    # scale_linetype_manual(values = c("dotted", "solid")) +  # assign specific line types
    theme_classic() +
    theme(
      axis.line = element_line(color = 'gray20'),
      axis.text = element_text(color = 'gray20'),
      axis.ticks = element_line(color = 'gray20'),
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 12),
      legend.position = 'bottom'
    )
  
}


meta_list <- list('meta_wb'=list(xlab='Mean temp. (%)', col='#D74B3C'),
                  'meta_tp'=list(xlab='Rainfall (%)', col='#1D98EB'),
                  'meta_rh'=list(xlab='Relative humidity (%)', col='#A6B727'),
                  'meta_fg10'=list(xlab='Wind speed (%)', col='#8968CD'),
                  'meta_pm25'=list(xlab='Air pollution (%)', col='#E569BF'))
elev_plots <- list()
popu_plots <- list()
for(v in names(meta_list)){
  elev_reg <- meta_reg(output[[v]], predictor='elev', perc=c(0.05, 0.95),
                       xlab=meta_list[[v]]$xlab, meta_list[[v]]$col)
  elev_plots[[v]] <- elev_reg
  popu_reg <- meta_reg(output[[v]], predictor='popu', perc=c(0.05, 0.95),
                       xlab=meta_list[[v]]$xlab, meta_list[[v]]$col)
  popu_plots[[v]] <- popu_reg
}

(elev_res <- wrap_plots(elev_plots, ncol=5)+plot_annotation(title='Exposure-resonse by elevation percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
(popu_res <- wrap_plots(popu_plots, ncol=5)+plot_annotation(title='Exposure-resonse by population density percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))



# # Get total tables and id relevant lags/exposure
# x <- meta_checks(output$meta_wb, 1.001)
# meta_checks(output$meta_tp, 1.001)
# meta_checks(output$meta_pm25, 1.000)

# Multiple outcomes
outcomes <- h_cols[1:8]
er14_wb <- list()
er0_tp <- list()
er14_tp <- list()
er10_rh <- list()
# er14_pm25 <- list()

for (i in outcomes){
  output <- run_model(df, locations, i)
  res <- format_results(output, i)
  er14_wb[[i]] <- res$lags$`14`$p_er_wb$p_blup + labs(title=replace_exact_matches(i, nm))
  er0_tp[[i]] <- res$lags$`0`$p_er_tp$p_blup + labs(title=replace_exact_matches(i, nm))
  er14_tp[[i]] <- res$lags$`14`$p_er_tp$p_blup + labs(title=replace_exact_matches(i, nm))
  er10_rh[[i]] <- res$lags$`10`$p_er_rh$p_blup + labs(title=replace_exact_matches(i, nm))
  # er14_pm25[[i]] <- res$lags$`14`$p_er_pm25$p_blup + labs(title=replace_exact_matches(i, nm))
}

(er14_wb <- (wrap_plots(er14_wb, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er0_tp <- (wrap_plots(er0_tp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=0)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er14_tp <- (wrap_plots(er14_tp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er10_rh <- (wrap_plots(er10_rh, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
# (er14_pm25 <- (wrap_plots(er14_pm25, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))

all_outcomes <- list(er5_wb=er5_wb, er10_wb=er10_wb, 
                     er0_tp=er0_tp, er14_rh=er14_rh, er14_pm25=er14_pm25)

# Specific location

#### EXPORT #####
width <- 10
height <- 6

# Save cumulative exposure-response (total_patient)
ggsave(filename=paste0(path_ch2, '/IV. Figures/MF_total_er_overall.png'), plot=total_er_overall, width=10, height=5, dpi=300)
ggsave(filename=paste0(path_report, '- figures/MF_total_er_overall.png'), plot=total_er_overall, width=10, height=5, dpi=300)

# Save lag-response by percentile (total_patient)
ggsave(filename=paste0(path_ch2, '/IV. Figures/MF_total_lr.png'), plot=total_lr, width=width, height=height, dpi=300)
ggsave(filename=paste0(path_report, '- figures/MF_total_lr.png'), plot=total_lr, width=width, height=height, dpi=300)

# Save exposure-response by lag (total_patient)
ggsave(filename=paste0(path_ch2, '/IV. Figures/MF_total_er.png'), plot=total_er, width=width, height=height, dpi=300)
ggsave(filename=paste0(path_report, '- figures/MF_total_er.png'), plot=total_er, width=width, height=height, dpi=300)

# Save cumulative exposure-response by city (total_patient)
for (n in names(all_cities)) {
  ggsave(filename=paste0(path_ch2, '/IV. Figures/SF_', n, '.png'), plot=all_cities[[n]], width=10, height=8, dpi=300)
  ggsave(filename=paste0(path_report, '- figures/SF_', n, '.png'), plot=all_cities[[n]], width=10, height=6, dpi=300)
}

# Population plots
names(all_outcomes) <- paste0("disag_", names(all_outcomes))
for (n in names(all_outcomes)) {
  ggsave(filename=paste0(path_ch2, '/IV. Figures/MF_', n, '.png'), plot=all_outcomes[[n]], width=10, height=8, dpi=300)
  ggsave(filename=paste0(path_report, '- figures/MF_', n, '.png'), plot=all_outcomes[[n]], width=10, height=6, dpi=300)
}








