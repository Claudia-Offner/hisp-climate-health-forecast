### ------------------------------------------------------------------------ ### 
###  DNLM - Multi-city meta-analysis                                         ###
### ------------------------------------------------------------------------ ### 

# LOAD THE PACKAGE
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork); 
library(colorspace); library(png); library(grid); library(reshape2); library(cowplot)


setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

#### FUNCTIONS ####

# Function to create city specific legend from dummy set
get_city_legend <- function(plot, pred='ALL'){
  
  legend_df <- data.frame(x = 1, y = 1, 
                          predictor = factor(c("wb","tp","rh", 'pm25', 'fg10', "Pooled","BLUP"), 
                                             levels = c("wb","tp","rh",'pm25', 'fg10',"Pooled","BLUP")), 
                          colour = c("#D74B3C","#1D98EB","#A6B727", '#E569BF', '#8968CD', "black","black"), 
                          linetype = c("solid","solid","solid","solid","solid","solid","dashed"), 
                          size = c(1.2,1.2,1.2,1.2,1.2,0.7,0.7))
  legend_df <- legend_df %>% filter(!predictor %in% c("pm25", "fg10"))
  
  if(pred != 'ALL'){
    legend_df <- subset(legend_df, predictor %in% c(pred, "Pooled", "BLUP"))
  }
  
  vec <- legend_df$predictor

  p <- ggplot(legend_df, aes(x, y, colour = predictor, linetype = predictor, size = predictor)) +
    geom_line() +
    scale_colour_manual(
      values = setNames(legend_df$colour, legend_df$predictor),
      breaks = c(vec[-((length(vec)-1):length(vec))]),
      labels = replace_exact_matches(as.data.frame(c(vec[-((length(vec)-1):length(vec))])), nm)[,1],
      guide = guide_legend(title = 'Raw', nrow = length(c(vec[-((length(vec)-1):length(vec))])), 
                           order = 1, 
                           override.aes = list(linewidth = rep(1, length(c(vec[-((length(vec)-1):length(vec))])))))) +
    scale_linetype_manual(
      values = setNames(legend_df$linetype, legend_df$predictor),
      breaks = c(vec[((length(vec)-1):length(vec))]),
      guide = guide_legend(title = 'Meta-analysis', nrow = 2, order = 2, 
                           override.aes = list(colour = "black", 
                                               linewidth = c(0.7, 0.7)))) +
    scale_size_manual(
      values = setNames(legend_df$size, legend_df$predictor),
      breaks = c(vec[((length(vec)-1):length(vec))]),
      guide = "none") +  theme_void() +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
  
  city_legend <- get_legend(p)
  
  res <- plot_grid(plot, city_legend, ncol = 2, rel_widths  = c(10, 2))
  
  return(res)
  
}

#### MAIN CODE ####

load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))

# Set location info
locations <- as.character(unique(df$location))
m <- length(locations)

# Run model 1 - Overall ----

# Total patients
outcome <- 'total_patients'
output <- run_model(df, locations, outcome)
output$meta_comp$sum_aic

# Main results
main_res <- format_results(output, outcome)
# Contours
# city_con_res(meta=output$meta_wb, xlab='Temperature (%)', title='Wet-bulb temperature')
# city_con_res(meta=output$meta_tp, xlab='Rainfall (%)', title='5-day cummulative precipitation')
# city_con_res(meta=output$meta_rh, xlab='Relative humidity (%)', title='Relative humidity')
# Overall
(total_er_overall <- main_res$p_er_overall + plot_annotation(title=replace_exact_matches(outcome, nm)))
(total_lr <- main_res$p_lr + plot_annotation(title=paste0('Lag-response: ', replace_exact_matches(outcome, nm))))
(total_er <- main_res$p_er + plot_annotation(title=paste0('Exposure-response: ', replace_exact_matches(outcome, nm))))

# Main exposure-response, by city
(city_wb_er <- get_city_legend(pred='wb', plot=main_res$p_er_wb_city + plot_annotation(title='Exposure-response, by city')))
(city_tp_er <- get_city_legend(pred='tp', plot=main_res$p_er_tp_city + plot_annotation(title='Exposure-response, by city')))
(city_rh_er <- get_city_legend(pred='rh', plot=main_res$p_er_rh_city + plot_annotation(title='Exposure-response, by city')))
# (city_pm25_er <- get_city_legend(pred='pm25', plot=main_res$p_er_pm25_city + plot_annotation(title='Exposure-response, by city')))
# (city_fg10_er <- get_city_legend(pred='fg10', plot=main_res$p_er_fg10_city + plot_annotation(title='Exposure-response, by city')))

# Main lag-response, by city
# (city_wb_lr <- main_res$p_lr_wb_city + plot_annotation(title=replace_exact_matches('wb', nm)))
# (city_tp_lr <- main_res$p_lr_tp_city + plot_annotation(title=replace_exact_matches('tp', nm)))
# (city_rh_lr <- main_res$p_lr_rh_city + plot_annotation(title=replace_exact_matches('rh', nm)))
# (city_fg10_lr <- main_res$p_lr_fg10_city + plot_annotation(title=replace_exact_matches('fg10', nm)))
# (city_pm25_lr <- main_res$p_lr_pm25_city + plot_annotation(title=replace_exact_matches('pm25', nm)))

# City case study: Una & Recking Peo
una_res <- format_results(output, outcome, city='Una')
(una_lr <- get_city_legend(wrap_plots(una_res$p_lr_wb_city, una_res$p_lr_tp_city, una_res$p_lr_rh_city,
                                       nrow=3)+plot_annotation(title='Lag-response: Una', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))))
(una_er <- get_city_legend(wrap_plots(una_res$p_er_wb_city, una_res$p_er_tp_city, una_res$p_er_rh_city,
                                        nrow=3)+plot_annotation(title='Exposure-response: Una', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))))

solan_res <- format_results(output, outcome, city='Solan')
(solan_lr <- get_city_legend(wrap_plots(solan_res$p_lr_wb_city, solan_res$p_lr_tp_city, solan_res$p_lr_rh_city,
                                      nrow=3)+plot_annotation(title='Lag-response: Solan', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))))
(solan_er <- get_city_legend(wrap_plots(solan_res$p_er_wb_city, solan_res$p_er_tp_city, solan_res$p_er_rh_city,
                                      nrow=3)+plot_annotation(title='Exposure-response: Solan', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))))


# Main percentile look up
cen_tables <- list()
datalist <- lapply(locations, function(region) df[df$location==region,])
variables <- c('wb', 'tp', 'rh', 'pm25', 'fg10')
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_tables[[v]] <- round(cen, 4)
}
names(cen_tables) <- c('Wet bulb temp', 'Precipitation', 'Relative humidity', 'Air pollution', 'Wind speed')

#### EXPORT #####

# FIGURES
exports <- list('MF_total_er_overall'=list(plot=total_er_overall, width=10, height=5),
                'MF_total_lr'=list(plot=total_lr, width=10, height=6),
                'MF_total_er'=list(plot=total_er, width=10, height=6),
                'MF_city_wb_er'=list(plot=city_wb_er, width=10, height=5),
                'MF_city_tp_er'=list(plot=city_tp_er, width=10, height=5),
                'MF_city_rh_er'=list(plot=city_rh_er, width=10, height=5),
                'SF_una_lr'=list(plot=una_lr, width=10, height=6),
                'SF_una_er'=list(plot=una_er, width=10, height=6),
                'SF_solan_lr'=list(plot=solan_lr, width=10, height=6),
                'SF_solan_er'=list(plot=solan_er, width=10, height=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_ch2, 'IV. Figures/', e, '.png'), plot=exports[[e]]$plot, width=exports[[e]]$width, height=exports[[e]]$height, dpi=300)
  ggsave(filename=paste0(path_report, '- figures/', e, '.png'), plot=exports[[e]]$plot, width=exports[[e]]$width, height=exports[[e]]$height, dpi=300)
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_ch2,"III. Tables/ST_main_perc.xlsx"), overwrite = TRUE)
saveWorkbook(wb, file=paste0(path_report,"- tables/ST_main_perc.xlsx"), overwrite = TRUE)

# Run model 2 - Seasonal differences ####

# Winter/pre-monsoon results (months=5)
df_win <- df %>%
  filter(season %in% c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)'))
win <- run_model(df_win, locations, outcome)
win$meta_comp$sum_aic
win_res <- format_results(win, outcome)
(win_er_overall <- win_res$p_er_overall + plot_annotation(title='Winter'))
(win_lr <- win_res$p_lr + plot_annotation(title='Lag-response: Winter'))
(win_er <- win_res$p_er + plot_annotation(title='Exposure-response: Winter'))
# Winter percentile look up
cen_win_tables <- list()
datalist <- lapply(locations, function(region) df_win[df_win$location==region,])
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_win_tables[[v]] <- round(cen, 4)
}
names(cen_win_tables) <- c('Wet bulb temp', 'Precipitation', 'Relative humidity', 'Air pollution', 'Wind speed')

# Monsoon/Post-monsoon, summer filter (months=7)
df_sum <- df %>%
  filter(season %in% c('Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)'))
sum <- run_model(df_sum, locations, outcome)
sum$meta_comp$sum_aic
sum_res <- format_results(sum, outcome)
(sum_er_overall <- sum_res$p_er_overall + plot_annotation(title='Summer'))
(sum_lr <- sum_res$p_lr + plot_annotation(title='Lag-response: Summer'))
(sum_er <- sum_res$p_er + plot_annotation(title='Exposure-response: Summer'))
# Summer percentile look up
cen_sum_tables <- list()
datalist <- lapply(locations, function(region) df_sum[df_sum$location==region,])
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_sum_tables[[v]] <- round(cen, 4)
}
names(cen_sum_tables) <- c('Wet bulb temp', 'Precipitation', 'Relative humidity', 'Air pollution', 'Wind speed')
  
#### EXPORT #####
exports <- list('MF_win_er_overall'=list(plot=win_er_overall, width=10, height=5),
                'MF_win_lr'=list(plot=win_lr, width=10, height=6),
                'MF_win_er'=list(plot=win_er, width=10, height=6),
                'MF_sum_er_overall'=list(plot=sum_er_overall, width=10, height=5),
                'MF_sum_lr'=list(plot=sum_lr, width=10, height=6),
                'MF_sum_er'=list(plot=sum_er, width=10, height=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_ch2, 'IV. Figures/', e, '.png'), plot=exports[[e]]$plot, width=exports[[e]]$width, height=exports[[e]]$height, dpi=300)
  ggsave(filename=paste0(path_report, '- figures/', e, '.png'), plot=exports[[e]]$plot, width=exports[[e]]$width, height=exports[[e]]$height, dpi=300)
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_win_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_win_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_ch2,"III. Tables/ST_winter_perc.xlsx"), overwrite = TRUE)
saveWorkbook(wb, file=paste0(path_report,"- tables/ST_winter_perc.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
for (name in names(cen_sum_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_sum_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_ch2,"III. Tables/ST_summer_perc.xlsx"), overwrite = TRUE)
saveWorkbook(wb, file=paste0(path_report,"- tables/ST_summer_perc.xlsx"), overwrite = TRUE)


# Run model 3 - Population differences ####

# # Get total tables and id relevant lags/exposure
# x <- meta_checks(output$meta_wb, 1.001)
# meta_checks(output$meta_tp, 1.001)
# meta_checks(output$meta_pm25, 1.000)

# Multiple outcomes
outcomes <- h_cols[1:8]
er14_wb <- list()
er5_tp <- list()
er10_rh <- list()

for (i in outcomes){
  output <- run_model(df, locations, i)
  res <- format_results(output, i)
  er14_wb[[i]] <- res$lags$`14`$p_er_wb$p_blup + labs(title=replace_exact_matches(i, nm))
  er5_tp[[i]] <- res$lags$`5`$p_er_tp$p_blup + labs(title=replace_exact_matches(i, nm))
  er10_rh[[i]] <- res$lags$`10`$p_er_rh$p_blup + labs(title=replace_exact_matches(i, nm))
}

(er14_wb <- (wrap_plots(er14_wb, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=14)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er5_tp <- (wrap_plots(er5_tp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=5)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er10_rh <- (wrap_plots(er10_rh, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Expsoure-response (Lag=10)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))

all_outcomes <- list(er14_wb=er14_wb, er5_tp=er5_tp, er10_rh=er10_rh)

#### EXPORT ####

# Population plots
names(all_outcomes) <- paste0("disag_", names(all_outcomes))
for (n in names(all_outcomes)) {
  ggsave(filename=paste0(path_ch2, '/IV. Figures/MF_', n, '.png'), plot=all_outcomes[[n]], width=10, height=8, dpi=300)
  ggsave(filename=paste0(path_report, '- figures/MF_', n, '.png'), plot=all_outcomes[[n]], width=10, height=6, dpi=300)
}

#################################

# # META REGRESSION: Elevation and population density
# meta_reg <- function(meta, xlab='Exposure (%)', predictor='elev', perc=c(0.1, 0.9), 
#                      col='darkred'){
#   
#   # Get unique predictors
#   pred <- sapply(locations, function(loc) {
#     if(predictor=='popu'){
#       mean(unique(df[[predictor]][df$location==loc]))
#     } else{
#       unique(df[[predictor]][df$location==loc])
#     }
#   })
#   
#   # Remove duplicates before predictions
#   rows_with_dup <- apply(meta$cen, 1, function(row) any(duplicated(row)))
#   cen <- meta$cen[!rows_with_dup, ]
#   
#   # Get pooled regression estimates
#   (mv_reg <- mvmeta(meta$ymat~pred,meta$Slist,method="reml"))
#   pred_perc <- round(quantile(pred,perc,1))
#   predict <- predict(mv_reg,data.frame(pred=pred_perc),vcov=T)
#   cp_p1 <- crosspred(meta$b_mean,coef=predict[[1]]$fit, vcov=predict[[1]]$vcov,model.link="log", 
#                      at=rowMeans(cen), cen=rowMeans(cen)["50%"])
#   cp_p2 <- crosspred(meta$b_mean,coef=predict[[2]]$fit, vcov=predict[[2]]$vcov,model.link="log", 
#                      at=rowMeans(cen), cen=rowMeans(cen)["50%"])
#   
#   # Format tables for plotting
#   dfp1 <- data.frame(exp=as.numeric(sub("%", "", names(rowMeans(cen)))),
#                      perc=paste0(perc[1]*100, "%"), RR =cp_p1$allRRfit, 
#                      low=cp_p1$allRRlow, high=cp_p1$allRRhigh)
#   dfp2 <- data.frame(exp=as.numeric(sub("%", "", names(rowMeans(cen)))),
#                      perc=paste0(perc[2]*100, "%"), RR =cp_p2$allRRfit, 
#                      low=cp_p2$allRRlow, high=cp_p2$allRRhigh)
#   df_plot <- rbind(dfp1, dfp2)
#   df_plot$cross_one <- ifelse(df_plot$low <= 1 & df_plot$high >= 1, "cross", "no_cross")
#   
#   # Plot
#   ggplot(df_plot, aes(x = exp, y = RR, color = perc)) + #linetype = perc
#     geom_hline(yintercept = 1, color = 'gray20') +
#     geom_line(size = 1.2) +
#     geom_ribbon(aes(ymin = low, ymax = high, fill = perc), alpha = 0.2, color = NA) +
#     labs(
#       x = xlab,
#       y = "RR",
#       title = "",
#       color = "",
#       fill  = "",
#       # linetype = ""
#     ) +
#     scale_color_manual(values = c('gray20', col)) +
#     scale_fill_manual(values = c('gray20', col)) +
#     # scale_linetype_manual(values = c("dotted", "solid")) +  # assign specific line types
#     theme_classic() +
#     theme(
#       axis.line = element_line(color = 'gray20'),
#       axis.text = element_text(color = 'gray20'),
#       axis.ticks = element_line(color = 'gray20'),
#       plot.title = element_text(face = 'bold', hjust = 0.5, size = 12),
#       legend.position = 'bottom'
#     )
#   
# }
# 
# 
# meta_list <- list('meta_wb'=list(xlab='Mean temp. (%)', col='#D74B3C'),
#                   'meta_tp'=list(xlab='Rainfall (%)', col='#1D98EB'),
#                   'meta_rh'=list(xlab='Relative humidity (%)', col='#A6B727'))
#                   # 'meta_fg10'=list(xlab='Wind speed (%)', col='#8968CD'),
#                   # 'meta_pm25'=list(xlab='Air pollution (%)', col='#E569BF'))
# elev_plots <- list()
# popu_plots <- list()
# for(v in names(meta_list)){
#   elev_reg <- meta_reg(output[[v]], predictor='elev', perc=c(0.51, 0.49),
#                        xlab=meta_list[[v]]$xlab, meta_list[[v]]$col)
#   elev_plots[[v]] <- elev_reg
#   popu_reg <- meta_reg(output[[v]], predictor='popu', perc=c(0.51, 0.49),
#                        xlab=meta_list[[v]]$xlab, meta_list[[v]]$col)
#   popu_plots[[v]] <- popu_reg
# }
# 
# (elev_res <- wrap_plots(elev_plots, ncol=3)+plot_annotation(title='Exposure-resonse by elevation percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# (popu_res <- wrap_plots(popu_plots, ncol=3)+plot_annotation(title='Exposure-resonse by population density percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# 







