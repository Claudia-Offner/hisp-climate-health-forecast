### ------------------------------------------------------------------------ ### 
###  DNLM: Multiple city meta- and uni-analysis
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')
source('2_DLNM/DLNM - Dependencies.R')

# Load packages
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork); library(openxlsx);
library(colorspace); library(png); library(grid); library(reshape2); library(cowplot)

#### FUNCTIONS ####

# Function to get empty plot with title, for nested wrap_plot work around
get_title <- function(title){
  
  ggplot() +
    labs(title = title) +param
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
}

# Function to create city specific legend from dummy set
get_city_legend <- function(plot, pred='ALL'){
  
  legend_df <- data.frame(x = 1, y = 1, 
                          predictor = factor(c("mx2t","tp","fg10", 'pm25', 'rh', "Pooled","BLUP"), 
                                             levels = c("mx2t","tp","fg10",'pm25', 'rh',"Pooled","BLUP")), 
                          colour = c("#D74B3C","#1D98EB","#8968CD", '#E569BF', '#A6B727', "black","black"), 
                          linetype = c("solid","solid","solid","solid","solid","solid","dashed"), 
                          size = c(1.2,1.2,1.2,1.2,1.2,0.7,0.7))
  legend_df <- legend_df %>% filter(!predictor %in% c('fg10', "pm25", "rh"))
  
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

# Function to easily assign parameters according to uni- or multi-site models
run_bestfit <- function(df, locations, outcome, params, pooled, city){
  
  # Get params
  longterm_df <- params$longterm
  autocorr_lag <- params$autocorr
  arglag_tmp <- params$arglag$mx2t
  arglag_tp <- params$arglag$tp
  arglag_fg10 <- params$arglag$fg10
  argvar_tmp <- params$argvar$mx2t
  argvar_tp <- params$argvar$tp
  argvar_fg10 <- params$argvar$fg10
  fam_type <- 'binom' #params$shape
  
  # Run model
  output <- run_model(df, locations, outcome, fam_type=fam_type, pooled=pooled,
                      longterm_df=longterm_df, autocorr_lag=autocorr_lag, 
                      arglag_tmp=arglag_tmp, arglag_tp=arglag_tp, arglag_fg10=arglag_fg10,
                      argvar_tmp=argvar_tmp, argvar_tp=argvar_tp, argvar_fg10=argvar_fg10)
  results <- format_results(output, outcome, city=city)
  
  print(paste0('Model AIC: ', output$meta_comp$sum_aic))
  
  return(list(output=output, results=results))
  
}


#### MAIN CODE ####

# Load data and parameters
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
load(paste0(path_project, '/ind_sens_params.RData'))

# Set location info by elev
locations <- as.character(df %>% distinct(location, elev) %>% arrange(elev) %>% pull(location))
m <- length(locations)

# Store results
results <- list()

# Set outcome outcome and exposures
outcome <- 'total_patients'
variables <- c('mx2t', 'tp')

# Analysis 1 - Overall (pooled params) ----

# Main results
res <- run_bestfit(df, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$overall <- res
output   <- res$output
main_res <- res$results

# Contours
# city_con_res(meta=output$meta_wb, xlab='Temperature (%)', title='Wet-bulb temperature')
# city_con_res(meta=output$meta_tp, xlab='Rainfall (%)', title='5-day cummulative precipitation')
# city_con_res(meta=output$meta_rh, xlab='Relative humidity (%)', title='Relative humidity')

# Overall
(total_er_overall <- main_res$p_er_overall + plot_annotation(title=replace_exact_matches(outcome, nm)))
(total_lr <- main_res$p_lr + plot_annotation(title='Lag-response'))
(total_er <- main_res$p_er + plot_annotation(title='Exposure-response'))

# Main exposure-response, by city
(city_tmp_er <- get_city_legend(pred='mx2t', plot=main_res$p_er_tmp_city + plot_annotation(title='Exposure-response, by city: Temperature')))
(city_tp_er <- get_city_legend(pred='tp', plot=main_res$p_er_tp_city + plot_annotation(title='Exposure-response, by city: Precipitation')))
# (city_fg10_er <- get_city_legend(pred='fg10', plot=main_res$p_er_fg10_city + plot_annotation(title='Exposure-response, by city: Wind speed')))
# (city_rh_er <- get_city_legend(pred='rh', plot=main_res$p_er_rh_city + plot_annotation(title='Exposure-response, by city')))
# (city_pm25_er <- get_city_legend(pred='pm25', plot=main_res$p_er_pm25_city + plot_annotation(title='Exposure-response, by city')))

# Main lag-response, by city
# (city_tmp_lr <- main_res$p_lr_tmp_city + plot_annotation(title=replace_exact_matches('wb', nm)))
# (city_tp_lr <- main_res$p_lr_tp_city + plot_annotation(title=replace_exact_matches('tp', nm)))
# (city_fg10_lr <- main_res$p_lr_fg10_city + plot_annotation(title=replace_exact_matches('fg10', nm)))
# (city_rh_lr <- main_res$p_lr_rh_city + plot_annotation(title=replace_exact_matches('rh', nm)))
# (city_pm25_lr <- main_res$p_lr_pm25_city + plot_annotation(title=replace_exact_matches('pm25', nm)))

# City case study: Chamba & Dharamshala

dharam_res <- format_results(output, outcome, city='Dharamshala')
dharam_lr <- get_city_legend(wrap_plots(dharam_res$p_lr_tmp_city, dharam_res$p_lr_tp_city,
                                        nrow=2)+plot_annotation(title='Lag-response: Dharamshala', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
dharam_er <- get_city_legend(wrap_plots(dharam_res$p_er_tmp_city, dharam_res$p_er_tp_city,
                                        nrow=2)+plot_annotation(title='Exposure-response: Dharamshala', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))

cham_res <- format_results(output, outcome, city='Chamba')
cham_lr <- get_city_legend(wrap_plots(cham_res$p_lr_tmp_city, cham_res$p_lr_tp_city, 
                                      nrow=2)+plot_annotation(title='Lag-response: Chamba', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
cham_er <- get_city_legend(wrap_plots(cham_res$p_er_tmp_city, cham_res$p_er_tp_city, 
                                      nrow=2)+plot_annotation(title='Exposure-response: Chamba', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))


# Main percentile look up
cen_tables <- list()
datalist <- lapply(locations, function(region) df[df$location==region,])

for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)-10
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_tables[[v]] <- round(cen, 4)
}
names(cen_tables) <- c('Maximum temperature', 'Precipitation') # Wind speed

#### EXPORT #####

# FIGURES
exports <- list('MF_1A_total_er_overall'=list(p=total_er_overall, w=9, h=5),
                'MF_1B_total_lr'=list(p=total_lr, w=8, h=6),
                'MF_1C_total_er'=list(p=total_er, w=8, h=6),
                'MF_1D_city_tmp_er'=list(p=city_tmp_er, w=9, h=5),
                'MF_1E_city_tp_er'=list(p=city_tp_er, w=9, h=5),
                # 'MF_1F_city_fg10_er'=list(p=city_fg10_er, w=9, h=5),
                'SF_1G_cham_lr'=list(p=cham_lr, w=8, h=6),
                'SF_1H_cham_er'=list(p=cham_er, w=8, h=6),
                'SF_1I_dharam_lr'=list(p=dharam_lr, w=8, h=6),
                'SF_1J_dharam_er'=list(p=dharam_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'figures/', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300, bg='white')
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"tables/ST_main_perc.xlsx"), overwrite = TRUE)

# Analysis 2 - City-specific models (city params) ####

# Get city-specific results
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- run_bestfit(df=sub, locations, outcome, params=final_params[[l]], pooled=FALSE, city=l)
  results[[l]] <- res
  print(l)
}

(spec_dharam_er <- results[['Dharamshala']]$results$p_er + plot_annotation(title=paste0('Exposure-response: Dharamshala')))
(spec_cham_er <- results[['Chamba']]$results$p_er + plot_annotation(title=paste0('Exposure-response: Chamba')))

#### EXPORT #####

# FIGURES
exports <- list('SF_2A_spec_dharam_er'=list(p=spec_dharam_er, w=8, h=6),
                'SF_2B_spec_cham_er'=list(p=spec_cham_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'figures/', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300, bg='white')
}


# Analysis 3 - Seasonal differences (pooled params) ####

# Winter/pre-monsoon model (months=5)
df_win <- df %>% filter(season %in% c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)'))
win <- run_bestfit(df_win, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$winter <- win
win_res <- win$results

### Winter - Overall
win_res$p_er_overall <- wrap_plots(lapply(win_res$p_er_overall, function(p) p + theme(legend.position = "none")), nrow=1)
win_er_overall <- wrap_plots(get_title('Winter'), win_res$p_er_overall, nrow=2, heights=c(1, 40))
win_lr <- wrap_plots(get_title('Lag-response: Winter'), win_res$p_lr, nrow=2, heights=c(1, 40))
win_er <- wrap_plots(get_title('Exposure-response: Winter'), win_res$p_er, nrow=2, heights=c(1, 40))

### Winter - Percentile look up
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
names(cen_tables) <- c('Maximum temperature', 'Precipitation') #, 'Wind speed'


# Monsoon/Post-monsoon model (months=7)
df_sum <- df %>% filter(season %in% c('Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)'))
sum <- run_bestfit(df_sum, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$summer <- sum
sum_res <- sum$results

### Summer - Overall
sum_er_overall <- wrap_plots(get_title('Summer'), sum_res$p_er_overall, nrow=2, heights=c(1, 40))
sum_lr<- wrap_plots(get_title('Lag-response: Summer'), sum_res$p_lr, nrow=2, heights=c(1, 40))
sum_er<- wrap_plots(get_title('Exposure-response: Summer'), sum_res$p_er, nrow=2, heights=c(1, 40))

### Summer - Percentile look up
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
names(cen_tables) <- c('Maximum temperature', 'Precipitation') #, 'Wind speed'


# MERGE PLOTS FOR COMPARISON
(seas_overall <- wrap_plots(win_er_overall, sum_er_overall, nrow=2))
(seas_lr <- wrap_plots(win_lr, plot_spacer(), sum_lr, nrow=1, widths = c(1, 0.08, 1)))
(seas_er <- wrap_plots(win_er, plot_spacer(), sum_er, nrow=1, widths = c(1, 0.08, 1)))


#### EXPORT #####
exports <- list('MF_3A_seas_overall'=list(p=seas_overall, w=8, h=6),
                'MF_3B_seas_lr'=list(p=seas_lr, w=15, h=6),
                'MF_3C_seas_er'=list(p=seas_er, w=15, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'figures/', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300)
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_win_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_win_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"tables/ST_winter_perc.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
for (name in names(cen_sum_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_sum_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"tables/ST_summer_perc.xlsx"), overwrite = TRUE)


# Analysis 4 - Population differences (pooled params) ####

# # Get total tables and id relevant lags/exposure
# x <- meta_checks(output$meta_tmp, 1.001)
# meta_checks(output$meta_tp, 1.001)
# meta_checks(output$meta_pm25, 1.000)

# Multiple outcomes
outcomes <- h_cols[1:8]
er0_tmp <- list()
er3_tmp <- list()
er10_tp <- list()
for (i in outcomes){
  print(i)
  res <- run_bestfit(df, locations, i, params=final_params$pooled, pooled=TRUE, city=NULL)
  results[[i]] <- res
  pop_res <- res$results
  er0_tmp[[i]] <- pop_res$lags$`0`$p_er_tmp$p_blup + labs(title=replace_exact_matches(i, nm))
  er3_tmp[[i]] <- pop_res$lags$`3`$p_er_tmp$p_blup + labs(title=replace_exact_matches(i, nm))
  er10_tp[[i]] <- pop_res$lags$`10`$p_er_tp$p_blup + labs(title=replace_exact_matches(i, nm))
}

(er0_tmp <- (wrap_plots(er0_tmp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Exposure-response (Lag=0)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er3_tmp <- (wrap_plots(er3_tmp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Exposure-response (Lag=3)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er10_tp <- (wrap_plots(er10_tp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Exposure-response (Lag=10)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))

all_outcomes <- list(er0_tmp=er0_tmp, er3_tmp=er3_tmp, er10_tp=er10_tp) #, er0_fg10=er0_fg10

#### EXPORT ####

# Population plots
names(all_outcomes) <- paste0("disag_", names(all_outcomes))
for (n in names(all_outcomes)) {
  ggsave(filename=paste0(path_report, 'figures/MF_4_', n, '.png'), plot=all_outcomes[[n]], width=8, height=5, dpi=300)
}

# Analysis 5 - Hamirpur emergency ward (emergency params) ####

# Set outcome and site
outcome <- 'emergency'
l <- 'Hamirpur'

# Main results
sub <- df %>% filter(location==l)
res <- run_bestfit(sub, locations, outcome, params=final_params$emergency, pooled=FALSE, city=l)
results$emergency <- res
emg_res <- res$results

# Overall
emg_res$p_er_overall <- wrap_plots(lapply(emg_res$p_er_overall, function(p) p + theme(legend.position = "none")), nrow=1)
(emg_er_overall <- emg_res$p_er_overall + plot_annotation(title='Overall: Hamirpur emergency ward', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
(emg_lr <- emg_res$p_lr + plot_annotation(title='Lag-response: Hamirpur emergency ward'))
(emg_er <- emg_res$p_er + plot_annotation(title='Exposure-response: Hamirpur emergency ward'))


#### EXPORT ####

exports <- list('MF_5A_emg_overall'=list(p=emg_er_overall, w=9, h=5),
                'MF_5B_emg_lr'=list(p=emg_lr, w=8, h=6),
                'MF_5C_emg_er'=list(p=emg_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'figures/', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300)
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







