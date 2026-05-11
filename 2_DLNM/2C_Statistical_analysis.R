### ------------------------------------------------------------------------ ### 
###  DNLM - Multiple city meta- and uni-analysis (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/offne/Documents/GitHub/hisp-climate-health-forecast/')


####  DEPENDENCIES  ####
source('2_DLNM/DLNM - Dependencies.R')

# Load packages
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork); library(openxlsx);
library(colorspace); library(png); library(grid); library(reshape2); library(cowplot)
library(TTR)

# Function to get empty plot with title, for nested wrap_plot work around
get_title <- function(title, size=16){
  
  ggplot() +
    labs(title = title) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5, size = size, face = "bold"))
  
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
load(paste0(path_project, '/ind_sens_params.RData'))
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData'))
df <- df_dlnm; rm(df_dlnm)

# Set location info by elev
locations <- as.character(df %>% distinct(location, elev) %>% arrange(elev) %>% pull(location))
m <- length(locations)

# Store results
results <- list()

# Set outcome outcome and exposures
outcome <- 'total_patients'
variables <- c('t2m', 'tp')

# Analysis 1 - Overall (pooled params) ----

# Main results
res <- run_bestfit(df, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$overall <- res
output   <- res$output
main_res <- res$results

# Overall exposure-response by lag
# (total_er_overall <- main_res$p_er_overall + plot_annotation(title=replace_exact_matches(outcome, nm)))
# (total_er <- main_res$p_er + plot_annotation(title='Exposure-response'))


# Total patient Contours #####

tot_contours <- function(sig=FALSE, max_dev_tmp=NULL, max_dev_tp=NULL){
  
  # CONTOURS
  plot_contours <- function(meta, locations, xlab, title, reference='50%', sig=sig, 
                            n_col=2, max_dev=NULL) {
    
    args <- list(meta=meta, xlab=xlab, reference=reference, sig=sig)
    
    # Shared limits
    dfs <- c(
      list(ALL=do.call(city_con_res, c(args, list(city='ALL', title='Overall')))$df),
      setNames(lapply(locations, function(l) do.call(city_con_res, c(args, list(city=l, title=l)))$df), locations)
    )
    
    if (is.null(max_dev)) max_dev <- max(abs(pretty(range(do.call(rbind, dfs)$z, na.rm=TRUE), n=15) - 1))
    shared_lims <- c(1 - max_dev, 1 + max_dev)
    shared_brks <- pretty(shared_lims, n=12)
    args <- c(args, list(shared_lims=shared_lims, shared_brks=shared_brks))
    
    # Plots
    contours <- list()
    for (l in locations) {
      contours[[l]] <- do.call(city_con_res, c(args, list(city=l, title=l)))$p
    }
    
    # Layout
    n_rows <- ceiling(length(contours) / n_col)
    contours_mod <- lapply(seq_along(contours), function(i) {
      row <- ceiling(i / n_col); col <- ((i - 1) %% n_col) + 1
      contours[[i]] +
        theme(
          axis.title.x=if(row == n_rows) element_text() else element_blank(),
          axis.title.y=if(col == 1)      element_text() else element_blank()
        )
    })
    
    wrap_plots(contours_mod, ncol=n_col) +
      plot_layout(guides='collect') +
      plot_annotation(title=title,
                      theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))) &
      theme(legend.position='none')  
  }
  
  # --- Compute shared max_dev across pooled AND city-specific ---
  compute_max_dev <- function(meta, locations, reference, args_base) {
    all_dfs <- c(
      list(do.call(city_con_res, c(args_base, list(city='ALL', title='Overall')))$df),
      lapply(locations, function(l) do.call(city_con_res, c(args_base, list(city=l, title=l)))$df)
    )
    max(abs(pretty(range(do.call(rbind, all_dfs)$z, na.rm=TRUE), n=15) - 1))
  }
  
  args_tmp <- list(meta=output$meta_tmp, xlab='Temperature (%)',   reference='50%', sig=sig)
  args_tp  <- list(meta=output$meta_tp,  xlab='Precipitation (%)', reference='20%', sig=sig)
  
  if (is.null(max_dev_tmp)) max_dev_tmp <- compute_max_dev(output$meta_tmp, locations, '50%', args_tmp)
  if (is.null(max_dev_tp)) max_dev_tp <- compute_max_dev(output$meta_tp, locations, '50%', args_tmp)

  # Shared limits for pooled plots
  make_shared <- function(max_dev) {
    lims <- c(1 - max_dev, 1 + max_dev)
    list(shared_lims=lims, shared_brks=pretty(lims, n=10))
  }
  
  # Overall pooled plots
  tmp <- do.call(city_con_res, c(args_tmp, make_shared(max_dev_tmp), list(city='ALL', title=NULL)))$p 
  tp  <- do.call(city_con_res, c(args_tp,  make_shared(max_dev_tp),  list(city='ALL', title=NULL)))$p
  
  # Extract legend from one plot using cowplot
  shared_legend <- get_legend(tmp + theme(legend.position='bottom', legend.text=element_text(angle=45, hjust=1)))
  
  # Get title
  tmp <- wrap_plots(get_title('Overall temperature'), tmp, nrow=2, heights=c(1, 40)) + 
    theme(legend.position='none')
  tp  <- wrap_plots(get_title('Overall precipitation'), tp, nrow=2, heights=c(1, 40)) + 
    theme(legend.position='none')
  
  # City-specific plots
  tmp_l <- plot_contours(meta=output$meta_tmp, locations=locations,
                         xlab='Temperature (%)', title='Maximum Temperature',
                         reference='50%', sig=sig, max_dev=max_dev_tmp) 
  tp_l  <- plot_contours(meta=output$meta_tp, locations=locations,
                         xlab='Precipitation (%)', title='Total Precipitation',
                         reference='20%', sig=sig, max_dev=max_dev_tp)
  tmp_l <- wrap_plots(get_title('Site-specific temperature'), tmp_l, nrow=2, heights=c(1, 40))
  tp_l  <- wrap_plots(get_title('Site-specific precipitation'), tp_l, nrow=2, heights=c(1, 40))
  
  con_tmp <- wrap_plots(list(tmp, tmp_l), ncol=1) + plot_layout(heights=c(1, 2))
  con_tp  <- wrap_plots(list(tp,  tp_l),  ncol=1) + plot_layout(heights=c(1, 2))
  combined <- wrap_plots(list(con_tmp, con_tp), nrow=1)
  
  plot_grid(combined, shared_legend, ncol=1, rel_heights=c(1, 0.08))
  
}

dev <- 0.16
contour_all <- tot_contours(sig=FALSE, max_dev_tmp=dev, max_dev_tp=dev)
contour_all_sig <- tot_contours(sig=TRUE, max_dev_tmp=dev, max_dev_tp=dev)

ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_ALL.png'), plot=contour_all, width=7, height=11, dpi=300, bg='white')
ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_ALL_sig.png'), plot=contour_all_sig, width=7, height=11, dpi=300, bg='white')


# Total patient Exposure response ####

tmp_lags <- list('ALL'= c(0,3), #c(0,3,7),
                 'Una'= c(0,3), #c(4),
                 'Hamirpur'=c(0,3), #c(0,2,8),
                 'Nahan'=c(0,3), #c(0,8),
                 'Chamba'=c(0,3), #c(0,4,8),
                 'Kullu'=c(0,3), #c(0,4,7),
                 'Dharamshala'=c(0,3), #c(0,8),
                 'Solan'=c(0,3), #c(0,4,7,10),
                 'Reckong Peo'=c(0,3)) #c(0,8))

tp_lags <- list('ALL'=c(3,10),
                'Una'=c(3,10), #c(),
                'Hamirpur'=c(3,10), #c(),
                'Nahan'=c(3,10), #c(),
                'Chamba'=c(3,10), #c(3,10),
                'Kullu'=c(3,10), #c(),
                'Dharamshala'=c(3,10), #c(1,9),
                'Solan'=c(3,10), #c(),
                'Reckong Peo'=c(3,10)) #c())



exp <- function(output, city, plot_key){
  
  city_tmp_lags <- tmp_lags[[city]]
  city_tp_lags  <- tp_lags[[city]]
  
  # Build tmp plots
  tmp_plots <- lapply(city_tmp_lags, function(lag) {
    city_exp_res(city=city, lag=lag, meta=output$meta_tmp,
                 col='#D74B3C', xlab='Temperature (%)', 
                 title=paste0('Lag=', lag), reference='50%')[[plot_key]]
  })
  # Build tp plots
  tp_plots <- lapply(city_tp_lags, function(lag) {
    city_exp_res(city=city, lag=lag, meta=output$meta_tp,
                 col='#1D98EB', xlab='Precipitation (%)', 
                 title=paste0('Lag=', lag), reference='0%')[[plot_key]]
  })
  
  all_plots <- c(tmp_plots, tp_plots)
  
  # Pad to exactly 5 with void placeholders
  while (length(all_plots) < 5) all_plots <- c(all_plots, list(ggplot() + theme_void()))
  
  # Strip legends
  all_plots <- lapply(all_plots, function(p) p + theme(legend.position='none'))
  
  return(all_plots)  # return list of exactly 5 plots
}

# --- Build label + 5 plots per row, assembled into one aligned grid ---
locations_order <- c('Chamba', 'Dharamshala', 'Solan', 'Hamirpur', 'Kullu', 
                     'Reckong Peo', 'Nahan', 'Una')

# Collect all rows as flat list of (label + 5 plots)
all_cells <- list()
for (l in locations_order) {
  
  # Centered city label as a single cell
  label_cell <- ggplot() +
    annotate('text', x=0.5, y=0.5, label=l, 
             fontface='bold', size=4.5, hjust=0.5, vjust=0.5) +
    theme_void()
  
  plot_cells <- exp(output, l, plot_key='p_city')  # list of exactly 5
  all_cells  <- c(all_cells, list(label_cell), plot_cells)  # 6 cells per row
}

# Assemble into grid: 6 cols (1 label + 5 plots), n rows
combined <- wrap_plots(all_cells, ncol=6, byrow=TRUE) +
  plot_layout(widths=c(10, rep(10, 5)))  # label col narrow, plot cols equal

# --- Pooled ---
exp_pool_plots <- exp(output, 'ALL', plot_key='p_blup')
exp_pool <- wrap_plots(exp_pool_plots, nrow=1, ncol=5)

# --- Unified legend ---
unified_legend <- ggplot(data.frame(x=1, y=1), aes(x, y)) +
  geom_line(aes(linetype='BLUP'),             colour='black') +
  geom_line(aes(linetype='Pooled'),           colour='black') +
  geom_line(aes(colour='Raw temperature'),    linewidth=1) +
  geom_line(aes(colour='Raw precipitation'),  linewidth=1) +
  scale_linetype_manual(name=NULL, values=c('BLUP'='dashed', 'Pooled'='solid')) +
  scale_colour_manual(name=NULL, values=c('Raw temperature'='#D74B3C', 'Raw precipitation'='#1D98EB')) +
  theme_void() +
  theme(legend.position='bottom',
        legend.direction='horizontal',
        legend.text=element_text(size=10))
shared_legend <- get_legend(unified_legend)


exp_pool <- plot_grid(exp_pool, shared_legend, ncol=1, rel_heights=c(40, 1))
exp_city <- plot_grid(combined, shared_legend, ncol=1, rel_heights=c(40, 1))


ggsave(filename=paste0(path_report, 'II. dlnm/fig_exp_pool.png'), 
       plot=exp_pool, width=15, height=4,  dpi=300, bg='white')
ggsave(filename=paste0(path_report, 'II. dlnm/fig_exp_city.png'), 
       plot=exp_city, width=10, height=15, dpi=300, bg='white')

# Population Contours & Exposure responses #####

pop_contours <- function(outcomes, output, tmp, tp, all_z_tmp, all_z_tp,max_dev=NULL, sig=FALSE){
  
  # --- Step 2: Compute shared limits across ALL outcomes ---
  make_shared <- function(all_z, max_dev=0.18) {
    if (is.null(max_dev)) max_dev <- max(abs(pretty(range(all_z, na.rm=TRUE), n=15) - 1))
    shared_lims <- c(1 - max_dev, 1 + max_dev)
    shared_brks <- pretty(shared_lims, n=12)
    list(shared_lims=shared_lims, shared_brks=shared_brks)
  }
  
  scale_tmp <- make_shared(all_z_tmp, max_dev=max_dev)
  scale_tp  <- make_shared(all_z_tp, max_dev=max_dev)
  
  # --- Step 3: Build plots using shared limits ---
  for (i in outcomes){
    args_tmp <- c(list(meta=output[[i]]$meta_tmp, xlab='Temperature (%)',   reference='50%', sig=sig), scale_tmp)
    args_tp  <- c(list(meta=output[[i]]$meta_tp,  xlab='Precipitation (%)', reference='20%', sig=sig), scale_tp)
    
    tmp[[i]] <- do.call(city_con_res, c(args_tmp, list(city='ALL', title=replace_exact_matches(i, nm))))$p
    tp[[i]]  <- do.call(city_con_res, c(args_tp,  list(city='ALL', title=replace_exact_matches(i, nm))))$p
  }
  
  # --- Step 4: Combine ---
  con_tmp <- wrap_plots(tmp, ncol=4) + 
    plot_layout(guides='collect') + 
    plot_annotation(title='Overall maximum temperature', 
                    theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))) &
    theme(legend.position='none',
          legend.text=element_text(angle=45, hjust=1))  
  con_tp <- wrap_plots(tp, ncol=4) + 
    plot_layout(guides='collect') + 
    plot_annotation(title='Overall total precipitation', 
                    theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))) &
    theme(legend.position='none',
          legend.text=element_text(angle=45, hjust=1))
  
  con_tmp <- wrap_plots(get_title('Pooled maximum temperature'), con_tmp, nrow=2, heights=c(1, 40))
  con_tp <- wrap_plots(get_title('Pooled 7-day cummulative precipitation'), con_tp, nrow=2, heights=c(1, 40))

  # Extract legend from one plot using cowplot
  shared_legend <- get_legend(
    tmp[[1]] + theme(legend.position='bottom',
                     legend.text=element_text(angle=45, hjust=1))
  )
  
  combined <- wrap_plots(list(con_tmp, con_tp), nrow=2) + plot_annotation(theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
  plot_grid(combined, shared_legend, ncol=1, rel_heights=c(1, 0.08))

}

outcomes <- h_cols[1:8]
output <- list()
tmp    <- list()
tp     <- list()
pop_exp <- list()

# --- Step 1: Run models and collect all z-values first ---
all_z_tmp <- c()
all_z_tp  <- c()

for (i in outcomes){
  print(i)
  res <- run_bestfit(df, locations, i, params=final_params$pooled, pooled=TRUE, city=NULL)
  output[[i]] <- res$output
  
  args_tmp <- list(meta=output[[i]]$meta_tmp, xlab='Temperature (%)',   reference='50%', sig=F)
  args_tp  <- list(meta=output[[i]]$meta_tp,  xlab='Precipitation (%)', reference='20%', sig=F)
  
  tp_lags$ALL <- c(3,10)
  if(i %in% c('M', 'age25_64', 'over_65')){
    tp_lags$ALL <- c(10)
  } else if (i == 'age5_14'){
    tp_lags$ALL <- c()
  }
  
  # exp_pool_plots <- exp(output[[i]] , 'ALL', plot_key='p_blup')
  # exp_pool_plots <- wrap_plots(exp_pool_plots, nrow=1)
  # pop_exp[[i]] <- wrap_plots(get_title(replace_exact_matches(i, nm), size=14), exp_pool_plots, nrow=2, heights=c(1, 40))
  # With:
  exp_pool_plots <- exp(output[[i]], 'ALL', plot_key='p_blup')
  
  # Pad to exactly 5 slots
  while (length(exp_pool_plots) < 5) exp_pool_plots <- c(exp_pool_plots, list(ggplot() + theme_void()))
  exp_pool_plots <- exp_pool_plots[1:5]  # ensure never more than 5
  
  exp_pool_plots <- wrap_plots(exp_pool_plots, nrow=1, ncol=5)
  pop_exp[[i]] <- wrap_plots(get_title(replace_exact_matches(i, nm), size=14), exp_pool_plots, nrow=2, heights=c(1, 40))
  
  # Collect z ranges across all outcomes
  all_z_tmp <- c(all_z_tmp, do.call(city_con_res, c(args_tmp, list(city='ALL', title=i)))$df$z)
  all_z_tp  <- c(all_z_tp,  do.call(city_con_res, c(args_tp,  list(city='ALL', title=i)))$df$z)
}

# Get contours
contour_pop <- pop_contours(outcomes, output, tmp, tp, all_z_tmp, all_z_tp, max_dev=0.18, sig=FALSE)
contour_pop_sig <- pop_contours(outcomes, output, tmp, tp, all_z_tmp, all_z_tp, max_dev=0.18, sig=TRUE)

# Get exp-response
pop_exp_all <- wrap_plots(pop_exp, ncol=2)
pop_exp_all <- plot_grid(pop_exp_all, shared_legend, ncol=1, rel_heights=c(40, 1))

ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_pop.png'), plot=contour_pop, width=10, height=12, dpi=300, bg='white')
ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_pop_sig.png'), plot=contour_pop_sig, width=10, height=12, dpi=300, bg='white')
ggsave(filename=paste0(path_report, 'II. dlnm/fig_pop_exp.png'), plot=pop_exp_all, width=15, height=10, dpi=300, bg='white')


#########################


# # Main exposure-response, by city
# (city_tmp_er <- get_city_legend(pred='mx2t', plot=main_res$p_er_tmp_city + plot_annotation(title='Exposure-response, by city: Temperature')))
# (city_tp_er <- get_city_legend(pred='tp', plot=main_res$p_er_tp_city + plot_annotation(title='Exposure-response, by city: Precipitation')))

# Main lag-response, by city
# (city_tmp_lr <- main_res$p_lr_tmp_city + plot_annotation(title=replace_exact_matches('wb', nm)))
# (city_tp_lr <- main_res$p_lr_tp_city + plot_annotation(title=replace_exact_matches('tp', nm)))

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
exports <- list('1A_total_er_overall'=list(p=total_er_overall, w=9, h=5),
                '1B_total_lr'=list(p=total_lr, w=8, h=6),
                '1C_total_er'=list(p=total_er, w=8, h=6),
                '1D_con_tmp'=list(p=con_tmp, w=5, h=4),
                '1D_con_tp'=list(p=con_tp, w=5, h=4),
                '1D_con_tmp_sig'=list(p=con_tmp_sig, w=5, h=4),
                '1D_con_tp_sig'=list(p=con_tp_sig, w=5, h=4),
                '1D_city_tmp_er'=list(p=city_tmp_er, w=9, h=5),
                '1E_city_tp_er'=list(p=city_tp_er, w=9, h=5),
                # '1F_city_fg10_er'=list(p=city_fg10_er, w=9, h=5),
                '1G_cham_lr'=list(p=cham_lr, w=8, h=6),
                '1H_cham_er'=list(p=cham_er, w=8, h=6),
                '1I_dharam_lr'=list(p=dharam_lr, w=8, h=6),
                '1J_dharam_er'=list(p=dharam_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300, bg='white')
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"II. dlnm/tab_main_perc.xlsx"), overwrite = TRUE)

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
exports <- list('2A_spec_dharam_er'=list(p=spec_dharam_er, w=8, h=6),
                '2B_spec_cham_er'=list(p=spec_cham_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300, bg='white')
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
exports <- list('3A_seas_overall'=list(p=seas_overall, w=8, h=6),
                '3B_seas_lr'=list(p=seas_lr, w=15, h=6),
                '3C_seas_er'=list(p=seas_er, w=15, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300)
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_win_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_win_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"II. dlnm/tab_winter_perc.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
for (name in names(cen_sum_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_sum_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"II. dlnm/tab_summer_perc.xlsx"), overwrite = TRUE)


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
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_4_', n, '.png'), plot=all_outcomes[[n]], width=8, height=5, dpi=300)
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

exports <- list('5A_emg_overall'=list(p=emg_er_overall, w=9, h=5),
                '5B_emg_lr'=list(p=emg_lr, w=8, h=6),
                '5C_emg_er'=list(p=emg_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300)
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







