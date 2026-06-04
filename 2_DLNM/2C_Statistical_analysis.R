### ------------------------------------------------------------------------ ### 
###  DNLM - Multiple city meta- and uni-analysis (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

####  DEPENDENCIES  ####
source('2_DLNM/2_Utilities.R')

# Load packages
library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
library(ggplot2); library(dplyr); library(purrr); library(patchwork); library(openxlsx);
library(colorspace); library(png); library(grid); library(reshape2); library(cowplot)
library(TTR); library(lmtest); library(flextable); library(officer)

# Function to create city specific legend from dummy set
get_city_legend <- function(plot, pred='ALL'){
  
  legend_df <- data.frame(x=1, y=1, 
                          predictor=factor(c("mx2t","tp","fg10", 'pm25', 'rh', "Pooled","BLUP"), 
                                             levels=c("mx2t","tp","fg10",'pm25', 'rh',"Pooled","BLUP")), 
                          colour=c("#D74B3C","#1D98EB","#8968CD", '#E569BF', '#A6B727', "black","black"), 
                          linetype=c("solid","solid","solid","solid","solid","solid","dashed"), 
                          size=c(1.2,1.2,1.2,1.2,1.2,0.7,0.7))
  legend_df <- legend_df %>% filter(!predictor %in% c('fg10', "pm25", "rh"))
  
  if(pred != 'ALL'){
    legend_df <- subset(legend_df, predictor %in% c(pred, "Pooled", "BLUP"))
  }
  
  vec <- legend_df$predictor

  p <- ggplot(legend_df, aes(x, y, colour=predictor, linetype=predictor, size=predictor)) +
    geom_line() +
    scale_colour_manual(
      values=setNames(legend_df$colour, legend_df$predictor),
      breaks=c(vec[-((length(vec)-1):length(vec))]),
      labels=replace_exact_matches(as.data.frame(c(vec[-((length(vec)-1):length(vec))])), nm)[,1],
      guide=guide_legend(title='Raw', nrow=length(c(vec[-((length(vec)-1):length(vec))])), 
                           order=1, 
                           override.aes=list(linewidth=rep(1, length(c(vec[-((length(vec)-1):length(vec))])))))) +
    scale_linetype_manual(
      values=setNames(legend_df$linetype, legend_df$predictor),
      breaks=c(vec[((length(vec)-1):length(vec))]),
      guide=guide_legend(title='Meta-analysis', nrow=2, order=2, 
                           override.aes=list(colour="black", 
                                               linewidth=c(0.7, 0.7)))) +
    scale_size_manual(
      values=setNames(legend_df$size, legend_df$predictor),
      breaks=c(vec[((length(vec)-1):length(vec))]),
      guide="none") +  theme_void() +
    theme(
      legend.position="right",
      legend.box="vertical"
    )
  
  city_legend <- get_legend(p)
  
  res <- plot_grid(plot, city_legend, ncol=2, rel_widths =c(10, 2))
  
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

# META REGRESSION: Elevation and population density
get_heterogeneity <- function(mv_fit) {
  s  <- summary(mv_fit)
  q  <- s$qstat$Q[".all"]
  df <- s$qstat$df[".all"]
  p  <- s$qstat$pvalue[".all"]
  i2 <- max(0, (q - df) / q) * 100
  list(
    AIC = round(s$AIC, 1),
    Q   = round(q, 1),
    df  = df,
    p   = ifelse(p < 0.001, "<0.001", round(p, 3)),
    I2  = round(i2, 1)
  )
}

get_lr_test <- function(mv_null, mv_mod) {
  
  likelihood_ratio_test <- lrtest(mv_null, mv_mod)
  lr_p <- likelihood_ratio_test$`Pr(>Chisq)`[2]
  list(
    stat = likelihood_ratio_test$LogLik[2],
    df   = likelihood_ratio_test$Df[2],
    p    = ifelse(lr_p < 0.001, "<0.001", round(lr_p, 3))
  )
}

make_meta_reg_table <- function(meta, predictors_list) {
  
  rows <- list()
  
  # Intercept-only model (already fitted and stored in meta$mv)
  het <- get_heterogeneity(meta$mv)
  
  rows[["Intercept only"]] <- data.frame(
    `Meta-predictors` = "Intercept only",
    LR_stat = NA, LR_df = NA, LR_p = NA,
    AIC    = het$AIC,
    Q_stat = het$Q,  Q_df = het$df, Q_p = het$p,
    I2     = het$I2,
    check.names = FALSE
  )
  
  # Single meta-predictor models
  for (predictor in names(predictors_list)) {
    
    pred <- sapply(locations, function(loc) {
      if (predictor != "elev") {
        mean(df[[predictor]][df$location == loc],na.rm = TRUE)
      } else {
        unique(df[[predictor]][df$location == loc])
      }
    })
    
    mv_mod <- tryCatch(
      mvmeta(meta$ymat ~ pred, meta$Slist, method = "ml"),
      error = function(e) NULL
    )
    if (is.null(mv_mod)) next
    
    lr  <- get_lr_test(meta$mv, mv_mod)
    het <- get_heterogeneity(mv_mod)
    
    rows[[predictor]] <- data.frame(
      `Meta-predictors` = predictors_list[[predictor]],
      LR_stat = lr$stat, LR_df = lr$df, LR_p = lr$p,
      AIC    = het$AIC,
      Q_stat = het$Q,  Q_df = het$df, Q_p = het$p,
      I2     = het$I2,
      check.names = FALSE
    )
  }
  
  results <-  do.call(rbind, rows)
  row.names(results) <- NULL
  
  return(results)
}

#### MAIN CODE ####

# Load data and parameters
load(paste0(path_project, 'Results/ind_sens_params.RData'))
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData'))
df <- df_dlnm; rm(df_dlnm)

# Set location info by elev
locations <- as.character(df %>% distinct(location, elev) %>% arrange(elev) %>% pull(location))
m <- length(locations)

# Store results
results <- list()

# Set outcome outcome and exposures
outcome <- 'total_patients'

# Part 1 - Pooled model ----

# Main results
res <- run_bestfit(df, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$overall <- res
output   <- res$output
main_res <- res$results

# Get contour plots
tot_contours <- function(sig=FALSE) {
  
  get_max_dev <- function(meta, reference) {
    args <- list(meta=meta, reference=reference, sig=sig)
    all_dfs <- c(
      list(do.call(city_con_res, c(args, list(city="ALL", title="Overall")))$df),
      lapply(locations, function(l) do.call(city_con_res, c(args, list(city=l, title=l)))$df)
    )
    combined <- do.call(rbind, all_dfs)
    brks <- pretty(range(combined$z, na.rm=TRUE), n=8)
    max(abs(brks - 1))
  }
  
  make_plots <- function(meta, xlab, reference, max_dev=NULL) {
    
    if(is.null(max_dev)){ max_dev <- get_max_dev(meta, reference) }
    shared_lims <- c(0.999 - max_dev, 1 + max_dev) * 1.001
    shared_brks <- pretty(shared_lims, n=8)
    
    args <- list(meta=meta, xlab=xlab, reference=reference, sig=sig,
                 shared_lims=shared_lims, shared_brks=shared_brks)
    
    plots <- list()
    plots[["ALL"]] <- do.call(city_con_res, c(args, list(city="ALL", title="Pooled")))$p
    for (l in locations) {
      plots[[l]] <- do.call(city_con_res, c(args, list(city=l, title=l)))$p
    }
    
    plots
  }
  
  list(
    tmp=make_plots(meta=output$meta_tmp, xlab="Temperature (%)", reference="50%", max_dev=0.2),
    tp =make_plots(meta=output$meta_tp,  xlab="Precipitation (%)", reference="20%", max_dev=0.2)
  )
}
cont <- tot_contours(FALSE)
cont_sig <- tot_contours(TRUE)

# Check significant lags before setting relevant exp-resp
wrap_plots(cont_sig$tmp, ncol=3) +
  plot_layout(guides='collect') +
  plot_annotation(title='Significant exposure-lag-response',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))) &
  theme(legend.position='bottom')

# Get non-linear plots
tot_exp <- function() {
  
  make_lag_plots <- function(meta, lags_list, col, xlab, reference) {
    
    plots <- list()
    
    all_locations <- c("ALL", locations)
    for (l in all_locations) {
      plot_key  <- 'p_blup' #if (l == "ALL") "p_blup" else "p_city"
      city_lags <- lags_list[[l]]
      plots[[l]] <- lapply(city_lags, function(lag) {
        city_exp_res(city=l, lag=lag, meta=meta,
                     col=col, xlab=xlab,
                     title=paste0("Lag=", lag),
                     reference=reference)[[plot_key]]
      })
    }
    
    plots
  }
  
  list(
    tmp=make_lag_plots(
      meta     =output$meta_tmp,
      lags_list=tmp_lags,
      col      ="grey30",
      xlab     ="Temperature (%)",
      reference="50%"
    ),
    tp=make_lag_plots(
      meta     =output$meta_tp,
      lags_list=tp_lags,
      col      ="grey30",
      xlab     ="Precipitation (%)",
      reference="0%"
    )
  )
}
tmp_lags <- list('ALL'= c(0,3),
                 'Una'= c(0,4), 
                 'Hamirpur'=c(0,10), 
                 'Nahan'=c(0,9),
                 'Chamba'=c(0,4),
                 'Kullu'=c(0,4),
                 'Dharamshala'=c(0,3),
                 'Solan'=c(0,10),
                 'Reckong Peo'=c(0,4))
tp_lags <- list('ALL'=c(3,10),
                'Una'=c(3,10), 
                'Hamirpur'=c(3,10),
                'Nahan'=c(0,10),
                'Chamba'=c(3,10),
                'Kullu'=c(3,10),
                'Dharamshala'=c(1,9), 
                'Solan'=c(3,10),
                'Reckong Peo'=c(3,7))

exp <- tot_exp()


# Export
saveRDS(output, file = "2_DLNM/Results/pooled_model.rds")
saveRDS(list(res=main_res, cont=cont, exp=exp), file = "2_DLNM/Results/pooled_plots.RData")

# Part 2 - Meta-predictor model ####

# Define predictors
predictors_list <- list(
  elev = "Altitude",
  popu = "Population density",
  rh   = "Relative humidity",
  pm25   = "Air pollution"
  
)
(tbl_tmp <- make_meta_reg_table(meta=output$meta_tmp, predictors_list))
(tbl_tp  <- make_meta_reg_table(output$meta_tp,  predictors_list))

# Export
saveRDS(list(tbl_tmp=tbl_tmp, tbl_tp=tbl_tp), file = "2_DLNM/Results/meta_tables.RData")

# Part 3 - Subgroup models ####

outcomes <- h_cols[1:8]
output   <- list()

# Shared z-value ranges across outcomes
all_z_tmp <- c()
all_z_tp  <- c()

# Run models and collect z-values
for (i in outcomes) {
  print(i)
  res <- run_bestfit(df, locations, i, params=final_params$pooled, pooled=TRUE, city=NULL)
  output[[i]] <- res$output
  
  args_tmp <- list(meta=output[[i]]$meta_tmp, xlab="Temperature (%)",   reference="50%", sig=FALSE)
  args_tp  <- list(meta=output[[i]]$meta_tp,  xlab="Precipitation (%)", reference="20%", sig=FALSE)
  
  all_z_tmp <- c(all_z_tmp, do.call(city_con_res, c(args_tmp, list(city="ALL", title=i)))$df$z)
  all_z_tp  <- c(all_z_tp,  do.call(city_con_res, c(args_tp,  list(city="ALL", title=i)))$df$z)
}

# Shared scales
make_shared <- function(all_z, max_dev=0.20) {
  if (is.null(max_dev)) max_dev <- max(abs(pretty(range(all_z, na.rm=TRUE), n=8) - 1))
  shared_lims <- c(0.999 - max_dev, 1 + max_dev) * 1.001
  shared_brks <- pretty(shared_lims, n=8)
  list(shared_lims=shared_lims, shared_brks=shared_brks)
}
scale_tmp <- make_shared(all_z_tmp, max_dev=0.20)
scale_tp  <- make_shared(all_z_tp,  max_dev=0.20)

# Collect contour and exposure plots per outcome
pooled_meta <- list(tmp = list(), tp = list())
cont_outcomes <- list(tmp=list(), tp=list())
exp_outcomes  <- list(tmp=list(), tp=list())

for (i in outcomes) {
  
  # Pooled results
  pooled_meta$tmp[[i]] <- output[[i]]$meta_tmp
  pooled_meta$tp[[i]]  <- output[[i]]$meta_tp
  
  # Contours
  cont_outcomes$tmp[[i]] <- do.call(city_con_res, c(
    list(meta=output[[i]]$meta_tmp, xlab="Temperature (%)", reference="50%", sig=FALSE),
    scale_tmp,
    list(city="ALL", title=replace_exact_matches(i, nm))
  ))$p
  cont_outcomes$tp[[i]] <- do.call(city_con_res, c(
    list(meta=output[[i]]$meta_tp, xlab="Precipitation (%)", reference="20%", sig=FALSE),
    scale_tp,
    list(city="ALL", title=replace_exact_matches(i, nm))
  ))$p
  
  # Exposure-lag-response plots
  exp_outcomes$tmp[[i]] <- lapply(tmp_lags$ALL, function(lag) {
    city_exp_res(city="ALL", lag=lag, meta=output[[i]]$meta_tmp,
                 col="grey30", xlab="Temperature (%)",
                 title=paste0("Lag=", lag), reference="50%")$p_blup
  })
  exp_outcomes$tp[[i]] <- lapply(tp_lags$ALL, function(lag) {
    city_exp_res(city="ALL", lag=lag, meta=output[[i]]$meta_tp,
                 col="grey30", xlab="Precipitation (%)",
                 title=paste0("Lag=", lag), reference="0%")$p_blup
  })
}

saveRDS(output, 
        file = "2_DLNM/Results/subgroup_models.rds")
saveRDS(list(pooled_meta=pooled_meta, cont_outcomes=cont_outcomes, exp_outcomes=exp_outcomes), 
        file = "2_DLNM/Results/subgroup_plots.RData")


# # Check significant lags before setting relevant exp-resp
# wrap_plots(cont_sig$tmp, ncol=3) +
#   plot_layout(guides='collect') +
#   plot_annotation(title='Significant exposure-lag-response',
#                   theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))) &
#   theme(legend.position='bottom')


