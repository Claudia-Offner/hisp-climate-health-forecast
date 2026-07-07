### ------------------------------------------------------------------------ ### 
###  DNLM – Pooled and subgroup models with meta predictors (HISP)
### ------------------------------------------------------------------------ ### 

rm(list=ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

####  DEPENDENCIES  ####
source('2_DLNM/2_Utilities.R')

library(dlnm); library(splines); library(mvmeta); library(MASS); library(TTR)
library(ggplot2); library(dplyr); library(purrr); library(patchwork)
library(openxlsx); library(colorspace); library(grid); library(reshape2)
library(cowplot); library(lmtest); library(flextable); library(officer)

## Legend panel for city-level plots
get_city_legend <- function(plot, pred='ALL') {
  legend_df <- data.frame(
    x        =1, y=1,
    predictor=factor(c("mx2t", "tp", "fg10", "pm25", "rh", "Pooled", "BLUP"),
                       levels=c("mx2t", "tp", "fg10", "pm25", "rh", "Pooled", "BLUP")),
    colour   =c("#D74B3C", "#1D98EB", "#8968CD", "#E569BF", "#A6B727", "black", "black"),
    linetype =c("solid", "solid", "solid", "solid", "solid", "solid", "dashed"),
    size     =c(1.2, 1.2, 1.2, 1.2, 1.2, 0.7, 0.7)
  ) |>
    filter(!predictor %in% c("fg10", "pm25", "rh"))
  
  if (pred != 'ALL')
    legend_df <- subset(legend_df, predictor %in% c(pred, "Pooled", "BLUP"))
  
  vec      <- legend_df$predictor
  raw_vec  <- vec[-((length(vec) - 1):length(vec))]
  meta_vec <- vec[(length(vec) - 1):length(vec)]
  
  p <- ggplot(legend_df, aes(x, y, colour=predictor, linetype=predictor, size=predictor)) +
    geom_line() +
    scale_colour_manual(
      values=setNames(legend_df$colour,   legend_df$predictor),
      breaks=raw_vec,
      labels=replace_exact_matches(as.data.frame(c(raw_vec)), nm)[, 1],
      guide =guide_legend(title='Raw', nrow=length(raw_vec), order=1,
                            override.aes=list(linewidth=rep(1, length(raw_vec))))
    ) +
    scale_linetype_manual(
      values=setNames(legend_df$linetype, legend_df$predictor),
      breaks=meta_vec,
      guide =guide_legend(title='Meta-analysis', nrow=2, order=2,
                            override.aes=list(colour="black", linewidth=c(0.7, 0.7)))
    ) +
    scale_size_manual(
      values=setNames(legend_df$size, legend_df$predictor),
      breaks=meta_vec,
      guide ="none"
    ) +
    theme_void() +
    theme(legend.position="right", legend.box="vertical")
  
  plot_grid(plot, get_legend(p), ncol=2, rel_widths=c(10, 2))
}

## Wrapper: fit best model for a given outcome
run_bestfit <- function(df, locations, outcome, params, pooled, city) {
  
  output  <- run_model(
    df, locations, outcome,
    fam_type    ='binom',
    pooled      =pooled,
    longterm_df =params$longterm,
    autocorr_lag=params$autocorr,
    arglag_tmp  =params$arglag_tmp,
    argvar_tmp  =params$argvar_tmp,
    arglag_tp   =params$arglag_tp,
    arglag_fg10 =params$arglag_fg10,
  )
  results <- format_results(output, outcome, city=NULL)
  message('Model AIC: ', output$meta_comp$sum_aic)
  list(output=output, results=results)
}

## Heterogeneity statistics from an mvmeta object
get_heterogeneity <- function(mv_fit) {
  s  <- summary(mv_fit)
  q  <- s$qstat$Q[".all"]
  df <- s$qstat$df[".all"]
  p  <- s$qstat$pvalue[".all"]
  list(
    AIC=round(s$AIC, 1),
    Q  =round(q, 1),
    df =df,
    p  =ifelse(p < 0.001, "<0.001", round(p, 3)),
    I2 =round(max(0, (q - df) / q) * 100, 1)
  )
}

## Likelihood-ratio test between a null and an alternative mvmeta model
get_lr_test <- function(mv_null, mv_mod) {
  lrt  <- lrtest(mv_null, mv_mod)
  lr_p <- lrt$`Pr(>Chisq)`[2]
  list(
    stat=lrt$LogLik[2],
    df  =lrt$Df[2],
    p   =ifelse(lr_p < 0.001, "<0.001", round(lr_p, 3))
  )
}

## Build meta-regression comparison table
make_meta_reg_table <- function(meta, predictors_list) {
  
  het_null <- get_heterogeneity(meta$mv)
  
  null_row <- data.frame(
    `Meta-predictors`="Intercept only",
    LR_stat=NA, LR_df=NA, LR_p=NA,
    AIC=het_null$AIC,
    Q_stat=het_null$Q, Q_df=het_null$df, Q_p=het_null$p,
    I2=het_null$I2,
    check.names=FALSE
  )
  
  pred_rows <- imap(predictors_list, function(label, predictor) {
    pred_vals <- sapply(locations, function(loc) {
      vals <- df[[predictor]][df$location == loc]
      if (predictor %in% c('lat', 'lon', "elev")) unique(vals) else mean(vals, na.rm=TRUE)
    })
    
    mv_mod <- tryCatch(
      mvmeta(meta$ymat ~ pred_vals, meta$Slist, method="ml"),
      error=function(e) NULL
    )
    if (is.null(mv_mod)) return(NULL)
    
    het <- get_heterogeneity(mv_mod)
    lr  <- get_lr_test(meta$mv, mv_mod)
    
    data.frame(
      `Meta-predictors`=label,
      LR_stat=lr$stat, LR_df=lr$df, LR_p=lr$p,
      AIC=het$AIC,
      Q_stat=het$Q, Q_df=het$df, Q_p=het$p,
      I2=het$I2,
      check.names=FALSE
    )
  }) |> compact()
  
  do.call(rbind, c(list(null_row), pred_rows)) |> `rownames<-`(NULL)
}

#### MAIN CODE ####

# Load data
load(paste0(path_project, 'Results/ind_final_params.RData'))
load('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData')
df        <- df_dlnm; rm(df_dlnm)

locations <- df |> distinct(location, elev) |> arrange(elev) |> pull(location) |> as.character()

m <- length(locations)
outcome   <- 'total_patients'

# Lag windows used for exposure plots
tmp_lags <- list(
  ALL         =c(0, 3, 4),  Una        =c(0, 3),
  Kangra= c(0, 3),       Nahan      =c(1, 8),
  Hamirpur    =c(1, 7),  Mandi= c(1, 4),
  Chamba      =c(1, 3),  Kullu      =c(1, 3),
  Dharamshala =c(0, 3),  Solan      =c(0, 3),
  'Reckong Peo'= c(0, 3)
)
tp_lags <- list(
  ALL         =c(0, 3, 9), Una        =c(0, 7),
  Kangra=c(0, 7),       Nahan      =c(0, 7),
  Hamirpur    =c(0, 7), Mandi      =c(0, 7),
  Chamba      =c(0, 7), Kullu      =c(0, 7),
  Dharamshala =c(0, 7),  Solan      =c(0, 7),
  'Reckong Peo'= c(0, 7)
)

# Part 1 – Pooled model ----

final_params <- final_params$pooled; gc()
final_params$shape <- 'quasi'

df$tp <- df$tp_roll_log
# df$mx2t <- df$mx2t_roll

res      <- run_bestfit(df, locations, outcome, params=final_params, pooled=TRUE, city=NULL)
output   <- res$output
main_res <- res$results


## Shared helper: build contour plots for all cities + pooled
make_contour_plots <- function(meta, xlab, reference, max_dev=0.20, sig=FALSE) {
  shared_lims <- c(0.999 - max_dev, 1 + max_dev) * 1.001
  shared_brks <- pretty(shared_lims, n=8)
  
  args <- list(meta=meta, xlab=xlab, reference=reference, sig=sig,
               shared_lims=shared_lims, shared_brks=shared_brks)
  
  setNames(
    lapply(c("ALL", locations), function(l)
      do.call(city_con_res, c(args, list(city=l, title=if (l == "ALL") "Pooled" else l)))$p),
    c("ALL", locations)
  )
}

cont_sig <- list(
  tmp=make_contour_plots(output$meta_tmp, "Temperature (%)",   "50%",  sig=TRUE),
  tp =make_contour_plots(output$meta_tp,  "Precipitation (%)", "20%",  sig=TRUE)
)

# Inspect significant lags
wrap_plots(cont_sig$tmp, nrow=3) +
  plot_layout(guides='collect') +
  plot_annotation(
    title='Significant exposure-lag-response',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
  ) & theme(legend.position='bottom')

wrap_plots(cont_sig$tp, nrow=3) +
  plot_layout(guides='collect') +
  plot_annotation(
    title='Significant exposure-lag-response',
    theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
  ) & theme(legend.position='bottom')


## Shared helper: build exposure-response plots for a given meta-object
make_exp_plots <- function(meta, lags_list, col="grey30", xlab, reference) {
  lapply(setNames(c("ALL", locations), c("ALL", locations)), function(l)
    lapply(lags_list[[l]], function(lag)
      city_exp_res(city=l, lag=lag, meta=meta, col=col, xlab=xlab,
                   title=paste0("Lag=", lag), reference=reference)$p_blup
    )
  )
}

exp <- list(
  tmp=make_exp_plots(output$meta_tmp, tmp_lags, xlab="Temperature (%)",   reference="50%"),
  tp =make_exp_plots(output$meta_tp,  tp_lags,  xlab="Precipitation (%)", reference="20%")
)


saveRDS(output, file="2_DLNM/Results/pooled_model.rds")
saveRDS(list(res=main_res, cont=cont_sig, exp=exp), file="2_DLNM/Results/pooled_plots.RData")

# Part 2 – Meta-predictor models ----

predictors_list <- list(
  lat         = 'Latitude',
  lon         = 'Longitude',
  elev        = "Altitude",
  t2m         = "Mean temperature",
  tp          = 'Total precipitation',
  rh          = "Relative humidity",
  fg10        = 'Wind speed',
  pm25        = "Air pollution",
  popu        = "Population density",
  hosp_n      = "Number of hospitals",
  hosp_teach  = "Teaching hospital",
  insurance   = "Health insurance coverage",
  decision    = "Healthcare decision-making",
  meducation  = "Maternal education",
  hsize       = "Household size",
  hwealth     = "Household wealth"
  # mcv1        = "MCV1 vaccination coverage",
  # dtp1        = "DTP1 vaccination coverage",
  # dtp3        = "DTP3 vaccination coverage",
  # anc         = "Antenatal care coverage",
  # malnutrition = "Malnutrition prevalence",
  # pnc         = "Postnatal care coverage",
  # sba         = "Skilled birth attendance"
  # tt          = "Tetanus toxoid coverage",
  )
tbl_tmp <- make_meta_reg_table(meta=output$meta_tmp, predictors_list)
tbl_tp  <- make_meta_reg_table(meta=output$meta_tp,  predictors_list)
print(tbl_tmp); print(tbl_tp)

saveRDS(list(tbl_tmp=tbl_tmp, tbl_tp=tbl_tp), file="2_DLNM/Results/meta_tables.RData")

# Part 3 – Subgroup models ----

load(paste0(path_project, 'Results/ind_final_params.RData'))


outcomes <- c('total_patients',
            'F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'general', 'ophthalmology') #, 'ent', 'gynaecology')
sub_out  <- list()

## Collect all contour z-values first to set shared colour scales
z_vals <- map(set_names(outcomes), function(i) {
  message(i)
  res        <- run_bestfit(df, locations, i, params=final_params$pooled, pooled=TRUE, city=NULL)
  sub_out[[i]] <<- res$output
  
  list(
    tmp=city_con_res(meta=res$output$meta_tmp, xlab="Temperature (%)",
                       reference="50%", sig=FALSE, city="ALL", title=i)$df$z,
    tp =city_con_res(meta=res$output$meta_tp,  xlab="Precipitation (%)",
                       reference="20%", sig=FALSE, city="ALL", title=i)$df$z
  )
})

make_scale <- function(all_z, max_dev=0.20) {
  lims <- c(0.999 - max_dev, 1 + max_dev) * 1.001
  list(shared_lims=lims, shared_brks=pretty(lims, n=8))
}

scale_tmp <- make_scale(unlist(map(z_vals, "tmp")))
scale_tp  <- make_scale(unlist(map(z_vals, "tp")))

## Contour and exposure plots per outcome
pooled_meta    <- list(tmp=list(), tp=list())
cont_outcomes  <- list(tmp=list(), tp=list())
exp_outcomes   <- list(tmp=list(), tp=list())

for (i in outcomes) {
  pooled_meta$tmp[[i]] <- sub_out[[i]]$meta_tmp
  pooled_meta$tp[[i]]  <- sub_out[[i]]$meta_tp
  
  cont_outcomes$tmp[[i]] <- do.call(city_con_res, c(
    list(meta=sub_out[[i]]$meta_tmp, xlab="Temperature (%)", reference="50%", sig=FALSE),
    scale_tmp, list(city="ALL", title=replace_exact_matches(i, nm))
  ))$p
  cont_outcomes$tp[[i]]  <- do.call(city_con_res, c(
    list(meta=sub_out[[i]]$meta_tp,  xlab="Precipitation (%)", reference="20%", sig=FALSE),
    scale_tp,  list(city="ALL", title=replace_exact_matches(i, nm))
  ))$p
  
  exp_outcomes$tmp[[i]] <- lapply(tmp_lags$ALL, function(lag)
    city_exp_res(city="ALL", lag=lag, meta=sub_out[[i]]$meta_tmp,
                 col="grey30", xlab="Temperature (%)",
                 title=paste0("Lag=", lag), reference="50%")$p_blup)
  exp_outcomes$tp[[i]]  <- lapply(tp_lags$ALL, function(lag)
    city_exp_res(city="ALL", lag=lag, meta=sub_out[[i]]$meta_tp,
                 col="grey30", xlab="Precipitation (%)",
                 title=paste0("Lag=", lag), reference="20%")$p_blup)
}

saveRDS(sub_out, file="2_DLNM/Results/subgroup_models.rds")
saveRDS(list(pooled_meta=pooled_meta, cont_outcomes=cont_outcomes, exp_outcomes=exp_outcomes), file="2_DLNM/Results/subgroup_plots.RData")

# Sub-group statistical tests ####

# Option 2: Post-hoc Z-test on pooled BLUP estimates (pragmatic, defensible)
# Since you already have mvmeta pooled estimates (mean vector + covariance), 
# you can compare the pooled coefficients across models using a Z-test or Wald test:

# compare_blup_coefs <- function(meta1, meta2, coef_idx = 1) {
#   # Extract pooled coef and its variance for each model
#   b1 <- coef(meta1$mv)[coef_idx]
#   b2 <- coef(meta2$mv)[coef_idx]
#   
#   v1 <- vcov(meta1$mv)[coef_idx, coef_idx]
#   v2 <- vcov(meta2$mv)[coef_idx, coef_idx]
#   
#   # Z-test assuming independence between models
#   z  <- (b1 - b2) / sqrt(v1 + v2)
#   p  <- 2 * pnorm(abs(z), lower.tail = FALSE)
#   
#   list(b1 = b1, b2 = b2, z = z, p = p)
# }