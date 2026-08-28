### ------------------------------------------------------------------------ ### 
###  DNLM – Pooled and subgroup models with meta predictors (HISP)
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

## Build meta-regression comparison table
make_meta_reg_table <- function(meta, predictors_list) {
  
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

export <- paste0(a, '_', model)
if(a=='t2m'){ dev <- 0.2 } else{ dev<-0.4 }

# Part 1 - Pooled model ----
final_params <- final_params$pooled; gc()

output  <- run_model(
  df, locations, outcome, pooled=TRUE,
  fam_type    =final_params$shape,
  longterm_df =final_params$longterm,
  autocorr_lag=final_params$autocorr,
  arglag_tmp  =final_params$arglag_tmp,
  argvar_tmp  =final_params$argvar_tmp,
  arglag_tp   =final_params$arglag_tp
)

# Get contour plots
make_contour_plots <- function(meta, xlab, reference, max_dev=dev, sig=FALSE) {
  shared_lims <- c(0.999 - max_dev, 1 + max_dev) * 1.001
  shared_brks <- pretty(shared_lims, n=8)
  
  args <- list(meta=meta, xlab=xlab, reference=reference, sig=sig,
               shared_lims=shared_lims, shared_brks=shared_brks)
  
  setNames(
    lapply(c("ALL", locations), function(l)
      do.call(city_con_res, c(args, list(city=l, title=if (l == "ALL") "Contour" else l)))$p),
    c("ALL", locations)
  )
}
cont_sig <- list(
  tmp=make_contour_plots(meta=output$meta_tmp, "Temperature (%)",   "50%",  sig=TRUE),
  tp =make_contour_plots(meta=output$meta_tp,  "Precipitation (%)", tp_ref,  sig=TRUE)
)

# Identify relevant lags
# wrap_plots(cont_sig$tmp, nrow=3) +
#   plot_layout(guides='collect') +
#   plot_annotation(
#     title='Significant exposure-lag-response',
#     theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
#   ) & theme(legend.position='bottom')
# 
# wrap_plots(cont_sig$tp, nrow=3) +
#   plot_layout(guides='collect') +
#   plot_annotation(
#     title='Significant exposure-lag-response',
#     theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))
#   ) & theme(legend.position='bottom')

if(any(grepl('immediate', export))){
  tmp_lags <- list(ALL = c(0, 4, 10))
  tp_lags <- list(ALL  = c(0, 4, 7))
} else if (any(grepl('cumulative', export))){
  tmp_lags <- list(ALL = c(0, 2, 6))
  tp_lags <- list(ALL  = c(0, 4, 10))
}

# Exposure-response plots
make_exp_plots <- function(meta, lags_list, col="grey30", xlab, reference) {
  lapply(setNames(c("ALL", locations), c("ALL", locations)), function(l) {
    
    # per-lag exposure-response curves
    lag_plots <- lapply(lags_list[[l]], function(lag)
      city_exp_res(city=l, lag=lag, meta=meta, col="grey30", xlab=xlab,
                   title=paste0("Lag=", lag), reference=reference)$p_blup
    )
    
    # overall (cumulative) exposure-response curve
    overall_plot <- city_exp_res(city=l, lag=NULL, meta=meta, col=col, xlab=xlab,
                                 title="Cumulative", reference=reference)$p_blup
    
    c(lag_plots, list(Overall = overall_plot))
  })
}

exp <- list(
  tmp=make_exp_plots(output$meta_tmp, tmp_lags, col='#8968CD', xlab="Temperature (%)",   reference="50%"),
  tp =make_exp_plots(output$meta_tp,  tp_lags,  col='#8968CD', xlab="Precipitation (%)", reference=tp_ref)
)

# EXPORT
saveRDS(output, file=paste0('2_DLNM/Results/', export, '_pooled_model.rds'))
saveRDS(list(cont=cont_sig, exp=exp), 
        file=paste0('2_DLNM/Results/', export, '_pooled_plots.RData'))


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
  hosp_teach  = "Teaching hospital",
  insurance   = "Health insurance coverage",
  decision    = "Healthcare decision-making",
  meducation  = "Maternal education",
  hsize       = "Household size",
  hwealth     = "Household wealth"
  # hosp_n      = "Number of hospitals",
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

saveRDS(list(tbl_tmp=tbl_tmp, tbl_tp=tbl_tp), 
        file=paste0('2_DLNM/Results/', export, '_meta_tables.RData'))

# Part 3 – Subgroup models ----

outcomes <- c('total_patients',
              'F','M','under_5','age5_14','age15_24','age25_64','over_65',
              'general', 'ophthalmology') #, 'ent', 'gynaecology')
sub_out  <- list()

# Run models and get z-values
z_vals <- map(set_names(outcomes), function(i) {
  message(i)
  output  <- run_model(
    df, locations, i, pooled=TRUE,
    fam_type    =final_params$shape,
    longterm_df =final_params$longterm,
    autocorr_lag=final_params$autocorr,
    arglag_tmp  =final_params$arglag_tmp,
    argvar_tmp  =final_params$argvar_tmp,
    arglag_tp   =final_params$arglag_tp
  )

  sub_out[[i]] <<- output

  list(
    tmp=city_con_res(meta=output$meta_tmp, xlab="Temperature (%)",
                     reference="50%", sig=FALSE, city="ALL", title=i)$df$z,
    tp =city_con_res(meta=output$meta_tp,  xlab="Precipitation (%)",
                     reference=tp_ref, sig=FALSE, city="ALL", title=i)$df$z
  )
})
make_scale <- function(all_z, max_dev=dev) {
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

  print(i)

  pooled_meta$tmp[[i]] <- sub_out[[i]]$meta_tmp
  pooled_meta$tp[[i]]  <- sub_out[[i]]$meta_tp

  cont_outcomes$tmp[[i]] <- do.call(city_con_res, c(
    list(meta=sub_out[[i]]$meta_tmp, xlab="Temperature (%)", reference="50%", sig=FALSE),
    scale_tmp, list(city="ALL", title=replace_exact_matches(i, nm))
  ))$p
  cont_outcomes$tp[[i]]  <- do.call(city_con_res, c(
    list(meta=sub_out[[i]]$meta_tp,  xlab="Precipitation (%)", reference=tp_ref, sig=FALSE),
    scale_tp,  list(city="ALL", title=replace_exact_matches(i, nm))
  ))$p

  exp_outcomes$tmp[[i]] <- city_exp_res(city="ALL", lag=NULL, meta=sub_out[[i]]$meta_tmp,
                                        col="grey30", xlab="Temperature (%)",
                                        title=replace_exact_matches(i, nm), reference="50%")$p_blup
  exp_outcomes$tp[[i]]  <- city_exp_res(city="ALL", lag=NULL, meta=sub_out[[i]]$meta_tp,
                                        col="grey30", xlab="Precipitation (%)",
                                        title=replace_exact_matches(i, nm), reference=tp_ref)$p_blup
}

saveRDS(sub_out, file=paste0('2_DLNM/Results/', export, '_subgroup_models.rds'))
saveRDS(list(pooled_meta=pooled_meta, cont_outcomes=cont_outcomes, exp_outcomes=exp_outcomes),
        file=paste0('2_DLNM/Results/', export, '_subgroup_plots.RData'))
