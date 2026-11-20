### ------------------------------------------------------------------------ ### 
###  DNLM – Pooled and subgroup models with meta predictors
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Build meta-regression comparison table
make_meta_reg_table <- function(meta, predictors_list) {

  ## Likelihood-ratio test between a null and an alternative mvmeta model
  get_lr_test <- function(mv_null, mv_mod) {
    lrt  <- lrtest(mv_null, mv_mod)
    lr_p <- lrt$`Pr(>Chisq)`[2]
    list(
      stat=lrt$Chisq[2],
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
    pred_vals <- sapply(rownames(meta$coefs), function(loc) {
      vals <- df[[predictor]][df$location == loc]
      if (predictor %in% c('lat', 'lon', "elev")) unique(vals) else mean(vals, na.rm=TRUE)
    })

    mv_mod <- tryCatch(
      mixmeta(meta$coefs ~ pred_vals, S=meta$vcovs, method="ml"),
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
if(a=='t2m'){ dev <- 0.3} else{ dev <- 0.4 }
if(model=='cumulative'){ ttl <- 'Seven-day'} else{ ttl <- 'Daily' }

# 1. Pooled model ----

final_params <- final_params$pooled; gc()
df$groups <- df$location
output  <- run_model(
  df, groups = unique(df$groups), outcome = "total_patients",
  group_meta = NULL, subgroup_var = NULL,
  fam_type    =final_params$shape,
  lag_window  = c(0, final_params$lag_window),
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
      do.call(city_con_res, c(args, list(city=l, title=if (l == "ALL") paste0(ttl, " risk summary") else l)))$p),
    c("ALL", locations)
  )
}
cont_sig <- list(
  tmp=make_contour_plots(meta=output$meta_tmp, "Temperature percentiles (%)",   "50%",  sig=TRUE),
  tp =make_contour_plots(meta=output$meta_tp,  "Precipitation percentiles (%)", tp_ref,  sig=TRUE)
)

# ggsave(paste0(path_report, 'test.png'), cont_sig$tmp$ALL, width = 5, height = 5, dpi = 300, bg = 'white')
# ggsave(paste0(path_report, 'test2.png'), cont_sig$tp$ALL, width = 5, height = 5, dpi = 300, bg = 'white')

# Get exposure curves
exp <- list(
  tmp=city_exp_res(meta=output$meta_tmp, col='#8968CD', xlab="Temperature percentiles (%)", title=paste0(ttl, " cumulative risk"), reference="50%")$p_blup,
  tp=city_exp_res(meta=output$meta_tp, col='#8968CD', xlab="Precipitation percentiles (%)", title=paste0(ttl, " cumulative risk"), reference=tp_ref)$p_blup
)

# Get lag curves
lag <- list(
  tmp=city_lag_res(meta=output$meta_tmp, percentiles=c(5, 95), col='grey30', title=paste0(ttl, ' lagged risk at')),
  tp=city_lag_res(meta=output$meta_tp, percentiles=c(95), col='grey30', title=paste0(ttl, ' lagged risk at'))
)

# EXPORT
saveRDS(output, file=paste0('outputs/', export, '_pooled_model.rds'))
saveRDS(list(cont=cont_sig, exp=exp, lag=lag), file=paste0('outputs/', export, '_pooled_plots.RData'))

# 2. Estimate attributable numbers from main model ####

if(a=='t2m'){

  # Function to get location and weather (i.e. scenario) specific attributable numbers (blup)
  attribution_multi <- function(output, scenarios, nsim = 2000) {

    by_exp <- split(names(scenarios), vapply(scenarios, `[[`, character(1), "exp"))
    summary_rows <- list()
    sims <- list()

    # Extract exposure-specific BLUP-predictions
    for (exp in names(by_exp)) {
      if (exp == 'tmp') {
        list2env(output$meta_tmp, envir = environment())
        argvar <- final_params$argvar_tmp; arglag <- final_params$arglag_tmp
      } else if (exp == 'tp') {
        list2env(output$meta_tp, envir = environment())
        argvar <- final_params$argvar_tp; arglag <- final_params$arglag_tp
      }
      argvar <- if (argvar == 1) list(fun="lin") else list(fun="ns", df=argvar)
      blup_all  <- blup(mv, vcov = TRUE)
      locations <- names(output$models)

      # Set location specific cross basis information
      for (loc in locations) {
        sub <- df %>% filter(location == loc) %>% mutate(total_patients = replace_na(total_patients, 1))
        bl  <- blup_all[[loc]]
        tot_patients <- sum(sub$total_patients, na.rm = TRUE)
        x_range <- range(sub[[exp]], na.rm = TRUE)
        b_mean <- crossbasis(sub[[exp]], lag = c(0, 10),
                             argvar = argvar,
                             arglag = list(fun = 'ns', df = arglag))
        # pred <- crosspred(b_mean, coef = bl$blup, vcov = bl$vcov, cen = cen['50%', loc])
        # in_range <- pred$predvar >=  cen['25%', loc] & pred$predvar <=  cen['75%', loc]
        # mmi <- pred$predvar[in_range][which.min(pred$allfit[in_range])]

        # Calculate estimations within relevant scenario range
        for (sn in by_exp[[exp]]) {
          sc <- scenarios[[sn]]
          ref_val <- cen[sc$reference, loc]
          # ref_val <- if(exp=='tmp') mmi else{ cen[sc$reference, loc]}
          rng <- NULL
          if (!is.null(sc$percentile)) {
            key <- paste0(sc$percentile, "%")
            thresh <- if (key %in% rownames(cen)) cen[key, loc] else
              quantile(sub[[exp]], probs = sc$percentile / 100, na.rm = TRUE)
            rng <- switch(sc$direction,
                          above = c(thresh, x_range[2]),
                          below = c(x_range[1], thresh))
          }

          an_t <- attrdl(sub[[exp]], b_mean, sub$total_patients,
                         coef = bl$blup, vcov = bl$vcov, type = "an",
                         cen = ref_val, dir = "forw", tot = FALSE, range = rng)
          # an_pos_est <- sum(an_t[an_t > 0], na.rm=TRUE)
          an_pos_est <- sum(an_t, na.rm=TRUE)


          an_sim <- attrdl(sub[[exp]], b_mean, sub$total_patients,
                           coef = bl$blup, vcov = bl$vcov, type = "an",
                           cen = ref_val, dir = "forw", tot = TRUE,
                           sim = TRUE, nsim = nsim, range = rng)

          summary_rows[[paste(sn, loc, sep = "_")]] <- data.frame(
            scenario = sn, location = loc, total_patients = tot_patients,
            an = an_pos_est
          )
          sims[[sn]][[loc]] <- an_sim
        }
      }

      res = do.call(rbind, summary_rows)

    }

    list(res = do.call(rbind, summary_rows), sims = sims)
  }
  # Function to calculate attributable fractions, CI's and weather total estimates
  build_panel <- function(sn='cold', sim_by_loc=sims$cold, ci = c(0.025, 0.975)) {
    d <- res_all %>% filter(scenario == sn, location %in% locations)

    # Get location specific an and af CIs
    loc_rows <- d %>%
      rowwise() %>%
      mutate(
        an_low  = quantile(sim_by_loc[[location]], ci[1], na.rm = TRUE),
        an_high = quantile(sim_by_loc[[location]], ci[2], na.rm = TRUE),
        af = an / total_patients,
        af_low  = quantile(sim_by_loc[[location]], ci[1], na.rm = TRUE) / total_patients,
        af_high = quantile(sim_by_loc[[location]], ci[2], na.rm = TRUE) / total_patients
      ) %>%
      ungroup() %>%
      dplyr::select(location, an, an_low, an_high, af, af_low, af_high)

    # Calculate total an and af
    overall_sim <- Reduce(`+`, sim_by_loc[locations])
    overall_tot <- sum(d$total_patients)
    overall_an  <- sum(d$an)
    q <- quantile(overall_sim, probs = ci, na.rm = TRUE)

    bind_rows(loc_rows, data.frame(location = "TOTAL",
                                   an = overall_an, an_low = q[1], an_high = q[2],
                                   af = overall_an / overall_tot, af_low = q[1] / overall_tot, af_high = q[2] / overall_tot))
  }

  # Define weather exposure attribution reference and range
  scenarios <- list(
    heat   = list(exp = "tmp", reference = "50%", percentile = 75, direction = "above"),
    cold   = list(exp = "tmp", reference = "50%", percentile = 25, direction = "below"),
    precip = list(exp = "tp",  reference = "0%", percentile = 75, direction = "above")
  )

  out     <- attribution_multi(output, scenarios, nsim = 2000)
  res_all <- out$res
  sims    <- out$sims
  rownames(res_all) <- NULL

  # Get CI's, af and total
  t <- if(model=='cumulative'){'Seven-day'}else{'Daily'}
  attr_df <- bind_rows(
    build_panel("cold",   sims$cold)   %>% mutate(panel = paste0(t, " Cold")),
    build_panel("heat",   sims$heat)   %>% mutate(panel = paste0(t, " Heat")),
    build_panel("precip", sims$precip) %>% mutate(panel = paste0(t, " Precipitation"))
  ) %>%
    mutate(
      panel = factor(panel),
      location = factor(location, levels = c(setdiff(locations, "TOTAL"), "TOTAL"))
    )

  saveRDS(attr_df, file=paste0('outputs/', export, '_attr.RData'))
}

# 3. Meta-predictor models ----

if(a %in% c('t2m', 'mn2t', 'mx2t')){

  predictors_list <- list(
    lat         = 'Latitude',
    lon         = 'Longitude',
    elev        = "Elevation",
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
  )

  tbl_tmp <- make_meta_reg_table(meta=output$meta_tmp, predictors_list)
  tbl_tp  <- make_meta_reg_table(meta=output$meta_tp,  predictors_list)
  print(tbl_tmp); print(tbl_tp)

  saveRDS(list(tbl_tmp=tbl_tmp, tbl_tp=tbl_tp),
          file=paste0('outputs/', export, '_meta_tables.RData'))

}

# 4. Subgroup models ----

if(a=='t2m'){

  make_subgroup_long <- function(df, value_cols, subgroup_var, locations) {

    x_full <- df %>%
      dplyr::select(-total_patients) %>%
      pivot_longer(cols      = dplyr::all_of(value_cols),
                   names_to  = "group",
                   values_to = "total_patients") %>%
      mutate(groups = paste(location, group, sep = "_"))

    group_meta <- x_full %>%
      dplyr::distinct(groups, group, location) %>%
      dplyr::rename(!!subgroup_var := group)

    x <- x_full %>%
      dplyr::select(groups, 2:13, total_patients, tmp, tp, fg10, rh, pm25)

    list(data = x, group_meta = group_meta)
  }

  extract_curves <- function(meta, reference) {

    list2env(meta, envir = environment())

    at_vals <- rowMeans(cen)
    names(at_vals) <- row.names(cen)
    at_vals <- at_vals[!duplicated(at_vals)]
    cp_list <- setNames(vector("list", length(levels)), levels)

    for (lv in levels) {
      term <- paste0("sg", lv)

      idx_coef <- grep(term, names(coef(mv)), fixed = TRUE)
      idx_vcov <- grep(term, rownames(vcov(mv)), fixed = TRUE)

      cp_list[[lv]] <- crosspred(b_mean,
                                 coef = coef(mv)[idx_coef],
                                 vcov = vcov(mv)[idx_vcov, idx_vcov],
                                 model.link = "log",
                                 at  = at_vals,
                                 cen = rowMeans(cen)[reference])
    }

    plot_df <- do.call(rbind, lapply(names(cp_list), function(lv) {
      cp <- cp_list[[lv]]
      x <- data.frame(perc  = names(at_vals), exp = cp$predvar, group = lv,
                      fit = cp$allRRfit, lo = cp$allRRlow, hi = cp$allRRhigh)
      row.names(x) <- NULL
      x
    }))
    names(plot_df)[names(plot_df) == "group"] <- var

    return(plot_df)
  }

  build_lrt_table <- function(res_list) {
    tidy_lrt <- function(lrt_obj, subgroup_var, exposure) {
      d <- as.data.frame(lrt_obj)

      aic_col   <- grep("AIC", names(d), value = TRUE)[1]
      chisq_col <- grep("Chisq$|Chi.?sq", names(d), value = TRUE, ignore.case = TRUE)[1]
      df_col    <- grep("^Df$|^df$", names(d), value = TRUE)[1]
      p_col     <- grep("Pr\\(|p.?val", names(d), value = TRUE, ignore.case = TRUE)[1]

      data.frame(
        subgroup_var = subgroup_var,
        exposure     = exposure,
        model        = rownames(d),
        AIC          = if (!is.na(aic_col))   d[[aic_col]]   else NA,
        Chisq        = if (!is.na(chisq_col)) d[[chisq_col]] else NA,
        Df           = if (!is.na(df_col))    d[[df_col]]    else NA,
        p_value      = if (!is.na(p_col))     d[[p_col]]     else NA,
        row.names = NULL
      )
    }
    do.call(rbind, lapply(names(res_list), function(sv) {
      res <- res_list[[sv]]
      rbind(
        tidy_lrt(res$meta_tmp$lrt, subgroup_var = sv, exposure = "tmp"),
        tidy_lrt(res$meta_tp$lrt,  subgroup_var = sv, exposure = "tp")
      )
    }))
  }

  plot_risk_curve <- function(plot_df, group_col, exposure_label = "Temperature",
                              pal = NULL, labels = NULL, level_order = NULL,
                              xlim = NULL, ylim = NULL, title = NULL,
                              stats = NULL) {

    # Default palettes matching your old gender_pal / age_pal, keyed by group_col
    default_palettes <- list(
      gender = setNames(mypalette[c(2, 1)], c("F", "M")),
      age    = setNames(mypalette[c(3, 6, 7)], c("under_5", "age5_64", "over_65"))
    )

    # Determine level order: explicit arg > names(pal) if supplied > default palette order > data order
    if (is.null(level_order)) {
      if (!is.null(pal) && !is.null(names(pal))) {
        level_order <- names(pal)
      } else if (group_col %in% names(default_palettes)) {
        level_order <- names(default_palettes[[group_col]])
      }
    }

    # Match old make_subgroup_plot() input shape: pred/RR/low/high/Subgroup
    plot_df <- plot_df %>%
      mutate(perc = as.numeric(sub("%", "", perc))) %>%
      rename(pred = perc, RR = fit, low = lo, high = hi) %>%
      rename(Subgroup = !!group_col)

    # Apply factor levels, falling back to whatever's present if level_order doesn't fully cover it
    if (!is.null(level_order) && all(unique(plot_df$Subgroup) %in% level_order)) {
      plot_df <- plot_df %>% mutate(Subgroup = factor(Subgroup, levels = level_order))
    } else {
      plot_df <- plot_df %>% mutate(Subgroup = as.factor(Subgroup))
    }

    if (is.null(pal)) {
      if (group_col %in% names(default_palettes) &&
          all(levels(plot_df$Subgroup) %in% names(default_palettes[[group_col]]))) {
        pal <- default_palettes[[group_col]]
      } else {
        pal <- scales::hue_pal()(nlevels(plot_df$Subgroup))
        names(pal) <- levels(plot_df$Subgroup)
      }
    }

    # Default legend-label dictionaries, keyed by group_col, mapping raw code -> display label
    if (is.null(labels)) {
      default_labels <- list(
        gender = c(F = "Women", M = "Men"),
        age    = c(under_5 = "Under 5", age5_64 = "Ages 5-64", over_65 = "Ages 65+")
      )
      if (group_col %in% names(default_labels) &&
          all(levels(plot_df$Subgroup) %in% names(default_labels[[group_col]]))) {
        labels <- default_labels[[group_col]]
      } else {
        labels <- setNames(levels(plot_df$Subgroup), levels(plot_df$Subgroup))
      }
    }
    labels <- labels[levels(plot_df$Subgroup)]

    p <- ggplot(plot_df, aes(x = pred, y = RR, colour = Subgroup, fill = Subgroup, group = Subgroup)) +
      geom_hline(yintercept = 1, color = "gray40") +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1, colour = NA) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(name = NULL, values = pal, labels = labels) +
      scale_fill_manual(name = NULL, values = pal, labels = labels) +
      labs(title = title, x = exposure_label, y = "RR") +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 9))

    if (!is.null(stats)) {
      p_fmt <- if (stats$p_value < 0.001) "<0.001" else formatC(stats$p_value, format = "f", digits = 3)
      stats_label <- bquote(chi^2 == .(round(stats$Chisq, 2)) * "," ~ df == .(stats$Df) * "," ~ italic(p) == .(p_fmt))
      p <- p + labs(subtitle = stats_label)
    }

    if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim, breaks = pretty(xlim, 5))
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, xlim = range(plot_df$pred, na.rm = TRUE)) +
      scale_y_continuous(breaks = pretty(ylim, 5))

    p
  }

  subgroup_specs <- list(
    gender = c("F", "M"),
    age    = c("under_5", "age5_64", "over_65")
  )

  res_list  <- list()
  res_data <- list()
  res_plots <- list()

  for (sv in names(subgroup_specs)) {

    prep <- make_subgroup_long(df, subgroup_specs[[sv]], sv, locations)

    res_list[[sv]] <- run_model(prep$data, groups=prep$group_meta$groups, outcome = "total_patients",
                                group_meta = prep$group_meta, subgroup_var = sv,
                                fam_type    =final_params$shape,
                                longterm_df =final_params$longterm,
                                autocorr_lag=0,
                                arglag_tmp  =final_params$arglag_tmp,
                                argvar_tmp  =final_params$argvar_tmp,
                                arglag_tp   =final_params$arglag_tp)

    res_data[[sv]]$tmp <- extract_curves(meta=res_list[[sv]]$meta_tmp, reference = "50%")
    res_data[[sv]]$tp <- extract_curves(meta=res_list[[sv]]$meta_tp, reference = tp_ref)
  }

  (lrt_table <- build_lrt_table(res_list) %>% filter(model==2))

  # Age, temperature curve
  res_plots[[paste0(model, '_tmp_gen')]] <- plot_risk_curve(res_data$gender$tmp, group_col = "gender", exposure_label = "Temperature (%)",
                                                            stats = lrt_table %>% filter(subgroup_var == "gender", exposure == "tmp"))
  res_plots[[paste0(model, '_tp_gen')]] <- plot_risk_curve(res_data$gender$tp, group_col = "gender", exposure_label = "Preciptiation (%)", xlim =c(as.numeric(sub("%", "", tp_ref)), 100),
                                                           stats = lrt_table %>% filter(subgroup_var == "gender", exposure == "tp"))
  res_plots[[paste0(model, '_tmp_age')]] <- plot_risk_curve(res_data$age$tmp, group_col = "age", exposure_label = "Temperature (%)",
                                                            stats = lrt_table %>% filter(subgroup_var == "age", exposure == "tmp"))
  res_plots[[paste0(model, '_tp_age')]]  <- plot_risk_curve(res_data$age$tp, group_col = "age", exposure_label = "Preciptiation (%)", xlim =c(as.numeric(sub("%", "", tp_ref)), 100),
                                                           stats = lrt_table %>% filter(subgroup_var == "age", exposure == "tp"))

  saveRDS(list(res_data=res_data, lrt_table=lrt_table, res_plots=res_plots),
          file=paste0('outputs/', export, '_subgroup_plots.RData'))

}

