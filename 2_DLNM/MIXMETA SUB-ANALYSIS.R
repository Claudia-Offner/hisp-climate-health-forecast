
final_params <- final_params$pooled; gc()

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
                            ref_line = 1) {
  
  plot_df <- plot_df %>% mutate(perc = as.numeric(sub("%", "", perc)))

  ggplot(plot_df, aes(x = perc, y = fit,
                      color = .data[[group_col]],
                      fill  = .data[[group_col]])) +
    geom_hline(yintercept = ref_line, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(trans = "log", 
                       breaks = scales::pretty_breaks(n = 6)) +
    labs(x = exposure_label, y = "RR",
         color = group_col, fill = group_col) +
    theme_classic() +
    theme(
      axis.line=element_line(color='gray20'),
      axis.text=element_text(color='gray20'),
      axis.ticks=element_line(color='gray20'),
      plot.title=element_text(face='bold', hjust=0.5, size=12),
      legend.position='bottom'
    )
}


subgroup_specs <- list(
  gender = c("F", "M"),
  age    = c("under_5", "age5_64", "over_65")
)

res_list  <- list()
res_data <- list()

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
plot_risk_curve(res_data$gender$tmp, group_col = "gender", exposure_label = "Temperature (°C)")
plot_risk_curve(res_data$gender$tp, group_col = "gender", exposure_label = "Preciptiation")
plot_risk_curve(res_data$age$tmp, group_col = "age", exposure_label = "Temperature (°C)")
plot_risk_curve(res_data$age$tp, group_col = "age", exposure_label = "Preciptiation")


