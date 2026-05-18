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
    labs(title=title) +
    theme_void() + 
    theme(plot.title=element_text(hjust=0.5, size=size, face="bold"))
  
}

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

# Model 1 - Pooled results ----

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


# Export ####
# saveRDS(output, file = "2_DLNM/pooled_model.rds")
# saveRDS(list(cont, exp), file = "2_DLNM/pooled_plots.RData")

# Meta-regression ####

# META REGRESSION: Elevation and population density
# meta_reg <- function(meta, xlab='Exposure (%)', predictor='elev', perc=c(0.1, 0.9),
#                      col='darkred'){
# 
#   # Get unique predictors
#   pred <- sapply(locations, function(loc) {
#     if(predictor != 'elev'){
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
# meta_list <- list('meta_tmp'=list(xlab='Temperature (%)', col='#D74B3C'),
#                   'meta_tp'=list(xlab='Precipitation (%)', col='#1D98EB'))
# elev_plots <- list()
# popu_plots <- list()
# rh_plots <- list()
# for(v in names(meta_list)){
#   elev_plots[[v]] <- meta_reg(meta=output[[v]], xlab=meta_list[[v]]$xlab, predictor='elev', col=meta_list[[v]]$col)
#   popu_plots[[v]] <- meta_reg(meta=output[[v]], xlab=meta_list[[v]]$xlab, predictor='popu', col=meta_list[[v]]$col)
#   rh_plots[[v]] <- meta_reg(meta=output[[v]], xlab=meta_list[[v]]$xlab, predictor='rh', col=meta_list[[v]]$col)
# }
# 
# (elev_res <- wrap_plots(elev_plots, nrow=1)+plot_annotation(title='Exposure-resonse by elevation percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# (popu_res <- wrap_plots(popu_plots, nrow=1)+plot_annotation(title='Exposure-resonse by population density percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# (rh_res <- wrap_plots(rh_plots, nrow=1)+plot_annotation(title='Exposure-resonse by humidity percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# 
library(lmtest)
library(flextable)
library(officer)

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

prep_display <- function(df) {
  data.frame(
    `Meta-predictors` = df[["Meta-predictors"]],
    `Stat`   = fmt_num(df$LR_stat),
    `df`     = fmt_int(df$LR_df),
    `p`      = fmt_p(df$LR_p),
    `AIC`    = fmt_num(df$AIC),
    `Stat `  = fmt_num(df$Q_stat),
    `df `    = fmt_int(df$Q_df),
    `p `     = fmt_p(df$Q_p),
    `I2 (%)` = fmt_num(df$I2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
fmt_p   <- function(x) ifelse(is.na(x), "\u2013", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_num <- function(x, d = 1) ifelse(is.na(x), "\u2013", formatC(x, format = "f", digits = d))
fmt_int <- function(x) ifelse(is.na(x), "\u2013", as.character(as.integer(x)))


section_row <- function() {
  data.frame(
    `Meta-predictors` = "",
    `Stat`   = "", `df` = "", `p` = "",
    `AIC`    = "",
    `Stat `  = "", `df ` = "", `p ` = "",
    `I2 (%)` = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Define predictors
predictors_list <- list(
  elev = "Altitude",
  popu = "Population density" #,
)
(tbl_tmp <- make_meta_reg_table(meta=output$meta_tmp, predictors_list))
(tbl_tp  <- make_meta_reg_table(output$meta_tp,  predictors_list))

# Format flextable
d_tmp <- prep_display(tbl_tmp)
d_tp  <- prep_display(tbl_tp)

combined <- rbind(section_row(), d_tmp, section_row(), d_tp)
tmp_rows   <- 2:4
tp_rows    <- 6:9
label_rows <- c(1, 5)
intercept_rows <- c(2, 6)   # italic "Intercept only" rows

# Build flextable
ft <- flextable(combined) |>
  
  # Spanning header
  add_header_row(
    values    = c("", "LR test", "", "", "Model fits", "Cochran Q test", "", "", "I\u00b2"),
    colwidths = c(1, 1, 1, 1, 1, 1, 1, 1, 1)
  ) |>
  merge_at(i = 1, j = 2:4, part = "header") |>
  merge_at(i = 1, j = 6:8, part = "header") |>
  set_header_labels(
    `Meta-predictors` = "Meta-predictors",
    `Stat`   = "Stat", `df` = "df", `p` = "p",
    `AIC`    = "AIC",
    `Stat `  = "Stat", `df ` = "df", `p ` = "p",
    `I2 (%)` = "(%)"
  ) |>
  
  # Section label text
  compose(i = 1, j = 1, part = "body",
          value = as_paragraph(as_chunk("Temperature", props = fp_text(bold = TRUE, font.family = "Times New Roman", font.size = 12)))) |>
  compose(i = 5, j = 1, part = "body",
          value = as_paragraph(as_chunk("Precipitation", props = fp_text(bold = TRUE, font.family = "Times New Roman", font.size = 12)))) |>
    compose(
    i = 1, j = 2, part = "header",
    value = as_paragraph(
      as_chunk("LR test", props = fp_text(font.family = "Times New Roman", font.size = 12)),
      as_chunk("a",       props = fp_text(font.family = "Times New Roman", font.size = 12,
                                          vertical.align = "superscript", italic=TRUE))
    )
  ) |>
  compose(
    i = 1, j = 6, part = "header",
    value = as_paragraph(
      as_chunk("Cochran Q test", props = fp_text(font.family = "Times New Roman", font.size = 12)),
      as_chunk("b",              props = fp_text(font.family = "Times New Roman", font.size = 12,
                                                 vertical.align = "superscript", italic=TRUE))
    )
  ) |>
  # Merge section label rows across all columns
  merge_h(i = label_rows, part = "body") |>
  
  # Italic for "Intercept only" rows
  italic(i = intercept_rows, part = "body") |>
  
  # Alignment
  align(j = 1,   align = "left",   part = "all") |>
  align(j = 2:9, align = "right",  part = "body") |>
  align(j = 2:9, align = "center", part = "header") |>
  align(i = label_rows, j = 1, align = "left", part = "body") |>
  
  # Font & size
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 12, part = "all") |>
  
  # Column widths
  width(j = 1,   width = 1.55) |>
  width(j = 2:4, width = 0.50) |>
  width(j = 5,   width = 0.55) |>
  width(j = 6:8, width = 0.50) |>
  width(j = 9,   width = 0.45) |>
  
  # Borders
  border_remove() |>
  hline_top(part = "header", border = fp_border(width = 1.5)) |>
  hline(i = 1, part = "header", border = fp_border(width = 0.75)) |>
  hline(i = 2, part = "header", border = fp_border(width = 0.75)) |>
  hline_bottom(part = "body",   border = fp_border(width = 1.5)) |>
  # Thin rule between the two sections
  hline(i = 4, part = "body",   border = fp_border(width = 0.75)) |>
  
  padding(padding.top = 2, padding.bottom = 2, part = "all") |>
  set_table_properties(layout = "fixed")

# ── Save ──────────────────────────────────────────────────────────────────────
ft


# Convert flextable to LaTeX
tab <- flextable:::knit_to_latex(
  ft,
  bookdown = FALSE,
  quarto = FALSE
)

# Optional footnote text
footnote_text <- paste(
  "All models were based on the relative scale measurement of temperature with ",
  "reference at 50\\textsuperscript{th} percentile for temperature and ",
  "20\\textsuperscript{th} percentile for precipitation.",
  "AIC: Akaike Information Criterion; $I^2$: percentage of residual heterogeneity. \\\\[4pt]",
  
  "\\textsuperscript{$a$} Likelihood ratio (LR) tests were used to test the significance",
  "of meta-predictors using the intercept-only model as reference. \\\\[2pt]",
  
  "\\textsuperscript{$b$} Cochran's Q test was used to test the significance of residual",
  "heterogeneity, where no heterogeneity is the null hypothesis. \\\\[2pt]"
)

# Wrap in a full LaTeX table environment with caption, label, and footnote
latex_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\fontsize{12}{14}\\selectfont",   # 12 pt font with 14 pt line spacing
  "\\caption{Meta-regression results for temperature and precipitation}",
  "\\label{tab:MT_meta_reg}",
  "",
  tab,
  "",
  "\\vspace{2mm}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  paste0("\\textit{Note:} ", footnote_text),
  "\\end{minipage}",
  "\\end{table}"
)

# Write to file
writeLines(
  latex_table,
  paste0(path_report, "II. dlnm/MT_meta_reg.tex")
)

# MF2: Pooled exposure-response ####

total_er_overall <- main_res$p_er_overall & theme(legend.position = "none", plot.title = element_blank())

legend_plot <- ggplot() +
  geom_line(aes(x=0, y=0, linetype="City-specific (BLUP)", colour="City-specific (BLUP)"), linewidth=1.2) +
  geom_line(aes(x=0, y=0, linetype="Pooled", colour="Pooled"), linewidth=1.2) +
  scale_linetype_manual(name="", values = c("City-specific (BLUP)" = "dashed", "Pooled" = "solid")) +
  scale_colour_manual(name="",   values = c("City-specific (BLUP)" = "grey70", "Pooled" = "black")) +
  guides(colour = guide_legend(override.aes = list(linewidth = c(1, 1.5)))) +
  theme_void() +
  theme(legend.position  = "top",
        legend.title     = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        legend.text      = element_text(size = 10))

fig <- plot_grid(total_er_overall, legend_plot, ncol=1, rel_heights=c(1, 0.05))
ggsave(paste0(path_report, 'II. dlnm/MF2_pool_exp.png'), fig, width=8, height=4, dpi=300, bg='white')

# MF3: Pooled exposure-lag-response #####

make_figure_overall <- function(contours, exposure) {
  
  strip_legend <- function(p) p + theme(legend.position="none")
  strip_title  <- function(p) p + theme(plot.title=element_blank())
  title_as_annotation <- function(p) {
    ttl <- ggplot_build(p)$plot$labels$title
    p + theme(plot.title=element_blank(),
              plot.tag.position=c(0.6, 0.95),
              plot.tag=element_text(size=12, face='bold', hjust=0.5, vjust=0.5)) +
      labs(tag=ttl)
  }
  row_title <- function(label) {
    ggplot() + 
      labs(title=label) +
      theme_void() +
      theme(plot.title=element_text(size=13, face="bold", hjust=-0.05))
  }
  
  # Strip legends and fix titles
  con_tmp  <- strip_title(strip_legend(contours$tmp$ALL))
  con_tp   <- strip_title(strip_legend(contours$tp$ALL))
  lag_tmp  <- lapply(exposure$tmp$ALL, function(p) title_as_annotation(strip_legend(p)))
  lag_tp   <- lapply(exposure$tp$ALL,  function(p) title_as_annotation(strip_legend(p)))
  
  # Extract legends 
  con_legend <- get_legend(contours$tmp$ALL + 
                             guides(fill = guide_colorbar(
                               barwidth  = unit(3.5, "cm"),
                               barheight = unit(0.3, "cm"),
                               frame.colour = "black",
                               ticks.colour = "black"
                             )) +
                             theme(legend.position  = "bottom",
                                   legend.text      = element_text(size = 7),
                                   legend.title     = element_text(size = 8, face='bold')))
  exp_legend <- get_legend(exposure$tmp$ALL[[1]] + theme(legend.position="bottom"))
  legends <- wrap_plots(list(plot_spacer(), con_legend, exp_legend), nrow=1, widths=c(0.15,1,2))
  
  # Build
  col_tmp <- wrap_plots(c(list(con_tmp), lag_tmp), nrow=1)
  col_tp  <- wrap_plots(c(list(con_tp),  lag_tp),  nrow=1)
  
  all <- wrap_plots(
    list(row_title("(A)"), col_tmp,
         row_title("(B)"), col_tp
    ),
    nrow=2,
    widths=c(0.05, 1))
  
  wrap_plots(list(all, legends), ncol=1, heights=c(1, 0.1))
  
}

fig <- make_figure_overall(cont, exp)
ggsave(filename=paste0(path_report, 'II. dlnm/MF3_pool_exp-lag.png'),  plot=fig, width=9, height=7,  dpi=300, bg='white')

# MF4: Site specific exposure-lag-response ####

make_figure_sites <- function(contours, exposure) {
  
  strip_legend <- function(p) p + theme(legend.position="none")
  strip_title  <- function(p) p + theme(plot.title=element_blank())
  title_as_annotation <- function(p) {
    ttl <- ggplot_build(p)$plot$labels$title
    p + theme(plot.title=element_blank(),
              plot.tag.position=c(0.6, 0.95),
              plot.tag=element_text(size=12, face='bold', hjust=0.5, vjust=0.5)) +
      labs(tag=ttl)
  }
  row_title <- function(label) {
    ggplot() + 
      labs(title=label) +
      theme_void() +
      theme(plot.title=element_text(size=13, face="bold", hjust=-0.05))
  }
  
  make_exp_fig <- function(con_list, lag_list, main_title) {
    
    # Extract legends from first location
    con_legend <- get_legend(contours$tmp$ALL + 
                               guides(fill = guide_colorbar(
                                 barwidth  = unit(7, "cm"),
                                 barheight = unit(0.5, "cm"),
                                 frame.colour = "black",
                                 ticks.colour = "black"
                               )) +
                               theme(legend.position  = "bottom",
                                     legend.text      = element_text(size = 12),
                                     legend.title     = element_text(size = 14, face='bold')))
    exp_legend <- get_legend(lag_list[[2]][[1]] + theme(legend.position="bottom"))
    legends    <- wrap_plots(list(plot_spacer(), con_legend, exp_legend, plot_spacer()), nrow=1, widths=c(1.5, 1, 1, 1.5))  &
      theme(plot.margin = margin(t=20, b=-10, unit="pt"))
    
    # One row per location
    rows <- lapply(locations, function(l) {
      con  <- strip_title(strip_legend(con_list[[l]]))
      lags <- lapply(lag_list[[l]], function(p) title_as_annotation(strip_legend(p)))
      wrap_plots(c(list(con), lags), nrow=1)
    })
    
    # Add row titles
    rows_titled <- lapply(seq_along(locations), function(i) {
      list(row_title(locations[i]), rows[[i]])
    })
    all_cells <- do.call(c, rows_titled)
    
    # Distinguish between altitudes
    low_t <- get_title('Lower altitude')
    low <- wrap_plots(all_cells[1:8], ncol=2, widths=c(0.2, 1))
    high_t <- get_title('Higher altitude')
    high <- wrap_plots(all_cells[9:16], ncol=2, widths=c(0.2, 1))
    all <- wrap_plots(list(low_t, high_t, low, high), ncol=2, heights=c(0.02,1))
    
    wrap_plots(list(all, legends), ncol=1, heights=c(1, 0.08)) +
      plot_annotation(title=main_title,
                      theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
  }
  
  list(
    tmp = make_exp_fig(contours$tmp, exposure$tmp, "Temperature"),
    tp  = make_exp_fig(contours$tp,  exposure$tp,  "Precipitation")
  )
}

fig <- make_figure_sites(cont, exp)
ggsave(paste0(path_report, 'II. dlnm/MF4_site_exp-lag_tmp.png'), fig$tmp, width=15.5, height=11, dpi=300, bg='white')
ggsave(paste0(path_report, 'II. dlnm/MF4_site_exp-lag_tp.png'),  fig$tp,  width=15.5, height=11, dpi=300, bg='white')


# Model 2 - Subgroup results ####

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
cont_outcomes <- list(tmp=list(), tp=list())
exp_outcomes  <- list(tmp=list(), tp=list())
pooled_meta <- list(tmp = list(), tp = list())

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


# Check significant lags before setting relevant exp-resp
wrap_plots(cont_sig$tmp, ncol=3) +
  plot_layout(guides='collect') +
  plot_annotation(title='Significant exposure-lag-response',
                  theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))) &
  theme(legend.position='bottom')




# MF5: Subgroup exposure-lag-response ONE PLOT ####

get_pooled_df <- function(meta, reference='50%') {
  
  list2env(meta, envir=environment())
  
  if (!is.null(meta$mv)) {
    coef <- coef(mv)
    vcov <- vcov(mv)
    # Remove duplicates row wise
    rows_with_dup <- apply(cen, 1, function(row) any(duplicated(row)))
    cen <- as.data.frame(cen[!rows_with_dup, ])
    perc <- as.numeric(sub("%", "", names(rowMeans(cen))))
    reference <- if(reference=='0%') paste0(min(perc), '%') else reference
  } else {
    city <- names(Slist)
    coef <- as.numeric(ymat[city, ])
    vcov <- Slist[[city]]
    m <- 1
    # Handle duplicates col wise
    cols_with_dup <- duplicated(cen) | duplicated(cen, fromLast=TRUE)
    cen <- cen[!cols_with_dup, , drop=FALSE]
    perc <- as.numeric(sub("%", "", names(rowMeans(cen))))
    reference <- if(reference=='0%') paste0(min(perc), '%') else reference
  }
  
  # Get overall estimates (for specific lag, if applicable)
  cp <- crosspred(b_mean, coef=coef, vcov=vcov, model.link='log', at=rowMeans(cen), cen=rowMeans(cen)[reference])
  
  data.frame(
    pred = as.numeric(sub("%","", names(rowMeans(cen)))),
    RR   = cp$allRRfit,
    low  = cp$allRRlow,
    high = cp$allRRhigh
  )
}

all_df <- map_dfr(outcomes, function(i) {
  
  tmp <- get_pooled_df(output[[i]]$meta_tmp, reference='50%') %>%
    mutate(Subgroup = as.factor(replace_exact_matches(i, nm)), exposure = "Temperature")
  
  tp <- get_pooled_df(output[[i]]$meta_tp, reference='20%') %>%
    mutate(Subgroup = as.factor(replace_exact_matches(i, nm)), exposure = "Precipitation")
  
  bind_rows(tmp, tp)
})

all_df <- all_df %>%
  mutate(pred = as.numeric(pred)) %>%
  filter(
    (exposure == "Temperature" & pred <= 20) |
      (exposure == "Precipitation" & pred >= 80)
  )
row.names(all_df) <- NULL

(p_temp <- ggplot(
  subset(all_df, exposure == "Temperature"),
  aes(x = pred, y = RR, color = Subgroup, group = Subgroup)
) +
    geom_hline(yintercept = 1, color = "gray40") +
    
    geom_ribbon(
      aes(ymin = low, ymax = high, fill = Subgroup),
      alpha = 0.1,
      colour = NA
    ) +
    
    geom_line(linewidth = 1.1) +
    
    scale_color_manual(
      name = NULL,
      values = c('black', scales::hue_pal()(length(unique(all_df$Subgroup))))
    ) +
    
    scale_fill_manual(
      name = NULL,
      values = c('black', scales::hue_pal()(length(unique(all_df$Subgroup))))
    ) +
    
    scale_x_continuous(limits = c(0, 20), breaks = pretty(c(0, 20), 5)) +
    
    labs(
      x = "Temperature (%)",
      y = "RR",
    ) +
    
    theme_classic() +
    theme(
      legend.position = "bottom"
    ))


(p_precip <- ggplot(
  subset(all_df, exposure == "Precipitation"),
  aes(x = pred, y = RR, color = Subgroup, group = Subgroup)
) +
    
    geom_hline(yintercept = 1, color = "gray20") +
    
    geom_ribbon(
      aes(ymin = low, ymax = high, fill = Subgroup),
      alpha = 0.1,
      colour = NA
    ) +
    
    geom_line(linewidth = 1.1) +
    
    scale_color_manual(
      name = NULL,
      values = c('black', scales::hue_pal()(length(unique(all_df$Subgroup))))
    ) +
    
    scale_fill_manual(
      name = NULL,
      values = c('black', scales::hue_pal()(length(unique(all_df$Subgroup))))
    ) +
    
    scale_x_continuous(
      breaks = pretty(subset(all_df, exposure == "Precipitation")$pred, n = 5)
    ) +
    
    scale_y_continuous(
      breaks = pretty(
        c(
          min(subset(all_df, exposure == "Precipitation")$low, na.rm = TRUE),
          max(subset(all_df, exposure == "Precipitation")$high, na.rm = TRUE)
        ),
        n = 5
      )
    ) +
    
    coord_cartesian(
      ylim = c(
        min(subset(all_df, exposure == "Precipitation")$low, na.rm = TRUE),
        max(subset(all_df, exposure == "Precipitation")$high, na.rm = TRUE)
      ),
      xlim = range(subset(all_df, exposure == "Precipitation")$pred, na.rm = TRUE)
    ) +
    
    labs(
      x = "Precipitation (%)",
      y = "RR"
    ) +
    
    theme_classic() +
    
    theme(
      axis.line = element_line(color = "gray20"),
      axis.text = element_text(color = "gray20"),
      axis.ticks = element_line(color = "gray20"),
      legend.position = "bottom"
    ))

legend <- get_legend(p_temp)
p_temp <- p_temp + theme(legend.position = "none")
p_precip <- p_precip + theme(legend.position = "none")

top <- wrap_plots(list(p_temp, p_precip), nrow = 1) 
final_plot <- wrap_plots(
  top,
  legend,
  ncol = 1,
  heights = c(1, 0.1)
)

ggsave(paste0(path_report, 'II. dlnm/MF5_pop_exp-lag.png'), final_plot, width=8, height=5, dpi=300, bg='white')


# MF5: Subgroup exposure-lag-response ####

# Figure function
make_figure_pop <- function(cont_outcomes, exp_outcomes, outcomes) {
  
  strip_legend <- function(p) p + theme(legend.position="none")
  strip_title  <- function(p) p + theme(plot.title=element_blank())
  title_as_annotation <- function(p) {
    ttl <- ggplot_build(p)$plot$labels$title
    p + theme(plot.title=element_blank(),
              plot.tag.position=c(0.6, 0.95),
              plot.tag=element_text(size=12, face='bold', hjust=0.5, vjust=0.5)) +
      labs(tag=ttl)
  }
  row_title <- function(label) {
    ggplot() +
      labs(title=label) +
      theme_void() +
      theme(plot.title=element_text(size=13, face="bold", hjust=-0.05))
  }
  
  make_exp_fig <- function(con_list, lag_list, main_title) {
    
    # Extract legends from first outcome
    # Extract legends from first location
    con_legend <- get_legend(con_list[[1]] + 
                               guides(fill = guide_colorbar(
                                 barwidth  = unit(7, "cm"),
                                 barheight = unit(0.5, "cm"),
                                 frame.colour = "black",
                                 ticks.colour = "black"
                               )) +
                               theme(legend.position  = "bottom",
                                     legend.text      = element_text(size = 12),
                                     legend.title     = element_text(size = 14, face='bold')))
    exp_legend <- get_legend(lag_list[[1]][[1]] + theme(legend.position="bottom"))
    legends    <- wrap_plots(list(plot_spacer(), con_legend, exp_legend, plot_spacer()), nrow=1, widths=c(1.5, 1, 1, 1.5))  &
      theme(plot.margin = margin(t=20, b=-10, unit="pt"))
    
    # One row per outcome
    rows <- lapply(outcomes, function(i) {
      con  <- strip_title(strip_legend(con_list[[i]]))
      lags <- lapply(lag_list[[i]], function(p) title_as_annotation(strip_legend(p)))
      wrap_plots(c(list(con), lags), nrow=1)
    })
    
    # Add row titles
    rows_titled <- lapply(seq_along(outcomes), function(i) {
      list(row_title(replace_exact_matches(outcomes[i], nm)), rows[[i]])
    })
    all_cells <- do.call(c, rows_titled)
    
    all <- wrap_plots(all_cells, ncol=4, widths=c(0.2, 1, 0.2, 1))
    
    wrap_plots(list(all, legends), ncol=1, heights=c(1, 0.05)) +
      plot_annotation(title=main_title,
                      theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
  }
  
  list(
    tmp = make_exp_fig(cont_outcomes$tmp, exp_outcomes$tmp, "Temperature"),
    tp  = make_exp_fig(cont_outcomes$tp,  exp_outcomes$tp,  "Precipitation")
  )
}

# Build and save
fig <- make_figure_pop(cont_outcomes, exp_outcomes, outcomes)

ggsave(paste0(path_report, 'II. dlnm/MF5_pop_exp-lag_tmp.png'), fig$tmp, width=15.5, height=11, dpi=300, bg='white')
ggsave(paste0(path_report, 'II. dlnm/MF5_pop_exp-lag_tp.png'), fig$tp,  width=15.5, height=11, dpi=300, bg='white')


# SFT1: RR with CI for all exposure-lags, across locations ####

make_stacked_excel <- function(tab_list_tmp, tab_list_tp, filepath) {
  
  wb <- openxlsx::createWorkbook()
  
  write_sheet <- function(tab_list, sheet_name) {
    openxlsx::addWorksheet(wb, sheet_name)
    current_row <- 1
    
    for (nm in names(tab_list)) {
      df      <- as.data.frame(tab_list[[nm]]$tab)
      sig_mat <- as.data.frame(tab_list[[nm]]$sig)
      
      # Write title row
      openxlsx::writeData(wb, sheet=sheet_name, x=nm, startRow=current_row, startCol=1)
      openxlsx::addStyle(wb, sheet=sheet_name,
                         style=openxlsx::createStyle(textDecoration="bold", fgFill="grey30", fontColour="white"),
                         rows=current_row, cols=1:ncol(df))
      current_row <- current_row + 1
      
      # Write table with header
      openxlsx::writeData(wb, sheet=sheet_name, x=df, startRow=current_row, startCol=1,
                          headerStyle=openxlsx::createStyle(textDecoration="bold",
                                                            border="TopBottom",
                                                            borderStyle="medium"))
      
      # Style data rows - center align
      openxlsx::addStyle(wb, sheet=sheet_name,
                         style=openxlsx::createStyle(halign="center"),
                         rows=(current_row + 1):(current_row + nrow(df)),
                         cols=2:ncol(df), gridExpand=TRUE)
      
      # Bold significant cells (skip Percentile col, hence col offset + 1)
      for (r in seq_len(nrow(sig_mat))) {
        for (c in 2:ncol(sig_mat)) {
          if (isTRUE(sig_mat[r, c])) {
            openxlsx::addStyle(wb, sheet=sheet_name,
                               style=openxlsx::createStyle(halign="center", textDecoration="bold"),
                               rows=current_row + r, cols=c, stack=TRUE)
          }
        }
      }
      
      current_row <- current_row + nrow(df) + 1
    }
    
    openxlsx::setColWidths(wb, sheet=sheet_name, cols=1, widths=10)
    openxlsx::setColWidths(wb, sheet=sheet_name, cols=2:13, widths=18)
  }
  
  write_sheet(tab_list_tmp, "Temperature")
  write_sheet(tab_list_tp,  "Precipitation")
  
  openxlsx::saveWorkbook(wb, filepath, overwrite=TRUE)
}

get_rr_table <- function(meta, city="ALL", reference="50%") {
  
  res <- city_con_res(meta=meta, city=city, reference=reference,
                      xlab="", title="")$df
  
  # --- Build overall (average across lags) using same logic as city_exp_res ---
  list2env(meta, envir=environment())
  
  if (!is.null(meta$mv)) {
    coef_use <- coef(mv)
    vcov_use <- vcov(mv)
    rows_with_dup <- apply(cen, 1, function(row) any(duplicated(row)))
    cen_clean <- as.data.frame(cen[!rows_with_dup, ])
    perc <- as.numeric(sub("%", "", names(rowMeans(cen_clean))))
    ref_idx <- if (reference == '0%') paste0(min(perc), '%') else reference
    cp_at  <- rowMeans(cen_clean)
    cp_cen <- rowMeans(cen_clean)[ref_idx]
  } else {
    city_name <- if (city == "ALL") rownames(ymat)[1] else city
    coef_use <- as.numeric(ymat[city_name, ])
    vcov_use <- Slist[[city_name]]
    cols_with_dup <- duplicated(cen) | duplicated(cen, fromLast=TRUE)
    cen_clean <- cen[!cols_with_dup, , drop=FALSE]
    perc <- as.numeric(sub("%", "", names(rowMeans(cen_clean))))
    ref_idx <- if (reference == '0%') paste0(min(perc), '%') else reference
    cp_at  <- rowMeans(cen_clean)
    cp_cen <- rowMeans(cen_clean)[ref_idx]
  }
  
  cp <- crosspred(b_mean, coef=coef_use, vcov=vcov_use,
                  model.link='log', at=cp_at, cen=cp_cen)
  
  overall_df <- data.frame(
    x      = perc,
    z      = cp$allRRfit,
    z_low  = cp$allRRlow,
    z_high = cp$allRRhigh,
    row.names = NULL
  )
  
  # --- Filter to same percentile rows used in the lag table ---
  filtered <- res %>%
    dplyr::mutate(
      rr_fmt = sprintf("%.2f (%.2f, %.2f)", z, z_low, z_high),
      sig = !(z_low <= 1 & z_high >= 1) & round(z, 2) != 1.00
    ) %>%
    dplyr::filter(
      x <= 5  |
        x >= 95 |
        (x >= 5  & x <= 20 & x %% 5  == 0) |
        (x >= 80 & x <= 95 & x %% 5  == 0) |
        (x >= 20 & x <= 80 & x %% 10 == 0)
    )
  
  # Filter overall to same percentile rows
  overall_filtered <- overall_df %>%
    dplyr::filter(x %in% filtered$x) %>%
    dplyr::mutate(
      rr_fmt = sprintf("%.2f (%.2f, %.2f)", z, z_low, z_high),
      sig = !(z_low <= 1 & z_high >= 1) & round(z, 2) != 1.00
    )
  
  # --- Formatted lag table ---
  tab <- filtered %>%
    dplyr::select(x, y, rr_fmt) %>%
    tidyr::pivot_wider(names_from=y, values_from=rr_fmt, names_prefix="Lag=") %>%
    dplyr::rename(Percentile=x) %>%
    dplyr::left_join(
      overall_filtered %>% dplyr::select(x, rr_fmt) %>% dplyr::rename(Percentile=x, Cumulative=rr_fmt),
      by="Percentile"
    )
  
  # --- Significance matrix ---
  sig_mat <- filtered %>%
    dplyr::select(x, y, sig) %>%
    tidyr::pivot_wider(names_from=y, values_from=sig, names_prefix="Lag=") %>%
    dplyr::rename(Percentile=x) %>%
    dplyr::left_join(
      overall_filtered %>% dplyr::select(x, sig) %>% dplyr::rename(Percentile=x, Cumulative=sig),
      by="Percentile"
    )
  
  list(tab=tab, sig=sig_mat)
}

# Get site specific temperature RR tables
tab_rr_site_tmp <- lapply(c("ALL", locations), function(l) {
  get_rr_table(output$total_patients$meta_tmp, city=l)
})
names(tab_rr_site_tmp) <- c("Pooled", locations)

# Get site specific precipitation RR tables
tab_rr_site_tp <- lapply(c("ALL", locations), function(l) {
  get_rr_table(output$total_patients$meta_tp, city=l)
})
names(tab_rr_site_tp) <- c("Pooled", locations)

# Export
make_stacked_excel(tab_rr_site_tmp, tab_rr_site_tp, filepath=paste0(path_report, "II. dlnm/ST1_rr_site.xlsx"))

# SFT2: RR with CI for all exposure-lags, across population groups ####

# Get outcome specific temperature RR tables
tab_rr_pop_tmp <- lapply(outcomes[-1], function(o) {
  get_rr_table(output[[o]]$meta_tmp, city='ALL')
}) 
names(tab_rr_pop_tmp) <- replace_exact_matches(outcomes[-1], nm)

# Get outcome specific precipitation RR tables
tab_rr_pop_tp <- lapply(outcomes[-1], function(o) {
  get_rr_table(output[[o]]$meta_tp, city='ALL')
}) 
names(tab_rr_pop_tp) <- replace_exact_matches(outcomes[-1], nm)

# Export
make_stacked_excel(tab_rr_pop_tmp, tab_rr_pop_tp, filepath=paste0(path_report, "II. dlnm/ST2_rr_pop.xlsx"))


