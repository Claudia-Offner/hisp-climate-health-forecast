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

# Load model 1 results ####
results <- readRDS(paste0('2_DLNM/pooled_model.RData'))
output <- results$output
cont <- results$cont
exp <- results$exp
rm(results)

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
  popu = "Population density",
  rh   = "Relative humidity"
)
(tbl_tmp <- make_meta_reg_table(meta=output$meta_tmp, predictors_list))
(tbl_tp  <- make_meta_reg_table(output$meta_tp,  predictors_list))

# Format flextable
d_tmp <- prep_display(tbl_tmp)
d_tp  <- prep_display(tbl_tp)

combined <- rbind(section_row(), d_tmp, section_row(), d_tp)
tmp_rows   <- 2:5
tp_rows    <- 7:10
label_rows <- c(1, 6)
intercept_rows <- c(2, 7)   # italic "Intercept only" rows

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
          value = as_paragraph(as_chunk("Temperature", props = fp_text(bold = TRUE, font.family = "Times New Roman", font.size = 9)))) |>
  compose(i = 6, j = 1, part = "body",
          value = as_paragraph(as_chunk("Precipitation", props = fp_text(bold = TRUE, font.family = "Times New Roman", font.size = 9)))) |>
  
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
  fontsize(size = 9, part = "all") |>
  
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
  hline(i = 5, part = "body",   border = fp_border(width = 0.75)) |>
  
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
  "LR = likelihood ratio test;",
  "AIC = Akaike information criterion;",
  "Q = Cochran's Q statistic;",
  "I$^2$ = percentage of residual heterogeneity."
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
  paste0("\\textit{Note.} ", footnote_text),
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
    
    all <- wrap_plots(all_cells, ncol=4, widths=c(0.2, 1, 0.2, 1))
    
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



################

# 
# tot_contours <- function(sig=FALSE, max_dev_tmp=NULL, max_dev_tp=NULL){
#   
#   # CONTOURS
#   plot_contours <- function(meta, locations, xlab, title, reference='50%', sig=sig, 
#                             n_col=2, max_dev=NULL) {
#     
#     args <- list(meta=meta, xlab=xlab, reference=reference, sig=sig)
#     
#     # Shared limits
#     dfs <- c(
#       list(ALL=do.call(city_con_res, c(args, list(city='ALL', title='Overall')))$df),
#       setNames(lapply(locations, function(l) do.call(city_con_res, c(args, list(city=l, title=l)))$df), locations)
#     )
#     
#     if (is.null(max_dev)) max_dev <- max(abs(pretty(range(do.call(rbind, dfs)$z, na.rm=TRUE), n=15) - 1))
#     shared_lims <- c(1 - max_dev, 1 + max_dev)
#     shared_brks <- pretty(shared_lims, n=12)
#     args <- c(args, list(shared_lims=shared_lims, shared_brks=shared_brks))
#     
#     # Plots
#     contours <- list()
#     for (l in locations) {
#       contours[[l]] <- do.call(city_con_res, c(args, list(city=l, title=l)))$p
#     }
#     
#     # Layout
#     n_rows <- ceiling(length(contours) / n_col)
#     contours_mod <- lapply(seq_along(contours), function(i) {
#       row <- ceiling(i / n_col); col <- ((i - 1) %% n_col) + 1
#       contours[[i]] +
#         theme(
#           axis.title.x=if(row == n_rows) element_text() else element_blank(),
#           axis.title.y=if(col == 1)      element_text() else element_blank()
#         )
#     })
#     
#     wrap_plots(contours_mod, ncol=n_col) +
#       plot_layout(guides='collect') +
#       plot_annotation(title=title,
#                       theme=theme(plot.title=element_text(hjust=0.5, face='bold', size=14))) &
#       theme(legend.position='none')  
#   }
#   
#   # --- Compute shared max_dev across pooled AND city-specific ---
#   compute_max_dev <- function(meta, locations, reference, args_base) {
#     all_dfs <- c(
#       list(do.call(city_con_res, c(args_base, list(city='ALL', title='Overall')))$df),
#       lapply(locations, function(l) do.call(city_con_res, c(args_base, list(city=l, title=l)))$df)
#     )
#     max(abs(pretty(range(do.call(rbind, all_dfs)$z, na.rm=TRUE), n=15) - 1))
#   }
#   
#   args_tmp <- list(meta=output$meta_tmp, xlab='Temperature (%)',   reference='50%', sig=sig)
#   args_tp  <- list(meta=output$meta_tp,  xlab='Precipitation (%)', reference='20%', sig=sig)
#   
#   if (is.null(max_dev_tmp)) max_dev_tmp <- compute_max_dev(output$meta_tmp, locations, '50%', args_tmp)
#   if (is.null(max_dev_tp)) max_dev_tp <- compute_max_dev(output$meta_tp, locations, '50%', args_tmp)
# 
#   # Shared limits for pooled plots
#   make_shared <- function(max_dev) {
#     lims <- c(1 - max_dev, 1 + max_dev)
#     list(shared_lims=lims, shared_brks=pretty(lims, n=10))
#   }
# 
#   # Overall pooled plots
#   tmp <- do.call(city_con_res, c(args_tmp, make_shared(max_dev_tmp), list(city='ALL', title=NULL)))$p 
#   tp  <- do.call(city_con_res, c(args_tp,  make_shared(max_dev_tp),  list(city='ALL', title=NULL)))$p
#   
#   # Extract legend from one plot using cowplot
#   # shared_legend <- get_legend(tmp + theme(legend.position='bottom', legend.text=element_text(angle=45, hjust=1)))
#   # 
#   # # Get title
#   # tmp <- wrap_plots(get_title('Overall temperature'), tmp, nrow=2, heights=c(1, 40)) + 
#   #   theme(legend.position='none')
#   # tp  <- wrap_plots(get_title('Overall precipitation'), tp, nrow=2, heights=c(1, 40)) + 
#   #   theme(legend.position='none')
#   # 
#   # # City-specific plots
#   # tmp_l <- plot_contours(meta=output$meta_tmp, locations=locations,
#   #                        xlab='Temperature (%)', title='Maximum Temperature',
#   #                        reference='50%', sig=sig, max_dev=max_dev_tmp) 
#   # tp_l  <- plot_contours(meta=output$meta_tp, locations=locations,
#   #                        xlab='Precipitation (%)', title='Total Precipitation',
#   #                        reference='20%', sig=sig, max_dev=max_dev_tp)
#   # tmp_l <- wrap_plots(get_title('Site-specific temperature'), tmp_l, nrow=2, heights=c(1, 40))
#   # tp_l  <- wrap_plots(get_title('Site-specific precipitation'), tp_l, nrow=2, heights=c(1, 40))
#   # 
#   # con_tmp <- wrap_plots(list(tmp, tmp_l), ncol=1) + plot_layout(heights=c(1, 2))
#   # con_tp  <- wrap_plots(list(tp,  tp_l),  ncol=1) + plot_layout(heights=c(1, 2))
#   # combined <- wrap_plots(list(con_tmp, con_tp), nrow=1)
#   # 
#   # plot_grid(combined, shared_legend, ncol=1, rel_heights=c(1, 0.08))
#   
# }
# 
# dev <- 0.16
# contour_all <- tot_contours(sig=FALSE, max_dev_tmp=dev, max_dev_tp=dev)
# contour_all_sig <- tot_contours(sig=TRUE, max_dev_tmp=dev, max_dev_tp=dev)
# 
# ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_ALL.png'), plot=contour_all, width=7, height=11, dpi=300, bg='white')
# ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_ALL_sig.png'), plot=contour_all_sig, width=7, height=11, dpi=300, bg='white')


# 
# exp <- function(output, city, plot_key){
#   
#   city_tmp_lags <- tmp_lags[[city]]
#   city_tp_lags  <- tp_lags[[city]]
#   
#   # Build tmp plots
#   tmp_plots <- lapply(city_tmp_lags, function(lag) {
#     city_exp_res(city=city, lag=lag, meta=output$meta_tmp,
#                  col='#D74B3C', xlab='Temperature (%)', 
#                  title=paste0('Lag=', lag), reference='50%')[[plot_key]]
#   })
#   # Build tp plots
#   tp_plots <- lapply(city_tp_lags, function(lag) {
#     city_exp_res(city=city, lag=lag, meta=output$meta_tp,
#                  col='#1D98EB', xlab='Precipitation (%)', 
#                  title=paste0('Lag=', lag), reference='0%')[[plot_key]]
#   })
#   
#   all_plots <- c(tmp_plots, tp_plots)
#   
#   # Pad to exactly 5 with void placeholders
#   while (length(all_plots) < 5) all_plots <- c(all_plots, list(ggplot() + theme_void()))
#   
#   # Strip legends
#   all_plots <- lapply(all_plots, function(p) p + theme(legend.position='none'))
#   
#   return(all_plots)  # return list of exactly 5 plots
# }
# 
# # --- Build label + 5 plots per row, assembled into one aligned grid ---
# locations_order <- c('Chamba', 'Dharamshala', 'Solan', 'Hamirpur', 'Kullu', 
#                      'Reckong Peo', 'Nahan', 'Una')
# 
# # Collect all rows as flat list of (label + 5 plots)
# all_cells <- list()
# for (l in locations_order) {
#   
#   # Centered city label as a single cell
#   label_cell <- ggplot() +
#     annotate('text', x=0.5, y=0.5, label=l, 
#              fontface='bold', size=4.5, hjust=0.5, vjust=0.5) +
#     theme_void()
#   
#   plot_cells <- exp(output, l, plot_key='p_city')  # list of exactly 5
#   all_cells  <- c(all_cells, list(label_cell), plot_cells)  # 6 cells per row
# }
# 
# # Assemble into grid: 6 cols (1 label + 5 plots), n rows
# combined <- wrap_plots(all_cells, ncol=6, byrow=TRUE) +
#   plot_layout(widths=c(10, rep(10, 5)))  # label col narrow, plot cols equal
# 
# # --- Pooled ---
# exp_pool_plots <- exp(output, 'ALL', plot_key='p_blup')
# exp_pool <- wrap_plots(exp_pool_plots, nrow=1)
# 
# # --- Unified legend ---
# unified_legend <- ggplot(data.frame(x=1, y=1), aes(x, y)) +
#   geom_line(aes(linetype='BLUP'),             colour='black') +
#   geom_line(aes(linetype='Pooled'),           colour='black') +
#   geom_line(aes(colour='Raw temperature'),    linewidth=1) +
#   geom_line(aes(colour='Raw precipitation'),  linewidth=1) +
#   scale_linetype_manual(name=NULL, values=c('BLUP'='dashed', 'Pooled'='solid')) +
#   scale_colour_manual(name=NULL, values=c('Raw temperature'='#D74B3C', 'Raw precipitation'='#1D98EB')) +
#   theme_void() +
#   theme(legend.position='bottom',
#         legend.direction='horizontal',
#         legend.text=element_text(size=10))
# shared_legend <- get_legend(unified_legend)
# 
# 
# exp_pool <- plot_grid(exp_pool, shared_legend, ncol=1, rel_heights=c(40, 1))
# exp_city <- plot_grid(combined, shared_legend, ncol=1, rel_heights=c(40, 1))

# Population Contours & Exposure responses

# pop_contours <- function(outcomes, output, tmp, tp, all_z_tmp, all_z_tp,max_dev=NULL, sig=FALSE){
#   
#   # --- Step 2: Compute shared limits across ALL outcomes ---
#   make_shared <- function(all_z, max_dev=0.18) {
#     if (is.null(max_dev)) max_dev <- max(abs(pretty(range(all_z, na.rm=TRUE), n=15) - 1))
#     shared_lims <- c(1 - max_dev, 1 + max_dev)
#     shared_brks <- pretty(shared_lims, n=12)
#     list(shared_lims=shared_lims, shared_brks=shared_brks)
#   }
#   
#   scale_tmp <- make_shared(all_z_tmp, max_dev=max_dev)
#   scale_tp  <- make_shared(all_z_tp, max_dev=max_dev)
#   
#   # --- Step 3: Build plots using shared limits ---
#   for (i in outcomes){
#     args_tmp <- c(list(meta=output[[i]]$meta_tmp, xlab='Temperature (%)',   reference='50%', sig=sig), scale_tmp)
#     args_tp  <- c(list(meta=output[[i]]$meta_tp,  xlab='Precipitation (%)', reference='20%', sig=sig), scale_tp)
#     
#     tmp[[i]] <- do.call(city_con_res, c(args_tmp, list(city='ALL', title=replace_exact_matches(i, nm))))$p
#     tp[[i]]  <- do.call(city_con_res, c(args_tp,  list(city='ALL', title=replace_exact_matches(i, nm))))$p
#   }
#   
#   # --- Step 4: Combine ---
#   con_tmp <- wrap_plots(tmp, ncol=4) + 
#     plot_layout(guides='collect') + 
#     plot_annotation(title='Overall maximum temperature', 
#                     theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))) &
#     theme(legend.position='none',
#           legend.text=element_text(angle=45, hjust=1))  
#   con_tp <- wrap_plots(tp, ncol=4) + 
#     plot_layout(guides='collect') + 
#     plot_annotation(title='Overall total precipitation', 
#                     theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))) &
#     theme(legend.position='none',
#           legend.text=element_text(angle=45, hjust=1))
#   
#   con_tmp <- wrap_plots(get_title('Pooled maximum temperature'), con_tmp, nrow=2, heights=c(1, 40))
#   con_tp <- wrap_plots(get_title('Pooled 7-day cummulative precipitation'), con_tp, nrow=2, heights=c(1, 40))
#   
#   # Extract legend from one plot using cowplot
#   shared_legend <- get_legend(
#     tmp[[1]] + theme(legend.position='bottom',
#                      legend.text=element_text(angle=45, hjust=1))
#   )
#   
#   combined <- wrap_plots(list(con_tmp, con_tp), nrow=2) + plot_annotation(theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold")))
#   plot_grid(combined, shared_legend, ncol=1, rel_heights=c(1, 0.08))
#   
# }
# 
# outcomes <- h_cols[1:8]
# output <- list()
# tmp    <- list()
# tp     <- list()
# pop_exp <- list()
# 
# # --- Step 1: Run models and collect all z-values first ---
# all_z_tmp <- c()
# all_z_tp  <- c()
# 
# for (i in outcomes){
#   print(i)
#   res <- run_bestfit(df, locations, i, params=final_params$pooled, pooled=TRUE, city=NULL)
#   output[[i]] <- res$output
#   
#   args_tmp <- list(meta=output[[i]]$meta_tmp, xlab='Temperature (%)',   reference='50%', sig=F)
#   args_tp  <- list(meta=output[[i]]$meta_tp,  xlab='Precipitation (%)', reference='20%', sig=F)
#   
#   tp_lags$ALL <- c(3,10)
#   if(i %in% c('M', 'age25_64', 'over_65')){
#     tp_lags$ALL <- c(10)
#   } else if (i == 'age5_14'){
#     tp_lags$ALL <- c()
#   }
#   
#   # exp_pool_plots <- exp(output[[i]] , 'ALL', plot_key='p_blup')
#   # exp_pool_plots <- wrap_plots(exp_pool_plots, nrow=1)
#   # pop_exp[[i]] <- wrap_plots(get_title(replace_exact_matches(i, nm), size=14), exp_pool_plots, nrow=2, heights=c(1, 40))
#   # With:
#   exp_pool_plots <- exp(output[[i]], 'ALL', plot_key='p_blup')
#   
#   # Pad to exactly 5 slots
#   while (length(exp_pool_plots) < 5) exp_pool_plots <- c(exp_pool_plots, list(ggplot() + theme_void()))
#   exp_pool_plots <- exp_pool_plots[1:5]  # ensure never more than 5
#   
#   exp_pool_plots <- wrap_plots(exp_pool_plots, nrow=1, ncol=5)
#   pop_exp[[i]] <- wrap_plots(get_title(replace_exact_matches(i, nm), size=14), exp_pool_plots, nrow=2, heights=c(1, 40))
#   
#   # Collect z ranges across all outcomes
#   all_z_tmp <- c(all_z_tmp, do.call(city_con_res, c(args_tmp, list(city='ALL', title=i)))$df$z)
#   all_z_tp  <- c(all_z_tp,  do.call(city_con_res, c(args_tp,  list(city='ALL', title=i)))$df$z)
# }
# 
# # Get contours
# contour_pop <- pop_contours(outcomes, output, tmp, tp, all_z_tmp, all_z_tp, max_dev=0.18, sig=FALSE)
# contour_pop_sig <- pop_contours(outcomes, output, tmp, tp, all_z_tmp, all_z_tp, max_dev=0.18, sig=TRUE)
# 
# # Get exp-response
# pop_exp_all <- wrap_plots(pop_exp, ncol=2)
# pop_exp_all <- plot_grid(pop_exp_all, shared_legend, ncol=1, rel_heights=c(40, 1))
# 
# ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_pop.png'), plot=contour_pop, width=10, height=12, dpi=300, bg='white')
# ggsave(filename=paste0(path_report, 'II. dlnm/fig_contour_pop_sig.png'), plot=contour_pop_sig, width=10, height=12, dpi=300, bg='white')
# ggsave(filename=paste0(path_report, 'II. dlnm/fig_pop_exp.png'), plot=pop_exp_all, width=15, height=10, dpi=300, bg='white')
# 

