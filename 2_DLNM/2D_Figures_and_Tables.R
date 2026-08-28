### ------------------------------------------------------------------------ ### 
###  DNLM – Tables & Figures (HISP)
### ------------------------------------------------------------------------ ### 


#### Dependencies ####

# LaTeX table builder
get_kbl_tab <- function(dat, title, note, colnames, label,
                        spanning_header = NULL, bold_rows = NULL,
                        italic_rows = NULL, hline_after_rows = NULL) {
  k <- kbl(dat, format = "latex", booktabs = TRUE, longtable = TRUE,
           linesep = "", col.names = colnames, escape = FALSE,
           caption = title, label = label) |>
    kable_styling(latex_options = c("repeat_header", "hold_position"),
                  font_size = 9, full_width = FALSE, position = "center")
  
  if (!is.null(spanning_header))  k <- k |> add_header_above(spanning_header, escape = FALSE)
  if (!is.null(bold_rows))        k <- k |> row_spec(bold_rows,   bold   = TRUE)
  if (!is.null(italic_rows))      k <- k |> row_spec(italic_rows, italic = TRUE)
  if (!is.null(hline_after_rows))
    for (r in hline_after_rows)   k <- k |> row_spec(r, extra_latex_after = "\\midrule")
  
  k |> footnote(general = note, general_title = "",
                footnote_as_chunk = TRUE, escape = FALSE)
}

# Escape Unicode characters to LaTeX math
escape_latex_unicode <- function(x) {
  x <- gsub("±", "$\\pm$",    x, fixed = TRUE)
  x <- gsub("°", "$^\\circ$", x, fixed = TRUE)
  x <- gsub("µ", "$\\mu$",    x, fixed = TRUE)
  x <- gsub("³", "$^{3}$",    x, fixed = TRUE)
  x
}

# Empty ggplot carrying only a bold title (patchwork spacer workaround)
get_title <- function(title, size = 16) {
  ggplot() + labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = size, face = "bold"))
}

# Formatting helpers for numeric table cells
fmt_p   <- function(x) ifelse(is.na(x), "\u2013", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_num <- function(x, d = 1) ifelse(is.na(x), "\u2013", formatC(x, format = "f", digits = d))
fmt_int <- function(x) ifelse(is.na(x), "\u2013", as.character(as.integer(x)))

#### MAIN CODE ####

export <-a

# Load data
load(paste0(path_data, '2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata'))

df_wide <- df_wide |>
  filter(!location %in% c('Bilaspur'),
         date >= as.Date('2015-02-14'),
         date <= as.Date('2020-02-14')) |>
  rename(rh = rh_mean)

h_cols <- c('closed', 'total_patients', 'location', 'F', 'M', 'under_5',
            'age5_14', 'age15_24', 'age25_64', 'over_65', 'general',
            'ophthalmology') #'ent', 'gynaecology')

#### FIGURES ####

## Shared plot-manipulation helpers
strip_legend    <- function(p) p + theme(legend.position = "none")
strip_title     <- function(p) p + theme(plot.title = element_blank())
tag_as_title    <- function(p) {
  ttl <- ggplot_build(p)$plot$labels$title
  strip_title(p) +
    theme(plot.tag.position = c(0.6, 0.95),
          plot.tag = element_text(size = 12, face = 'bold', hjust = 0.5, vjust = 0.5)) +
    labs(tag = ttl)
}
row_label       <- function(label) {
  ggplot() + labs(title = label) + theme_void() +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = -0.05))
}
contour_legend  <- function(ref_plot, bar_w = 3.5, bar_h = 0.3, txt = 7, ttl = 8) {
  get_legend(ref_plot +
               guides(fill = guide_colorbar(barwidth  = unit(bar_w, "cm"), barheight = unit(bar_h, "cm"),
                                            frame.colour = "black", ticks.colour = "black")) +
               theme(legend.position = "bottom",
                     legend.text  = element_text(size = txt),
                     legend.title = element_text(size = ttl, face = 'bold')))
}
sig_legend      <- function(){
  
  leg <- get_legend(ggplot(data.frame(x = 1, y = 1, grp = "Significance"), aes(x = x, y = y)) +
    geom_tile_pattern(
      aes(pattern_fill = grp), 
      pattern = "crosshatch",
      fill = "white", pattern_colour = NA, color = NA, 
      pattern_density = 0.1, pattern_spacing = 0.03, 
      pattern_angle = 45, pattern_alpha=NA
    ) +
    scale_pattern_fill_manual(name = NULL, values = c("Significance" = "grey50")) +
    theme_void() +
    theme(legend.position = "bottom"))
  
  return(leg)

  }
title_only      <- function(title){
  ggplot() +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size=16, face='bold'))
}
make_panel <- function(plots, title=''){
  wrap_plots(plots, nrow = 2) +
    plot_layout(guides = 'collect') +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 14))
    ) & theme(legend.position = 'none')
}
prep_display <- function(df) {
  data.frame(
    `Meta-predictors` = df[["Meta-predictors"]],
    `Stat`   = fmt_num(df$LR_stat), `df`  = fmt_int(df$LR_df),  `p`    = fmt_p(df$LR_p),
    `AIC`    = fmt_num(df$AIC),
    `Stat `  = fmt_num(df$Q_stat),  `df ` = fmt_int(df$Q_df),   `p `   = fmt_p(df$Q_p),
    `I2 (%)` = fmt_num(df$I2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
get_meta_plot <- function(table){
  
  booktabs_theme <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0, x = 0.05, fontsize = 10.5),
      bg_params = list(fill = NA, col = NA)          # no shading, no cell borders
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10.5),
      bg_params = list(fill = NA, col = NA)
    ),
    rowhead = list(fg_params = list(fontsize = 10.5))
  )
  
  tab_grob <- tableGrob(table, rows = NULL, theme = booktabs_theme)
  span_labels <- c("", "LR Test", 'Fit', 'Cochran Q test', "I2")
  span_widths <- c(1, 3, 1, 3, 1)
  stopifnot(sum(span_widths) == ncol(table))
  
  # explicit, generous row height instead of borrowing from a data cell
  tab_grob <- gtable_add_rows(tab_grob, heights = grid::unit(1.2, "cm"), pos = 0)
  
  col_ends   <- cumsum(span_widths)
  col_starts <- c(1, head(col_ends, -1) + 1)
  
  for (i in seq_along(span_labels)) {
    if (span_labels[i] == "") next
    tab_grob <- gtable_add_grob(
      tab_grob,
      grobs = textGrob(
        span_labels[i],
        gp = gpar(fontface = "bold", fontsize = 10.5),
        just = "center"
      ),
      t = 1, l = col_starts[i], r = col_ends[i],
      name = paste0("span-", i)
    )
    tab_grob <- gtable_add_grob(
      tab_grob,
      grobs = segmentsGrob(x0 = 0.05, x1 = 0.95, y0 = 0.1, y1 = 0.1, gp = gpar(lwd = 1.5)),
      t = 1, l = col_starts[i], r = col_ends[i],
      name = paste0("span-line-", i)
    )
  }
  return(tab_grob)
}
results1 <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_pooled_plots.RData'))
results2 <- readRDS(paste0('2_DLNM/Results/',  export, '_', model2, '_pooled_plots.RData'))

# Get all exposure-response results (t2m, tp, mn2t, mx2t) ----

make_exp_fig <- function(cont1, cont2, exp1, exp2) {
  
  imm_ov <- strip_legend(exp1$ALL$Overall)
  cum_ov <- strip_legend(exp2$ALL$Overall)
  exp1$ALL$Overall <- NULL
  exp2$ALL$Overall <- NULL
  con_imm <- strip_legend(cont1$ALL)
  con_cum  <- strip_legend(cont2$ALL)
  lag_imm <- lapply(exp1$ALL, function(p) tag_as_title(strip_legend(p)))
  lag_cum  <- lapply(exp2$ALL,  function(p) tag_as_title(strip_legend(p)))
  
  legends <- wrap_plots(
    list(plot_spacer(),
         contour_legend(cont1$ALL),
         # sig_legend(),
         get_legend(exp1$ALL[[1]] + theme(legend.position = "bottom"))),
    nrow = 1, widths = c(0.25, 1, 4)
  )
  
  all_rows <- wrap_plots(
    list(row_label("(A)"), wrap_plots(c(list(con_imm), imm_ov, lag_imm), nrow = 1),
         row_label("(B)"), wrap_plots(c(list(con_cum), cum_ov, lag_cum),  nrow = 1)),
    nrow = 2, widths = c(0.05, 1)
  )
  wrap_plots(list(all_rows, legends), ncol = 1, heights = c(1, 0.1))
}

if(export == 't2m'){
  
  cont1 <- results1$cont$tmp; exp1 <- results1$exp$tmp
  cont2 <- results2$cont$tmp; exp2 <- results2$exp$tmp
  fig <- make_exp_fig(cont1, cont2, exp1, exp2)
  ggsave(paste0(path_report, 'II. dlnm/MF_',  export, '_pool.png'), fig, width = 14, height = 7, dpi = 300, bg = 'white')
  
  cont1 <- results1$cont$tp; exp1 <- results1$exp$tp
  cont2 <- results2$cont$tp; exp2 <- results2$exp$tp
  fig <- make_exp_fig(cont1, cont2, exp1, exp2)
  ggsave(paste0(path_report, 'II. dlnm/MF_tp_pool.png'), fig, width = 14, height = 7, dpi = 300, bg = 'white')
  
  # SITES
  results1$cont$tmp$ALL <- NULL
  results1$cont$tp$ALL <- NULL
  results2$cont$tmp$ALL <- NULL
  results2$cont$tp$ALL <- NULL
  site_imm_tmp <- make_panel(results1$cont$tmp)
  site_cum_tmp <- make_panel(results2$cont$tmp)
  site_imm_tp <- make_panel(results1$cont$tp)
  site_cum_tp <- make_panel(results2$cont$tp)
  
  # META TESTS
  output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_meta_tables.RData'))
  meta_imm_tmp <- get_meta_plot(prep_display(output$tbl_tmp))
  meta_imm_tp <- get_meta_plot(prep_display(output$tbl_tp))
  output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model2, '_meta_tables.RData'))
  meta_cum_tmp <- get_meta_plot(prep_display(output$tbl_tmp))
  meta_cum_tp <- get_meta_plot(prep_display(output$tbl_tp))

  all_rows <- wrap_plots(
    list(row_label("(A)"), wrap_plots(c(site_imm_tmp, wrap_elements(full = meta_imm_tmp)), nrow = 1, widths=c(2, 1)),
         row_label("(B)"), wrap_plots(c(site_cum_tmp, wrap_elements(full = meta_cum_tmp)), nrow = 1, widths=c(2, 1)),
         row_label("(C)"), wrap_plots(c(site_imm_tp, wrap_elements(full = meta_imm_tp)), nrow = 1, widths=c(2, 1)),
         row_label("(D)"), wrap_plots(c(site_cum_tp, wrap_elements(full = meta_cum_tp)), nrow = 1, widths=c(2, 1))),
    nrow = 4, widths = c(0.05, 1)
  )
  
  legends <- wrap_plots(list(plot_spacer(), contour_legend(cont1$ALL), plot_spacer()), nrow = 1, widths = c(0.05, 2, 1))
  fig <- wrap_plots(list(legends, all_rows), ncol = 1, heights = c( 0.025, 1))
  
  ggsave(
    paste0(path_report, 'II. dlnm/SF_', export, '_site.png'),
    fig, width = 20, height = 22, dpi = 100, bg = 'white'
  )
  
} else {
  
  cont1 <- results1$cont$tmp; exp1 <- results1$exp$tmp
  cont2 <- results2$cont$tmp; exp2 <- results2$exp$tmp
  
  # POOLED
  imm_ov <- strip_legend(exp1$ALL$Overall)
  cum_ov <- strip_legend(exp2$ALL$Overall)
  exp1$ALL$Overall <- NULL
  exp2$ALL$Overall <- NULL
  con_imm <- strip_legend(cont1$ALL)
  con_cum  <- strip_legend(cont2$ALL)
  lag_imm <- lapply(exp1$ALL, function(p) tag_as_title(strip_legend(p)))
  lag_cum  <- lapply(exp2$ALL,  function(p) tag_as_title(strip_legend(p)))
  
  # SITES
  results1$cont$tmp$ALL <- NULL
  results1$cont$tp$ALL <- NULL
  results2$cont$tmp$ALL <- NULL
  results2$cont$tp$ALL <- NULL
  site_imm_tmp <- make_panel(results1$cont$tmp)
  site_cum_tmp <- make_panel(results2$cont$tmp)
  
  # META TESTS
  output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_meta_tables.RData'))
  imm_tmp   <- prep_display(output$tbl_tmp)
  meta_imm <- get_meta_plot(imm_tmp)
  output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model2, '_meta_tables.RData'))
  cum_tmp   <- prep_display(output$tbl_tmp)
  meta_cum <- get_meta_plot(cum_tmp)
  
  all_rows <- wrap_plots(
    list(row_label("(A)"), wrap_plots(c(list(con_imm), imm_ov, lag_imm), nrow = 1),
         row_label("(B)"), wrap_plots(c(list(con_cum), cum_ov, lag_cum),  nrow = 1),
         row_label("(C)"), wrap_plots(c(site_imm_tmp, wrap_elements(full = meta_imm)), nrow = 1, widths=c(2, 1)),
         row_label("(D)"), wrap_plots(c(site_cum_tmp, wrap_elements(full = meta_cum)), nrow = 1, widths=c(2, 1))),
    nrow = 4, widths = c(0.05, 1), heights = c(1, 1, 2, 2)
  )
  
  legends <- wrap_plots(
    list(plot_spacer(),
         contour_legend(cont1$ALL),
         # sig_legend(),
         get_legend(exp1$ALL[[1]] + theme(legend.position = "bottom"))),
    nrow = 1, widths = c(0.25, 1, 4)
  )
  fig <- wrap_plots(list(legends, all_rows), ncol = 1, heights = c( 0.025, 1))
  
  ggsave(
    paste0(path_report, 'II. dlnm/SF_', export, '_res.png'),
    fig, width = 20, height = 20, dpi = 100, bg = 'white'
  )
  
}








# # Time series visual ####
# 
# if(export == 't2m'){
#   
#   get_time_series <- function(df_weekly, title=NULL, month_breaks='4 months', hjust.left=0.25, hjust.right=0.2,
#                               hosp_bottom=400, band_bottom=800, weather_min=-10, weather_max=40){
#     
#     # Weather occupies the band from band_bottom (i.e. anchor) up to hosp_top
#     hosp_top <- round(band_bottom + (band_bottom/2.667))
#     
#     band_top <- hosp_top
#     
#     # Set breaks
#     hosp_breaks <- seq(hosp_bottom, hosp_top_full, by = 100)
#     hosp_labels <- ifelse(hosp_breaks > hosp_top_visible, "", hosp_breaks)
#     
#     transform_weather <- function(x) {
#       band_bottom + (x - weather_min) / (weather_max - weather_min) * (band_top - band_bottom)
#     }
#     untransform_weather <- function(y) {
#       (y - band_bottom) / (band_top - band_bottom) * (weather_max - weather_min) + weather_min
#     }
#     
#     # Generate breaks every 4 months, then drop anything past the last real data point
#     max_date <- max(df_weekly$date, na.rm = TRUE)
#     all_breaks <- seq(min(df_weekly$date, na.rm = TRUE), max_date + months(4), by = month_breaks)
#     valid_breaks <- all_breaks[all_breaks <= max_date]
#     
#     weat_breaks <- seq(weather_min, weather_max, by = 10)
#     if (length(weat_breaks) > 6) {
#       weat_breaks <- seq(weather_min, weather_max, by = 20)
#     }
#     
#     hosp_breaks <- seq(hosp_bottom, band_bottom, by = 100)
#     if (length(hosp_breaks) > 6) {
#       hosp_breaks <- seq(hosp_bottom, band_bottom, by = 200)
#     }
#     
#     
#     ggplot(df_weekly, aes(x = date)) +
#       geom_col(aes(y = hosp_month), fill = "grey70", alpha = 0.7, width = 7) +
#       geom_smooth(aes(y = hosp_month, color = "Smoothed trend"), method = "loess", span = 0.07,
#                   se = FALSE, linewidth = 0.7) +
#       geom_line(aes(y = transform_weather(tp_month), color = "Total Precipitation (mm)"), linewidth = 0.7) +
#       geom_line(aes(y = transform_weather(mn2t_month), color = "Min Temp (°C)"), linewidth = 0.7) +
#       geom_line(aes(y = transform_weather(mx2t_month), color = "Max Temp (°C)"), linewidth = 0.7) +
#       scale_x_date(breaks = valid_breaks, date_labels = "%b %Y") +
#       scale_y_continuous(
#         name = "Total patients",
#         breaks = hosp_breaks,
#         sec.axis = sec_axis(~ untransform_weather(.), name = "Weather (°C/mm)",
#                             breaks=weat_breaks)
#       ) +
#       coord_cartesian(ylim = c(hosp_bottom, hosp_top)) +
#       scale_color_manual(values = c("Min Temp (°C)" = '#1D98EB',
#                                     "Max Temp (°C)" = '#D74B3C',
#                                     "Total Precipitation (mm)" = "#8968CD",
#                                     "Smoothed trend" = "black"),
#                          breaks = c("Max Temp (°C)", "Min Temp (°C)",
#                                     "Total Precipitation (mm)", "Smoothed trend")) +
#       labs(title=title, color = "") +
#       theme_classic() +
#       theme(legend.position = "bottom",
#             plot.title = element_text(hjust=0.5),
#             axis.text.x = element_text(angle = 45, hjust = 1),
#             axis.title.x = element_blank(),
#             axis.title.y = element_text(hjust = hjust.left),
#             axis.title.y.right = element_text(hjust = hjust.right),
#             panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
#             panel.grid.major.y = element_blank(),
#             panel.grid.minor = element_blank())
#   }
#   
#   df_weekly <- df %>%
#     mutate(date = floor_date(date, "week")) %>%
#     group_by(date) %>%
#     summarise(
#       mn2t_month = mean(mn2t, na.rm = TRUE),
#       mx2t_month = mean(mx2t, na.rm = TRUE),
#       tp_month = mean(tp, na.rm = TRUE),
#       hosp_month = mean(total_patients, na.rm = TRUE)
#     )
#   
#   (fig <- get_time_series(df_weekly))
#   
#   ggsave(
#     paste0(path_report, 'II. dlnm/SF_pooled_time_series.png'),
#     fig, width = 8, height = 4, dpi = 300, bg = 'white'
#   )
#   
#   # args <- list(
#   #   Una = list(
#   #     hosp_bottom = 200, band_bottom = 725, weather_min = 0, weather_max = 40
#   #   ),
#   #   Kangra = list(
#   #     hosp_bottom = 300, band_bottom = 800, weather_min = 0, weather_max = 40
#   #   ),
#   #   Nahan = list(
#   #     hosp_bottom = 200, band_bottom = 750, weather_min = 0, weather_max = 40
#   #   ),
#   #   Hamirpur = list(
#   #     hosp_bottom = 500, band_bottom = 1050, weather_min = 0, weather_max = 40
#   #   ),
#   #   Mandi = list(
#   #     hosp_bottom = 400, band_bottom = 1400, weather_min = 0, weather_max = 40
#   #   ),
#   #   Chamba = list(
#   #     hosp_bottom = 400, band_bottom = 900, weather_min = -10, weather_max = 35
#   #   ),
#   #   Dharamshala = list(
#   #     hosp_bottom = 100, band_bottom = 1100, weather_min = 0, weather_max = 80
#   #   ),
#   #   Solan = list(
#   #     hosp_bottom = 500, band_bottom = 1000, weather_min = -10, weather_max = 35
#   #   ),
#   #   Kullu = list(
#   #     hosp_bottom = 300, band_bottom = 800, weather_min = -10, weather_max = 35
#   #   ),
#   #   'Reckong Peo' = list(
#   #     hosp_bottom = 100, band_bottom = 600, weather_min = -10, weather_max = 35
#   #   )
#   # )
#   #
#   # city_ts <- list()
#   # for (l in locations){
#   #
#   #   df_weekly <- df %>%
#   #     filter(location==l) %>%
#   #     mutate(date = floor_date(date, "week")) %>%
#   #     group_by(date) %>%
#   #     summarise(
#   #       mn2t_month = mean(mn2t, na.rm = TRUE),
#   #       mx2t_month = mean(mx2t, na.rm = TRUE),
#   #       tp_month = mean(tp, na.rm = TRUE),
#   #       hosp_month = mean(total_patients, na.rm = TRUE)
#   #     )
#   #
#   #   city_ts[[l]] <- get_time_series(df_weekly, title=l, month_breaks='6 months',
#   #                                   hjust.left=0.5, hjust.right=0.5,
#   #                                   hosp_bottom=args[[l]]$hosp_bottom,
#   #                                   band_bottom=args[[l]]$band_bottom,
#   #                                   weather_min=args[[l]]$weather_min,
#   #                                   weather_max=args[[l]]$weather_max)
#   #
#   # }
#   # fig <- wrap_plots(city_ts, nrow = 2) +
#   #   plot_layout(guides = "collect", axis_titles  = "collect") &
#   #   theme(legend.position = "bottom")
#   #
#   # ggsave(
#   #   paste0(path_report, 'II. dlnm/SF_city_time_series.png'),
#   #   fig, width = 15, height = 7, dpi = 300, bg = 'white'
#   # )
#   
#   
# }
# 
# Subgroup exposure-response ----

if(export == 't2m'){
  
  get_pooled_df <- function(meta, reference = '50%') {
    list2env(meta, envir = environment())
    
    if (!is.null(meta$mv)) {
      coef_use  <- coef(mv)
      vcov_use  <- vcov(mv)
      cen_clean <- as.data.frame(cen[!apply(cen, 1, anyDuplicated), ])
    } else {
      city_name <- names(Slist)
      coef_use  <- as.numeric(ymat[city_name, ])
      vcov_use  <- Slist[[city_name]]
      cen_clean <- cen[!duplicated(cen) & !duplicated(cen, fromLast = TRUE), , drop = FALSE]
    }
    
    perc    <- as.numeric(sub("%", "", names(rowMeans(cen_clean))))
    ref_idx <- if (reference == '0%') paste0(min(perc), '%') else reference
    cp_at   <- rowMeans(cen_clean)
    cp      <- crosspred(b_mean, coef = coef_use, vcov = vcov_use, model.link = 'log',
                         at = cp_at, cen = cp_at[ref_idx])
    
    data.frame(pred = perc, RR = cp$allRRfit, low = cp$allRRlow, high = cp$allRRhigh)
  }
  
  p_temp_imm   <- make_subgroup_plot(data=temp_imm,   x_lab="Temperature (%)",   pal=pal)
  make_subgroup_plot <- function(data, x_lab, pal, xlim = NULL, ylim = NULL, title=NULL) {
    p <- ggplot(data, aes(x = pred, y = RR, colour = Subgroup, fill = Subgroup, group = Subgroup)) +
      geom_hline(yintercept = 1, color = "gray40") +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1, colour = NA) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(name = NULL, values = pal) +
      scale_fill_manual(name = NULL, values = pal) +
      labs(title=title, x = x_lab, y = "RR") +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    
    if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim, breaks = pretty(xlim, 5))
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, xlim = range(data$pred, na.rm = TRUE)) +
        scale_y_continuous(breaks = pretty(ylim, 5))
    p
  }
  get_curves <- function(output, outcomes){
    
    all_df <- map_dfr(outcomes, function(i) {
      label <- as.factor(replace_exact_matches(i, nm))
      bind_rows(
        get_pooled_df(output$tmp[[i]], '50%') |> mutate(Subgroup = label, exposure = "Temperature"),
        get_pooled_df(output$tp[[i]],  tp_ref) |> mutate(Subgroup = label, exposure = "Precipitation")
      )
    }) |>
      mutate(pred = as.numeric(pred)) #|>
    # filter((exposure == "Temperature"   & pred <= 20) |
    #          (exposure == "Precipitation" & pred >= 80))
    
    return(all_df)
  }
  
  # Load model outputs
  sub1 <- readRDS(paste0('2_DLNM/Results/',  export, '_', model1, '_subgroup_plots.RData'))
  sub2 <- readRDS(paste0('2_DLNM/Results/',  export, '_', model2, '_subgroup_plots.RData'))
  outcomes <- h_cols[c(2,4:length(h_cols))]
  
  sub1_df <- get_curves(sub1$pooled_meta, outcomes)
  sub2_df <- get_curves(sub2$pooled_meta, outcomes)
  
  
  gender_groups <- c("Women", "Men")
  age_groups    <- c("Under 5", "Ages 5-14", "Ages 65+")
  opd_groups    <- c('General', 'Ophthalmology', 'ENT', 'Gynaecology')
  gender_pal <- mypalette[c(2,1)] #c(scales::hue_pal()(2))
  age_pal    <- mypalette[c(3,6,7)] #c(scales::hue_pal()(5))
  opd_pal    <- mypalette[c(5,11)] #c(scales::hue_pal()(7))
  
  make_panel <- function(groups=gender_groups, pal=gender_pal, title=TRUE) {
    
    imm_t   <- filter(sub1_df, exposure == 'Temperature',   Subgroup %in% groups) |> droplevels()
    cum_t  <- filter(sub2_df, exposure == 'Temperature',   Subgroup %in% groups) |> droplevels()
    imm_p   <- filter(sub1_df, exposure == 'Precipitation',   Subgroup %in% groups) |> droplevels()
    cum_p  <- filter(sub2_df, exposure == 'Precipitation',   Subgroup %in% groups) |> droplevels()
    
    p_imm_t   <- make_subgroup_plot(imm_t, "Temperature (%)",   pal)
    p_cum_t   <- make_subgroup_plot(cum_t, "Temperature (%)",   pal)
    p_imm_p   <- make_subgroup_plot(imm_p, "Precipitation (%)",   pal)
    p_cum_p   <- make_subgroup_plot(cum_p, "Precipitation (%)",   pal)
    
    if(title==TRUE){
      wrap_plots(
        wrap_plots(title_only('Daily'), title_only('Seven-day')),
        wrap_plots(list(
          p_imm_t + theme(legend.position = "none"),
          p_imm_p + theme(legend.position = "none"),
          plot_spacer(),
          p_cum_t + theme(legend.position = "none"),
          p_cum_p + theme(legend.position = "none")), nrow = 1, widths = c(1,1,0.2,1,1)),
        get_legend(p_cum_p),
        ncol = 1, heights = c(0.1, 1, 0.1)
      )
      
    } else {
      wrap_plots(
        wrap_plots(list(
          p_imm_t + theme(legend.position = "none"),
          p_imm_p + theme(legend.position = "none"),
          plot_spacer(),
          p_cum_t + theme(legend.position = "none"),
          p_cum_p + theme(legend.position = "none")), nrow = 1, widths = c(1,1,0.2,1,1)),
        get_legend(p_cum_p),
        ncol = 1, heights = c(1, 0.1)
      )
    }
  }
  
  panel_gender <- make_panel(gender_groups, gender_pal, title=TRUE)
  panel_age   <- make_panel(age_groups,    age_pal, title=FALSE)
  panel_opd    <- make_panel(opd_groups,    opd_pal, title=FALSE)
  
  final_plot <- wrap_plots(
    list(row_label("(A)"), wrap_plots(panel_gender, nrow = 1),
         row_label("(B)"), wrap_plots(panel_age, nrow = 1),
         row_label("(C)"), wrap_plots(panel_opd,  nrow = 1)),
    ncol = 2, widths = c(0.05, 1)
  )
  
  ggsave(paste0(path_report, 'II. dlnm/SF_subgroups.png'), final_plot, width = 10, height = 8, dpi = 300, bg = 'white')
  
}

#### TABLES ####
# ST: Descriptive statistics ----

get_desc <- function(df, h_cols, w_cols, p_cols=NULL) {
  
  h_cols <- h_cols[h_cols != "location"]
  check  <- c('F', 'under_5', 'general', 't2m', 'elev')
  labels <- c('Gender', 'Age', 'Ward', 'Weather Characteristics', 'Socio-geographic Characteristics')
  
  locs <- unique(df$location)
  
  summarise_col <- function(col, sub_df) {
    x         <- sub_df[[col]]
    is_binary <- all(x %in% c(0, 1, NA), na.rm = TRUE)
    
    total_visits <- sum(x, na.rm = TRUE)
    raw_mean <- if (is_binary) NA else mean(x, na.rm = TRUE)
    mean_x <- if (is_binary) NA else round(raw_mean)
    sd_x   <- if (is_binary) NA else round(sd(x,   na.rm = TRUE))
    min_x  <- if (is_binary) NA else round(min(x,  na.rm = TRUE))
    max_x  <- if (is_binary) NA else round(max(x,  na.rm = TRUE))
    
    avg <- if (is_binary || is.na(sd_x)) {
      NA
    } else if (col %in% p_cols) {
      if (!is.na(raw_mean) && raw_mean >= 0 && raw_mean <= 1) {
        as.character(round(raw_mean, 2) * 100)
      } else {
        as.character(mean_x)
      }
    } else {
      paste0(mean_x, ifelse(col == 'tp', ' + ', ' ± '), sd_x, ' (', min_x, ', ', max_x, ')')
    }
    
    data.frame(
      Variable     = sub('_b$', '', col),
      total_visits = total_visits,
      Average      = avg
    )
  }
  
  all_cols <- c(h_cols, w_cols, p_cols)
  
  overall <- bind_rows(lapply(all_cols, summarise_col, sub_df = df)) |>
    dplyr::select(Variable, total_visits, Overall = Average)
  
  per_loc <- bind_rows(lapply(locs, function(loc) {
    bind_rows(lapply(all_cols, summarise_col, sub_df = df[df$location == loc, ])) |>
      mutate(location = loc)
  })) |>
    dplyr::select(Variable, location, Average) |>
    tidyr::pivot_wider(names_from = location, values_from = Average)
  
  base_desc <- overall |>
    left_join(per_loc, by = "Variable") |>
    dplyr::select(Variable, total_visits, Overall, all_of(as.character(locs)))
  
  # Sort ward rows by descending total visits
  ward_cols <- h_cols[11:28]
  ward_idx  <- which(base_desc$Variable %in% ward_cols)
  if (length(ward_idx) > 0) {
    base_desc[ward_idx, ] <- base_desc[ward_idx, ][
      order(base_desc$total_visits[ward_idx], decreasing = TRUE, na.last = TRUE), ]
  }
  
  base_desc <- base_desc |> dplyr::select(-total_visits)
  
  # Insert section headings
  blank_row <- function(v) setNames(
    as.list(c(v, rep(NA, ncol(base_desc) - 1))), names(base_desc)
  ) |> as_tibble()
  
  base_desc <- bind_rows(blank_row('Hospital Characteristics'), base_desc)
  positions <- which(base_desc$Variable %in% check)
  new_rows  <- bind_rows(lapply(labels, blank_row))
  for (i in rev(seq_along(positions))) {
    pos       <- positions[i]
    base_desc <- bind_rows(base_desc[1:(pos - 1), ], new_rows[i, ], base_desc[pos:nrow(base_desc), ])
  }
  
  base_desc
}

# replace nm % with //%

base <- get_desc(df_wide, h_cols[-1], w_cols, p_cols) |> replace_exact_matches(nm)
value_cols   <- setdiff(names(base), "Variable")
header_rows  <- which(grepl("Characteristics", base$Variable))
section_rows <- c(setdiff(which(rowSums(!is.na(base[value_cols])) == 0), header_rows),
                  2, 3, header_rows[2]:nrow(base))

base_display <- base |>
  mutate(
    across(all_of(value_cols), trimws),
    across(all_of(value_cols), ~gsub("%", "\\%", .x, fixed = TRUE)),
    across(everything(), escape_latex_unicode),
    across(everything(), ~ifelse(is.na(.), "", .)),
    Variable  = case_when(
      row_number() %in% header_rows  ~ paste0("\\textbf{", Variable, "}"),
      row_number() %in% section_rows ~ paste0("\\hspace{1em}", Variable),
      TRUE                           ~ paste0("\\hspace{2em}\\textit{", Variable, "}")
    )
  )

# split columns into two groups, e.g. first half vs second half
split_pt <- ceiling(length(value_cols) / 2)
cols1 <- value_cols[1:split_pt]
cols2 <- value_cols[(split_pt + 1):length(value_cols)]

build_tab <- function(disp, cols, title, label) {
  t <- get_kbl_tab(
    disp[, c("Variable", cols)],
    title    = NULL,
    note     = "Averages are presented as daily total registrations.
    Precipitation values are reported with positive SDs only because of its
    inverse exponential distribution.",
    colnames = c("Variable", cols),
    label    = label
  ) |>
    column_spec(1, width = "4.5cm", latex_column_spec = "p{4.5cm}")
  
  for (i in seq_along(cols)) {
    t <- t |> column_spec(i + 1, width = "2.5cm", latex_column_spec = "p{2.5cm}")
  }
  t
}

tab1 <- build_tab(base_display, cols1, "Descriptive statistics of hospital and weather characteristics", "MT_base_desc_t1")
tab2 <- build_tab(base_display, cols2, "Descriptive statistics of hospital and weather characteristics", "MT_base_desc_t2")

writeLines(as.character(tab1), paste0(path_report, "II. dlnm/ST_base_desc_t1.tex"))
writeLines(as.character(tab2), paste0(path_report, "II. dlnm/ST_base_desc_t2.tex"))


# # ST: Meta-regression results ----
# 
# prep_display <- function(df) {
#   data.frame(
#     `Meta-predictors` = df[["Meta-predictors"]],
#     `Stat`   = fmt_num(df$LR_stat), `df`  = fmt_int(df$LR_df),  `p`    = fmt_p(df$LR_p),
#     `AIC`    = fmt_num(df$AIC),
#     `Stat `  = fmt_num(df$Q_stat),  `df ` = fmt_int(df$Q_df),   `p `   = fmt_p(df$Q_p),
#     `I2 (%)` = fmt_num(df$I2),
#     check.names = FALSE, stringsAsFactors = FALSE
#   )
# }
# 
# blank_row <- data.frame(
#   `Meta-predictors` = "", `Stat` = "", `df` = "", `p` = "", `AIC` = "",
#   `Stat ` = "", `df ` = "", `p ` = "", `I2 (%)` = "",
#   check.names = FALSE, stringsAsFactors = FALSE
# )
# 
# output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_meta_tables.RData'))
# d_tmp   <- prep_display(output$tbl_tmp)
# d_tp    <- prep_display(output$tbl_tp)
# immediate <- rbind(blank_row, d_tmp, blank_row, d_tp)
# immediate[1, 1] <- paste0(str_to_title(model1), ' Temperature')
# immediate[19, 1] <- paste0(str_to_title(model1), ' Precipitation')
# immediate <- immediate |>
#   mutate(`Meta-predictors` = ifelse(
#     `Meta-predictors` %in% c(paste0(str_to_title(model1), ' Temperature'), 
#                              paste0(str_to_title(model1), ' Precipitation'), ''),
#     `Meta-predictors`,
#     paste0("\\hspace{1em}", `Meta-predictors`)
#   ))
# 
# 
# output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_meta_tables.RData'))
# d_tmp   <- prep_display(output$tbl_tmp)
# d_tp    <- prep_display(output$tbl_tp)
# cumulative <- rbind(blank_row, d_tmp, blank_row, d_tp)
# cumulative[1, 1] <- paste0(str_to_title(model2), ' Temperature')
# cumulative[19, 1] <- paste0(str_to_title(model2), ' Precipitation')
# cumulative <- cumulative |>
#   mutate(`Meta-predictors` = ifelse(
#     `Meta-predictors` %in% c(paste0(str_to_title(model2), ' Temperature'), 
#                              paste0(str_to_title(model2), ' Precipitation'), ''),
#     `Meta-predictors`,
#     paste0("\\hspace{1em}", `Meta-predictors`)
#   ))
# 
# combined <- rbind(immediate, cumulative)
# 
# tab <- get_kbl_tab(
#   combined,
#   title           = "Meta-regression results for immediate and cumulative weather models",
#   note            = "LR test was used to test the significance of meta-predictors with the intercept-only model as reference. Cochran Q test was used to test the significance of residual heterogeneity with the null hypothesis as no heterogeneity. LR = likelihood ratio; AIC = Akaike information criterion; $I^2$ = heterogeneity statistic.",
#   colnames        = c("Meta-predictors", "Stat", "df", "\\textit{p}", "AIC",
#                       "Stat", "df", "\\textit{p}", "(\\%)"),
#   label           = "MT_meta_reg",
#   spanning_header = c(" " = 1, "LR test" = 3, "Model fits" = 1, "Cochran Q test" = 3, "$I^2$" = 1),
#   bold_rows       = c(1, 19, 37, 55),
#   italic_rows     = c(2, 20, 38, 56),
#   hline_after_rows = c(18, 36, 54)
# )
# 
# tab <- sub(
#   "\\\\multicolumn\\{9\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}",
#   "\\\\multicolumn{9}{p{0.8\\\\textwidth}}{\\\\rule{0pt}{1em}",
#   tab
# )
# 
# writeLines(as.character(tab), paste0(path_report, "II. dlnm/ST_meta_reg_", export, ".tex"))
# 
