### ------------------------------------------------------------------------ ### 
###  DNLM – Tables & Figures
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
contour_legend  <- function(ref_plot, bar_w = 5, bar_h = 0.3, txt = 12, ttl = 13) {
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
#### MAIN CODE ####

export <-a

# Load data
load(paste0(path_data, '2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata'))
df_wide <- df_wide |>
  filter(!location %in% c('Bilaspur'),
         date >= as.Date('2015-02-14'),
         date <= as.Date('2020-02-14')) |>
  rename(rh = rh_mean)

h_cols <- c('closed', 'total_patients', 'location', 'F', 'M', 'under_5', 'age5_64', 'over_65') 
w_cols <- c("t2m","mn2t", "mx2t", "tp","fg10", "rh",   "pm25") 
p_cols <- c('elev', 'popu', 'hosp_n', 'insurance', 'decision', 'meducation', 'hsize', 'hwealth')

#### FIGURES ####

results1 <- readRDS(paste0('outputs/', export, '_', model1, '_pooled_plots.RData'))
results2 <- readRDS(paste0('outputs/',  export, '_', model2, '_pooled_plots.RData'))

# MF/SF: Time series visual ####

if(export == 't2m'){

  get_time_series <- function(df_weekly, title=NULL, month_breaks='4 months', hjust.left=0.25, hjust.right=0.2,
                              hosp_bottom=400, band_bottom=800, weather_min=0, weather_max=30){

    # Weather occupies the band from band_bottom (i.e. anchor) up to hosp_top
    hosp_top <- round(band_bottom + (band_bottom/2.667))
    band_top <- hosp_top

    # Set breaks
    hosp_breaks <- seq(hosp_bottom, hosp_top, by = 100)
    hosp_labels <- ifelse(hosp_breaks > hosp_top, "", hosp_breaks)

    transform_weather <- function(x) {
      band_bottom + (x - weather_min) / (weather_max - weather_min) * (band_top - band_bottom)
    }
    untransform_weather <- function(y) {
      (y - band_bottom) / (band_top - band_bottom) * (weather_max - weather_min) + weather_min
    }

    # Generate breaks every 4 months, then drop anything past the last real data point
    max_date <- max(df_weekly$date, na.rm = TRUE)
    all_breaks <- seq(min(df_weekly$date, na.rm = TRUE), max_date + months(4), by = month_breaks)
    valid_breaks <- all_breaks[all_breaks <= max_date]

    weat_breaks <- seq(weather_min, weather_max, by = 10)
    if (length(weat_breaks) > 6) {
      weat_breaks <- seq(weather_min, weather_max, by = 20)
    }

    hosp_breaks <- seq(hosp_bottom, band_bottom, by = 100)
    if (length(hosp_breaks) > 6) {
      hosp_breaks <- seq(hosp_bottom, band_bottom, by = 200)
    }


    ggplot(df_weekly, aes(x = date)) +
      geom_col(aes(y = hosp_week), fill = "grey70", alpha = 0.7, width = 7) +
      geom_smooth(aes(y = hosp_week, color = "Smoothed trend"), method = "loess", span = 0.07,
                  se = FALSE, linewidth = 0.7) +
      geom_line(aes(y = transform_weather(tp_week), color = "Total Precipitation (mm)"), linewidth = 0.7) +
      geom_line(aes(y = transform_weather(t2m_week), color = "Mean Temperature (°C)"), linewidth = 0.7) +
      # geom_line(aes(y = transform_weather(mx2t_week), color = "Max Temp (°C)"), linewidth = 0.7) +
      scale_x_date(breaks = valid_breaks, date_labels = "%b %Y") +
      scale_y_continuous(
        name = "Total patients",
        breaks = hosp_breaks,
        sec.axis = sec_axis(~ untransform_weather(.), name = "Weather variables",
                            breaks=weat_breaks)
      ) +
      coord_cartesian(ylim = c(hosp_bottom, hosp_top)) +
      scale_color_manual(values = c("Mean Temperature (°C)" = '#D74B3C',
                                    "Total Precipitation (mm)" = "#1D98EB",
                                    "Smoothed trend" = "black"),
                         breaks = c("Mean Temperature (°C)", "Total Precipitation (mm)", "Smoothed trend")) +
      labs(title=title, color = "") +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust=0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(hjust = hjust.left, size=12),
            axis.title.y.right = element_text(hjust = hjust.right-0.1, size=12),
            panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank())
  }

  df_weekly <- df %>%
    mutate(date = floor_date(date, "week")) %>%
    group_by(date) %>%
    summarise(
      t2m_week = mean(t2m, na.rm = TRUE),
      # mx2t_week = mean(mx2t, na.rm = TRUE),
      tp_week = mean(tp, na.rm = TRUE),
      hosp_week = mean(total_patients, na.rm = TRUE)
    )

  (fig <- get_time_series(df_weekly))

  ggsave(
    paste0(path_report, 'MF2.png'),
    fig, width = 8, height = 4, dpi = 300, bg = 'white'
  )

  # args <- list(
  #   Una = list(
  #     hosp_bottom = 200, band_bottom = 725, weather_min = 0, weather_max = 40
  #   ),
  #   Kangra = list(
  #     hosp_bottom = 300, band_bottom = 800, weather_min = 0, weather_max = 40
  #   ),
  #   Nahan = list(
  #     hosp_bottom = 200, band_bottom = 750, weather_min = 0, weather_max = 40
  #   ),
  #   Hamirpur = list(
  #     hosp_bottom = 500, band_bottom = 1050, weather_min = 0, weather_max = 40
  #   ),
  #   Mandi = list(
  #     hosp_bottom = 400, band_bottom = 1400, weather_min = 0, weather_max = 40
  #   ),
  #   Chamba = list(
  #     hosp_bottom = 400, band_bottom = 900, weather_min = -10, weather_max = 35
  #   ),
  #   Dharamshala = list(
  #     hosp_bottom = 100, band_bottom = 1100, weather_min = 0, weather_max = 80
  #   ),
  #   Solan = list(
  #     hosp_bottom = 500, band_bottom = 1000, weather_min = -10, weather_max = 35
  #   ),
  #   Kullu = list(
  #     hosp_bottom = 300, band_bottom = 800, weather_min = -10, weather_max = 35
  #   ),
  #   'Reckong Peo' = list(
  #     hosp_bottom = 100, band_bottom = 600, weather_min = -10, weather_max = 35
  #   )
  # )
  #
  # city_ts <- list()
  # for (l in locations){
  #
  #   df_weekly <- df %>%
  #     filter(location==l) %>%
  #     mutate(date = floor_date(date, "week")) %>%
  #     group_by(date) %>%
  #     summarise(
  #       mn2t_week = mean(mn2t, na.rm = TRUE),
  #       mx2t_week = mean(mx2t, na.rm = TRUE),
  #       tp_week = mean(tp, na.rm = TRUE),
  #       hosp_week = mean(total_patients, na.rm = TRUE)
  #     )
  #
  #   city_ts[[l]] <- get_time_series(df_weekly, title=l, month_breaks='6 months',
  #                                   hjust.left=0.5, hjust.right=0.5,
  #                                   hosp_bottom=args[[l]]$hosp_bottom,
  #                                   band_bottom=args[[l]]$band_bottom,
  #                                   weather_min=args[[l]]$weather_min,
  #                                   weather_max=args[[l]]$weather_max)
  #
  # }
  # fig <- wrap_plots(city_ts, nrow = 2) +
  #   plot_layout(guides = "collect", axis_titles  = "collect") &
  #   theme(legend.position = "bottom")
  #
  # ggsave(
  #   paste0(path_report, 'supplement/SF_city_time_series.png'),
  #   fig, width = 15, height = 7, dpi = 300, bg = 'white'
  # )


}

# MF/SF: All exposure-lag-response results (t2m, tp, mn2t, mx2t) ----


if(export == 't2m'){
  
  make_exp_fig <- function(cont1, cont2, exp1, exp2, width=c(1,3), center=FALSE) {
    
    # Daily weather
    imm_con <- strip_legend(cont1$ALL)
    imm_cum <- strip_legend(exp1)
    imm_lagP5 <- strip_legend(lag1$P5$p_blup)
    imm_lagP95 <- strip_legend(lag1$P95$p_blup)
    
    # Seven-day weather
    cum_con <- strip_legend(cont2$ALL)
    cum_cum <- strip_legend(exp2)
    cum_lagP5 <- strip_legend(lag2$P5$p_blup)
    cum_lagP95 <- strip_legend(lag2$P95$p_blup)
    
    if (center == TRUE) {
      
      plots <- list(imm_cum, imm_lagP95, cum_cum, cum_lagP95)
      
      # Combined y range across all plots
      y_rng <- range(unlist(lapply(plots, function(p) {
        ggplot_build(p)$layout$panel_params[[1]]$y.range
      })))
      
      # Shared pretty breaks, and snap the limits to them
      y_brks <- pretty(y_rng, n = 5)
      y_lims <- range(y_brks)
      
      # Apply the same scale to every plot
      fix_y <- function(p) {
        p +
          scale_y_continuous(
            breaks = y_brks,
            labels = scales::label_number(accuracy = 0.01)  # same decimals everywhere
          ) +
          coord_cartesian(ylim = y_lims, expand = FALSE)
      }
      
      imm_cum    <- fix_y(imm_cum)
      imm_lagP95 <- fix_y(imm_lagP95)
      cum_cum    <- fix_y(cum_cum)
      cum_lagP95 <- fix_y(cum_lagP95)
    }

    # Get legend
    legends <- wrap_plots(
      list(#plot_spacer(),
           contour_legend(cont2$ALL),
           get_legend(exp1 + theme(legend.position = "bottom"))),
      nrow = 1, widths = width
    )
    
    # Set abc labels
    plots <- c(imm_con, imm_cum, imm_lagP5, imm_lagP95,
               cum_con, cum_cum, cum_lagP5, cum_lagP95)
    labels <- paste0("(", LETTERS[seq_along(plots)], ") ")
    plots <- Map(function(p, lab) {
      p + labs(title = paste0(lab, p$labels$title))
    }, plots, labels)
    
    # Plot
    all_rows <- wrap_plots(plots, nrow = 2)
    wrap_plots(list(all_rows, legends), ncol = 1, heights = c(1, 0.15))
    
  }
  
  cont1 <- results1$cont$tmp; exp1 <- results1$exp$tmp; lag1 <- results1$lag$tmp
  cont2 <- results2$cont$tmp; exp2 <- results2$exp$tmp; lag2 <- results2$lag$tmp
  fig <- make_exp_fig(cont1, cont2, exp1, exp2, width=c(1,3))
  ggsave(paste0(path_report, 'MF3.png'), fig, width = 14, height = 7, dpi = 300, bg = 'white')

  cont1 <- results1$cont$tp; exp1 <- results1$exp$tp; lag1 <- results1$lag$tp
  cont2 <- results2$cont$tp; exp2 <- results2$exp$tp; lag2 <- results2$lag$tp
  fig <- make_exp_fig(cont1, cont2, exp1, exp2, width=c(1, 2), center=TRUE)
  ggsave(paste0(path_report, 'MF4.png'), fig, width = 10, height = 7, dpi = 300, bg = 'white')

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
  output  <- readRDS(paste0('outputs/', export, '_', model1, '_meta_tables.RData'))
  meta_imm_tmp <- get_meta_plot(prep_display(output$tbl_tmp))
  meta_imm_tp <- get_meta_plot(prep_display(output$tbl_tp))
  output  <- readRDS(paste0('outputs/', export, '_', model2, '_meta_tables.RData'))
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
    paste0(path_report, 'supplement/S7.png'),
    fig, width = 20, height = 22, dpi = 100, bg = 'white'
  )

} else if (export %in% c('mx2t', 'mn2t')) {

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
  output  <- readRDS(paste0('outputs/', export, '_', model1, '_meta_tables.RData'))
  imm_tmp   <- prep_display(output$tbl_tmp)
  meta_imm <- get_meta_plot(imm_tmp)
  output  <- readRDS(paste0('outputs/', export, '_', model2, '_meta_tables.RData'))
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
    paste0(path_report, 'supplement/SF_', export, '_res.png'),
    fig, width = 20, height = 20, dpi = 100, bg = 'white'
  )

}

# MF: Attributable estimations ####

if (export=='t2m'){
  
  attr1 <- readRDS(paste0('outputs/', export, '_', model1, '_attr.RData'))
  attr2 <- readRDS(paste0('outputs/',  export, '_', model2, '_attr.RData'))
  plot_df <- rbind(attr1, attr2)
  
  # Set legend colours
  fill_values <- unlist(lapply(levels(plot_df$panel), function(p) {
    panel_colors <- c(cold = "#1D98EB", heat = "#D74B3C", precip = "#8968CD")
    base_col <- function(panel) {
      hit <- names(panel_colors)[sapply(names(panel_colors), function(k) grepl(k, panel, ignore.case = TRUE))]
      if (length(hit) == 0) "grey50" else panel_colors[[hit[1]]]
    }
    setNames(c(base_col(p), "grey50"), paste(p, c("FALSE", "TRUE"), sep = "."))
  }))
  get_color <- function(pattern) {
    hit <- grep(paste0("^", pattern, "\\."), names(fill_values), value = TRUE)
    hit <- hit[!grepl("TRUE$", hit)]  # exclude the "Total" (TRUE) entries
    unname(fill_values[hit[1]])
  }
  legend_colors <- c(
    Cold = get_color(".*Cold"),
    Heat = get_color(".*Heat"),
    Precipitation = get_color(".*Precipitation"),
    Total = "grey50"
  )
  
  panel_labels <- setNames(
    paste0("(", LETTERS[seq_along(levels(plot_df$panel))], ") ", levels(plot_df$panel)),
    levels(plot_df$panel))
  
  # in the ggplot call, replace the facet_wrap line with:
  
  fig <- ggplot(
    plot_df,
    aes(x = location, y = af*100, fill = interaction(panel, location == "TOTAL"))) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_errorbar(aes(ymin = af_low*100, ymax = af_high*100), width = 0.25) +
    facet_wrap(~panel, nrow = 2, labeller = labeller(panel = panel_labels)) +
    scale_fill_manual(values = fill_values) +
    new_scale_fill() +
    geom_point(
      data = data.frame(legend_group = names(legend_colors)),
      aes(x = 1, y = 0, fill = factor(legend_group, levels = names(legend_colors))),
      shape = 22, size = 0, inherit.aes = FALSE
    ) +
    scale_fill_manual(name = NULL, values = legend_colors,
                      guide = guide_legend(override.aes = list(size = 5))) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    labs(x = NULL, y = "Attributable fraction (%)") +
    # coord_cartesian(ylim = c(-0.06, 0.03)) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "pt")
    )
  
  ggsave(paste0(path_report, 'MF5.png'), fig, width = 8, height = 5, dpi = 300, bg = 'white')
  
}

# SF: Subgroup exposure-response ----

if(export == 't2m'){
 
  # Load model outputs
  sub1 <- readRDS(paste0('outputs/',  export, '_', model1, '_subgroup_plots.RData'))
  sub2 <- readRDS(paste0('outputs/',  export, '_', model2, '_subgroup_plots.RData'))
  
  # Get legend
  leg1 <- get_legend(sub1$res_plots$immediate_tmp_gen)
  leg2 <- get_legend(sub1$res_plots$immediate_tmp_age)
  plt_imm <- lapply(sub1$res_plots, function(p) p + theme(legend.position = "none"))
  plt_cum <- lapply(sub2$res_plots, function(p) p + theme(legend.position = "none"))
  
  fig <- wrap_plots(
    wrap_plots(
      list(
        wrap_plots(list(plt_imm[[1]] + labs(x = "Temperature Percentile (%)",       title = "(A) Daily temp: Gender"),
                        plt_cum[[1]] + labs(x = "Temperature Percentile (%)",   title = "(B) Seven-day temp: Gender"),
                        plt_imm[[3]] + labs(x = "Temperature Percentile (%)",       title = "(C) Daily temp: Age"),
                        plt_cum[[3]] + labs(x = "Temperature Percentile (%)",   title = "(D) Seven-day temp: Age")), nrow = 1),
        wrap_plots(list(plt_imm[[2]] + labs(x = "Precipitation Percentile (%)",       title = "(E) Daily precip: Gender"),
                        plt_cum[[2]] + labs(x = "Precipitation Percentile (%)",   title = "(F) Seven-day precip: Gender"),
                        plt_imm[[4]] + labs(x = "Precipitation Percentile (%)",       title = "(G) Daily precip: Age"),
                        plt_cum[[4]] + labs(x = "Precipitation Percentile (%)",   title = "(H) Seven-day precip: Age")), nrow = 1)
      ),
      nrow = 2
    ),
    wrap_plots(list(leg1, leg2), nrow = 1),
    ncol = 1, heights = c(1, 0.1)
  ) &
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
  
  ggsave(paste0(path_report, 'supplement/S6.png'), fig, width = 10, height = 5, dpi = 300, bg = 'white')
  
}

# SF: Multi-day windows ####

if (grepl('roll6', export)) {
  
  roll_files <- list.files('outputs/', pattern = "roll.*\\.rds$", full.names = TRUE)
  roll_files <- roll_files[grepl("roll[2-6](?!\\d)", basename(roll_files), perl = TRUE)]
  roll_models <- lapply(roll_files, readRDS)
  names(roll_models) <- sub("_cumulative_pooled_model\\.rds$", "", basename(roll_files))

  get_curves <- function(x='meta_tmp', xlab='', ref='50%'){
    
    curves <- bind_rows(lapply(names(roll_models), function(i) {
      cp  <- meta_predict(roll_models[[i]][[x]], city='ALL', reference=ref, blup=TRUE)
      idx <- match(cp$at, cp$cp$predvar)
      res <- data.frame(
        test = i,
        perc = cp$perc,
        tmp  = cp$at,
        rr   = cp$cp$allRRfit[idx],
        lo   = cp$cp$allRRlow[idx],
        hi   = cp$cp$allRRhigh[idx]
      )
      row.names(res) <- NULL
      return(res)
    }))
    curves <- curves %>% mutate(test = gsub("^roll(\\d+)$", "\\1 days", test))
    # labels <- labels[order(as.numeric(gsub("-days", "", )))]
    labels <- unique(curves$test)
    curves$test <- factor(curves$test, labels)
    cols <- setNames(scales::hue_pal()(length(labels)), labels)
    
    fig <- ggplot() +
      geom_ribbon(data = curves, aes(x = perc, ymin = lo, ymax = hi, fill = test),
                  alpha = 0.2, show.legend = FALSE) +
      geom_line(data = curves, aes(x = perc, y = rr, color = test, linetype = test),
                linewidth = 1.2) +
      geom_hline(yintercept = 1, color = 'gray20') +
      scale_color_manual(name = NULL, values = cols) +
      scale_fill_manual(name = NULL, values = cols) +
      scale_linetype_manual(name = NULL, values = rep('solid', length(labels))) +
      scale_x_continuous(breaks = pretty(curves$perc, n = 5)) +
      scale_y_continuous(breaks = pretty(c(min(curves$lo), max(curves$hi)), n = 5)) +
      labs(x = xlab, y = 'RR', title = '') +
      coord_cartesian(ylim = c(min(curves$lo), max(curves$hi)), xlim = range(curves$perc)) +
      guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
      theme_classic() +
      theme(
        axis.line = element_line(color = 'gray20'),
        axis.text = element_text(color = 'gray20'),
        axis.ticks = element_line(color = 'gray20'),
        plot.title = element_blank(),
        legend.position = 'bottom'
      )
  
    return(fig)
    
  }
  
  tmp_curv <- get_curves('meta_tmp', 'Temperature percentiles (%)')
  tp_curv <- get_curves(x='meta_tp', 'Percipitation percentiles (%)', ref=tp_ref)
  
  fig <- wrap_plots(
    list(row_label("(A)"), tmp_curv, row_label("(B)"), tp_curv),
    nrow = 1, widths = c(0.05, 1, 0.05, 1)
  ) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom')
  
  ggsave(paste0(path_report, 'supplement/S8.png'), fig, width = 8, height = 4, dpi = 300, bg = 'white')
  
}

# SF: Different temperature measures ####

if (grepl('mn2t', export)) {
  
  roll_files <- list.files('outputs/', pattern = "^(mn2t|mx2t|t2m).*_pooled_model\\.rds$", full.names = TRUE)
  roll_models <- lapply(roll_files, readRDS)
  names(roll_models) <- sub("_pooled_model\\.rds$", "", basename(roll_files))
  
  get_curves <- function(type = 'immediate', x = 'meta_tmp', xlab = '', ttl=NULL, ref = '50%'){
    
    var_labs <- c(mn2t = 'Minimum temperature', mx2t = 'Maximum temperature', t2m = 'Mean temperature')
    nms <- grep(paste0("_", type, "$"), names(roll_models), value = TRUE)
    
    curves <- bind_rows(lapply(nms, function(i) {
      cp  <- meta_predict(roll_models[[i]][[x]], city = 'ALL', reference = ref, blup = TRUE)
      idx <- match(cp$at, cp$cp$predvar)
      data.frame(
        test = var_labs[[sub(paste0("_", type, "$"), "", i)]],
        perc = cp$perc,
        tmp  = cp$at,
        rr   = cp$cp$allRRfit[idx],
        lo   = cp$cp$allRRlow[idx],
        hi   = cp$cp$allRRhigh[idx]
      )
    }))
    
    labels <- unname(var_labs)
    curves$test <- factor(curves$test, labels)
    cols <- setNames(scales::hue_pal()(length(labels)), labels)
    
    ggplot(curves, aes(x = perc)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = test), alpha = 0.2, show.legend = FALSE) +
      geom_line(aes(y = rr, color = test), linewidth = 1.2) +
      geom_hline(yintercept = 1, color = 'gray20') +
      scale_color_manual(name = NULL, values = cols) +
      scale_fill_manual(name = NULL, values = cols) +
      scale_x_continuous(breaks = pretty(curves$perc, n = 5)) +
      scale_y_continuous(breaks = pretty(c(min(curves$lo), max(curves$hi)), n = 5)) +
      labs(x = xlab, y = 'RR', title=ttl) +
      coord_cartesian(ylim = c(min(curves$lo), max(curves$hi)), xlim = range(curves$perc)) +
      guides(color = guide_legend(nrow = 1)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5, face='bold', size=13),
        axis.line = element_line(color = 'gray20'),
        axis.text = element_text(color = 'gray20'),
        axis.ticks = element_line(color = 'gray20'),
        
        legend.position = 'bottom'
      )
  }
  
  imm_curv <- get_curves('immediate', 'meta_tmp', 'Temperature percentiles (%)', ttl='(A) Daily effects')
  cum_curv  <- get_curves('cumulative', 'meta_tmp', 'Temperature percentiles (%)', ttl='(B) Seven-day effects')
  
  fig <- wrap_plots(list(imm_curv, cum_curv), nrow = 1) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom')
  
  ggsave(paste0(path_report, 'supplement/S9.png'), fig, width = 8, height = 4, dpi = 300, bg = 'white')

  }

#### TABLES ####

# ST: Descriptive statistics ----

if(export == 't2m'){
  
  get_desc <- function(df, h_cols, w_cols, p_cols=NULL) {
    
    h_cols <- h_cols[h_cols != "location"]
    check  <- c('F', 'under_5', 'general', 't2m', 'elev')
    labels <- c('Gender', 'Age', 'Weather Characteristics', 'Socio-geographic Characteristics')
    
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
        paste0(mean_x, ' (', min_x, ', ', max_x, ')')
      }
      
      data.frame(
        Variable     = sub('_b$', '', col),
        total_visits = total_visits,
        Average      = avg
      )
    }
    
    all_cols <- c(h_cols, w_cols, p_cols)
    
    overall <- bind_rows(lapply(all_cols, summarise_col, sub_df = df)) |>
      dplyr::select(Variable, total_visits, Total = Average)
    
    per_loc <- bind_rows(lapply(locs, function(loc) {
      bind_rows(lapply(all_cols, summarise_col, sub_df = df[df$location == loc, ])) |>
        mutate(location = loc)
    })) |>
      dplyr::select(Variable, location, Average) |>
      tidyr::pivot_wider(names_from = location, values_from = Average)
    
    base_desc <- overall |>
      left_join(per_loc, by = "Variable") |>
      dplyr::select(Variable, total_visits, Total, all_of(as.character(locs))) |> 
      dplyr::select(-total_visits)
    
    # ---- Missing total_patients row, one column per location + Total ----
    fmt_missing <- function(sub_df) {
      n_missing <- sum(is.na(sub_df$total_patients))
      n_total   <- nrow(sub_df)
      sprintf("%d (%.1f%%)", n_missing, n_missing / n_total * 100)
    }
    
    missing_row <- tibble(
      Variable = "Missing days",
      Total    = fmt_missing(df)
    ) |>
      bind_cols(
        setNames(as.list(sapply(locs, function(loc) fmt_missing(df[df$location == loc, ]))), as.character(locs))
      )
    
    # if (length(tp_idx) > 0) {
    #   base_desc <- bind_rows(missing_row, base_desc[1:nrow(base_desc), ])
    # }
    
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
  
  base <- get_desc(df_wide, h_cols[-1], w_cols, p_cols) |> replace_exact_matches(nm)
  value_cols   <- setdiff(names(base), "Variable")
  header_rows  <- which(grepl("Characteristics", base$Variable))
  section_rows <- c(setdiff(which(rowSums(!is.na(base[value_cols])) == 0), header_rows),
                    2, 3, header_rows[2]:nrow(base))
  
  base_display <- base |>
    mutate(
      across(all_of(value_cols), trimws),
      # across(all_of(value_cols), ~gsub("%", "\%", .x, fixed = TRUE)),
      across(everything(), escape_latex_unicode),
      across(everything(), ~ifelse(is.na(.), "", .)),
      Variable  = case_when(
        row_number() %in% header_rows  ~ paste0("\\textbf{", Variable, "}"),
        row_number() %in% section_rows ~ paste0("\\hspace{1em}", Variable),
        TRUE                           ~ paste0("\\hspace{2em}\\textit{", Variable, "}")
      )
    )
  
  # Split columns into two groups, e.g. first half vs second half
  split_pt <- ceiling(length(value_cols) / 2)
  cols1 <- value_cols[1:split_pt]
  cols2 <- value_cols[(split_pt + 1):length(value_cols)]
  
  build_tab <- function(disp, cols, title, label) {
    t <- get_kbl_tab(
      disp[, c("Variable", cols)],
      title    = NULL,
      note     =  NULL, #"Table cells present daily averages with minimum and maximum values in brackets, except missing days, which display the proportion in brackets. No ranges are provided for static socio-graphic characteristics.",
      colnames = c("Variable", cols),
      label    = label
    ) |>
      column_spec(1, width = "4.5cm", latex_column_spec = "p{4.5cm}")
    
    for (i in seq_along(cols)) {
      t <- t |> column_spec(i + 1, width = "2.25cm", latex_column_spec = "p{2.25cm}")
    }
    
    t <- gsub("%", "\\%", t, fixed = TRUE)
    
    # Force footnote/note rows to wrap instead of overflowing the page:
    # replace \multicolumn{n}{l}{...} with a p{} column sized to the
    # table's actual line width, so long notes wrap within the table.
    t <- gsub(
      "\\\\multicolumn\\{([0-9]+)\\}\\{l\\}",
      "\\\\multicolumn{\\1}{p{\\\\dimexpr\\\\linewidth-2\\\\tabcolsep\\\\relax}}",
      t
    )
    
    t
  }
  tab1 <- build_tab(base_display, cols1, "Descriptive statistics of hospital and weather characteristics", "MT_base_desc_t1")
  tab2 <- build_tab(base_display, cols2, "Descriptive statistics of hospital and weather characteristics", "MT_base_desc_t2")
  
  writeLines(as.character(tab1), paste0(path_report, "supplement/S4A.tex"))
  writeLines(as.character(tab2), paste0(path_report, "supplement/S4B.tex"))

}

# ST: Percentile lookup ####

if(export == 't2m'){
  
  datalist <- lapply(locations, function(region) df[df$location == region, ])
  
  # Percentiles to show (lower tail, then upper tail)
  lower_pct   <- c(0.01, 0.5, 1, 2, 3, 4, 5, 10, 15, 20, 25)
  upper_pct   <- c(75, 80, 85, 90, 95, 96, 97, 98, 99, 99.5, 99.99)
  percentiles <- c(lower_pct, upper_pct)
  
  # Clean labels: max 2 dp, no trailing zeros (e.g. "0.01", "25", "99.95")
  pct_label <- function(p) sub("\\.?0+$", "", formatC(p, format = "f", digits = 2))
  
  make_centile <- function(data, var, locations, probs, transform = identity, digits = 2) {
    tab <- sapply(data, function(x) round(transform(quantile(x[[var]], probs, na.rm = TRUE)), digits))
    tab <- setNames(as.data.frame(tab), locations)
    rownames(tab) <- paste0(pct_label(probs * 100), "%")
    tab
  }
  
  var_specs <- tibble::tribble(
    ~name,        ~col,            ~transform,
    "t2m_sing",   "t2m",           identity,
    "t2m_multi",  "t2m_roll7",     identity,
    "tp_sing",    "tp_log",        \(q) exp(q) - 2,
    "tp_multi",   "tp_roll7_log",  \(q) exp(q) - 10
  )
  
  cen <- var_specs |>
    pmap(\(name, col, transform)
         make_centile(datalist, col, locations, probs = percentiles / 100, transform = transform)) |>
    set_names(var_specs$name) |>
    map(rowMeans) |>
    as.data.frame()
  
  # Rename: split at last "_", replace variable part, relabel type part
  split_names <- str_split_fixed(names(cen), "_(?=[^_]+$)", 2)
  type_label  <- recode(split_names[, 2], sing = "Daily", multi = "Seven-day")
  names(cen)  <- paste0(type_label, " ", replace_exact_matches(split_names[, 1], nm))
  
  # Table data frame, all values to 2 dp
  tab_df <- cbind(Percentile = rownames(cen), cen)
  rownames(tab_df) <- NULL
  
  value_cols <- setdiff(names(tab_df), "Percentile")
  tab_df[value_cols] <- lapply(tab_df[value_cols], \(col) sprintf("%.2f", col))
  
  tab <- get_kbl_tab(tab_df, title = NULL,
                     note = 'Please note that daily and seven-day precipitation have different log transformations, which skew the zero-values.',
                     label = 'ST_perc',
                     colnames = names(tab_df),
                     hline_after_rows = length(lower_pct)) |>   # line between lower and upper tails
    gsub("%", "\\%", x = _, fixed = TRUE) |>
    gsub("\\midrule\\\\", "\\midrule", x = _, fixed = TRUE)
  
  writeLines(as.character(tab), paste0(path_report, "supplement/S5.tex"))
  
}
# ST: RR tables (Excel) ----

keep_pred <- function(pred) pred <= 5 | pred >= 95 |
  (pred >= 5  & pred <= 20 & pred %% 5  == 0) |
  (pred >= 80 & pred <= 95 & pred %% 5  == 0) |
  (pred >= 20 & pred <= 80 & pred %% 10 == 0)

fmt_rr <- function(rr, lo, hi) sprintf("%.2f (%.2f, %.2f)", rr, lo, hi)
is_sig <- function(lo, hi, rr) !(lo <= 1 & hi >= 1) & round(rr, 2) != 1.00

tot1 <- readRDS('outputs/t2m_immediate_pooled_model.rds')
tot2 <- readRDS('outputs/t2m_cumulative_pooled_model.rds')
sub1 <- readRDS('outputs/t2m_immediate_subgroup_plots.RData')
sub2 <- readRDS('outputs/t2m_cumulative_subgroup_plots.RData')

tp_ref_imm <- '55%'
tp_ref_cum <- '20%'

# One row per table: label, window (1 = daily, 2 = 7-day), exposure, reference
specs <- list(
  list(label = "Daily Temp",   win = 1, exp = "tmp", ref = "50%"),
  list(label = "7-Day Temp",   win = 2, exp = "tmp", ref = "50%"),
  list(label = "Daily Precip", win = 1, exp = "tp",  ref = tp_ref_imm),
  list(label = "7-Day Precip", win = 2, exp = "tp",  ref = tp_ref_cum)
)
tot_list <- list(tot1, tot2)   # indexed by spec$win
sub_list <- list(sub1, sub2)

get_rr_table <- function(meta, city = "ALL", reference = "50%") {
  res  <- city_con_res(meta = meta, city = city, reference = reference, xlab = "", title = "")$df
  pred <- meta_predict(meta, city = city, reference = reference, blup = FALSE)
  list2env(pred, envir = environment())
  
  idx <- match(pred$at, cp$predvar)
  overall_df <- data.frame(pred = perc,
                           RR   = cp$allRRfit[idx],
                           low  = cp$allRRlow[idx],
                           high = cp$allRRhigh[idx])
  
  add_fmt <- function(d) mutate(d, rr_fmt = fmt_rr(RR, low, high), sig = is_sig(low, high, RR))
  
  filtered         <- res |> filter(keep_pred(pred)) |> add_fmt()
  overall_filtered <- overall_df |> filter(pred %in% filtered$pred) |> add_fmt()
  
  make_wide <- function(val_col) {
    filtered |>
      dplyr::select(pred, lag, {{ val_col }}) |>
      pivot_wider(names_from = lag, values_from = {{ val_col }}, names_prefix = "Lag=") |>
      rename(Percentile = pred) |>
      left_join(
        overall_filtered |> dplyr::select(pred, {{ val_col }}) |>
          rename(Percentile = pred, Cumulative = {{ val_col }}),
        by = "Percentile")
  }
  list(tab = make_wide(rr_fmt), sig = make_wide(sig))
}

get_site_tables <- function(meta, locs, reference = "50%") {
  setNames(lapply(c("ALL", locs), function(l) get_rr_table(meta, city = l, reference = reference)),
           c("Pooled", locs))
}

get_sub_table <- function(df, filter_perc = TRUE) {
  level_col <- setdiff(names(df), c("perc", "exp", "fit", "lo", "hi"))
  
  long <- df |>
    rename(level = all_of(level_col)) |>
    mutate(Percentile = as.numeric(sub("%", "", perc)),
           rr_fmt     = fmt_rr(fit, lo, hi),
           sig        = is_sig(lo, hi, fit))
  if (filter_perc) long <- filter(long, keep_pred(Percentile))
  
  wide <- function(val_col) {
    long |>
      dplyr::select(Percentile, level, {{ val_col }}) |>
      pivot_wider(names_from = level, values_from = {{ val_col }})
  }
  list(tab = wide(rr_fmt), sig = wide(sig))
}

get_group_tables <- function(group, filter_perc = TRUE) {
  setNames(lapply(specs, function(s) get_sub_table(sub_list[[s$win]]$res_data[[group]][[s$exp]], filter_perc)),
           vapply(specs, `[[`, "", "label"))
}

write_stacked_excel <- function(sheets, filepath, col_width = 18) {
  wb <- createWorkbook()
  
  title_style  <- createStyle(textDecoration = "bold", fgFill = "grey30", fontColour = "white")
  header_style <- createStyle(textDecoration = "bold", border = "TopBottom", borderStyle = "medium")
  centre_style <- createStyle(halign = "center")
  sig_style    <- createStyle(halign = "center", textDecoration = "bold")
  
  for (sheet_name in names(sheets)) {
    addWorksheet(wb, sheet_name)
    tab_list    <- sheets[[sheet_name]]
    current_row <- 1
    
    for (nm in names(tab_list)) {
      df      <- as.data.frame(tab_list[[nm]]$tab)
      sig_mat <- as.data.frame(tab_list[[nm]]$sig)
      
      writeData(wb, sheet_name, x = nm, startRow = current_row, startCol = 1)
      addStyle(wb, sheet_name, title_style, rows = current_row, cols = 1:ncol(df))
      current_row <- current_row + 1
      
      writeData(wb, sheet_name, x = df, startRow = current_row, startCol = 1, headerStyle = header_style)
      addStyle(wb, sheet_name, centre_style,
               rows = (current_row + 1):(current_row + nrow(df)),
               cols = 2:ncol(df), gridExpand = TRUE)
      
      for (r in seq_len(nrow(sig_mat))) {          # bold significant cells
        for (c in which(unlist(sig_mat[r, ]) %in% TRUE)) {
          addStyle(wb, sheet_name, sig_style, rows = current_row + r, cols = c, stack = TRUE)
        }
      }
      current_row <- current_row + nrow(df) + 2    # blank row between tables
    }
    setColWidths(wb, sheet_name, cols = 1,     widths = 12)
    setColWidths(wb, sheet_name, cols = 2:13,  widths = col_width)
  }
  saveWorkbook(wb, filepath, overwrite = TRUE)
}

sheets_site <- setNames(
  lapply(specs, function(s) get_site_tables(tot_list[[s$win]][[paste0("meta_", s$exp)]],
                                            locations, reference = s$ref)),
  vapply(specs, `[[`, "", "label")
)
write_stacked_excel(sheets_site, paste0(path_report, "II. dlnm/ST_rr_site_", export, ".xlsx"))

sheets_sub <- list(Gender = get_group_tables("gender"),
                   Age    = get_group_tables("age"))
write_stacked_excel(sheets_sub, paste0(path_report, "II. dlnm/ST_rr_pop.xlsx"), col_width = 20)