### ------------------------------------------------------------------------ ###
###  DNLM – Tables & Figures
### ------------------------------------------------------------------------ ###

#### DEPENDENCIES ####

# --- Table functions ####

# LaTeX longtable builder (booktabs). `hline_after_rows` inserts a \midrule.
get_kbl_tab <- function(dat, colnames, label, note = NULL, hline_after_rows = NULL) {
  k <- kbl(dat, format = "latex", booktabs = TRUE, longtable = TRUE, linesep = "",
           col.names = colnames, escape = FALSE, label = label) |>
    kable_styling(latex_options = c("repeat_header", "hold_position"),
                  font_size = 9, full_width = FALSE, position = "center")
  for (r in hline_after_rows) k <- k |> row_spec(r, extra_latex_after = "\\midrule")
  k |> footnote(general = note, general_title = "",
                footnote_as_chunk = TRUE, escape = FALSE)
}

# Unicode -> LaTeX math
escape_latex_unicode <- function(x) {
  x <- gsub("±", "$\\pm$",    x, fixed = TRUE)
  x <- gsub("°", "$^\\circ$", x, fixed = TRUE)
  x <- gsub("µ", "$\\mu$",    x, fixed = TRUE)
  x <- gsub("³", "$^{3}$",    x, fixed = TRUE)
  x
}

# Numeric cell formatting (NA -> en dash)
dash    <- function(x, out) ifelse(is.na(x), "\u2013", out)
fmt_p   <- function(x)         dash(x, ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_num <- function(x, d = 1)  dash(x, formatC(x, format = "f", digits = d))
fmt_int <- function(x)         dash(x, as.character(as.integer(x)))

prep_display <- function(df) {
  data.frame(
    `Meta-predictors` = df[["Meta-predictors"]],
    `Stat`   = fmt_num(df$LR_stat), `df`  = fmt_int(df$LR_df), `p`  = fmt_p(df$LR_p),
    `AIC`    = fmt_num(df$AIC),
    `Stat `  = fmt_num(df$Q_stat),  `df ` = fmt_int(df$Q_df),  `p ` = fmt_p(df$Q_p),
    `I2 (%)` = fmt_num(df$I2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# --- Figure functions ####

path_out     <- function(model, suffix){paste0("outputs/", export, "_", model, "_", suffix)}
save_fig     <- function(file, fig, w, h, dpi = 300){
  ggsave(paste0(path_report, file), fig, width = w, height = h, dpi = dpi, bg = "white")
}
strip_legend <- function(p){p + theme(legend.position = "none")}
strip_title  <- function(p){p + theme(plot.title = element_blank())}

# Move the title into a centred plot tag
tag_as_title <- function(p) {
  ttl <- ggplot_build(p)$plot$labels$title
  strip_title(p) +
    theme(plot.tag.position = c(0.6, 0.95),
          plot.tag = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    labs(tag = ttl)
}

# Empty ggplot used as a row label
row_label <- function(label) {
  ggplot() + labs(title = label) + theme_void() +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = -0.05))
}

contour_legend <- function(ref_plot, bar_w = 5, bar_h = 0.3, txt = 12, ttl = 13) {
  get_legend(
    ref_plot +
      guides(fill = guide_colorbar(barwidth = unit(bar_w, "cm"), barheight = unit(bar_h, "cm"),
                                   frame.colour = "black", ticks.colour = "black")) +
      theme(legend.position = "bottom",
            legend.text  = element_text(size = txt),
            legend.title = element_text(size = ttl, face = "bold")))
}

make_panel <- function(plots, title = "") {
  wrap_plots(plots, nrow = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))) &
    theme(legend.position = "none")
}

# Grob table with spanning headers for the meta-analysis tests
get_meta_plot <- function(table) {
  ths <- ttheme_default(
    core    = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 10.5),
                   bg_params = list(fill = NA, col = NA)),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 10.5),
                   bg_params = list(fill = NA, col = NA)),
    rowhead = list(fg_params = list(fontsize = 10.5))
  )
  tab_grob <- tableGrob(table, rows = NULL, theme = ths) |>
    gtable_add_rows(heights = grid::unit(1.2, "cm"), pos = 0)
  
  spans <- data.frame(label = c("", "LR Test", "Fit", "Cochran Q test", "I2"),
                      width = c(1, 3, 1, 3, 1))
  stopifnot(sum(spans$width) == ncol(table))
  spans$end   <- cumsum(spans$width)
  spans$start <- spans$end - spans$width + 1
  spans <- spans[spans$label != "", ]
  
  for (i in seq_len(nrow(spans))) {
    s <- spans[i, ]
    tab_grob <- tab_grob |>
      gtable_add_grob(textGrob(s$label, gp = gpar(fontface = "bold", fontsize = 10.5)),
                      t = 1, l = s$start, r = s$end, name = paste0("span-", i)) |>
      gtable_add_grob(segmentsGrob(x0 = 0.05, x1 = 0.95, y0 = 0.1, y1 = 0.1, gp = gpar(lwd = 1.5)),
                      t = 1, l = s$start, r = s$end, name = paste0("span-line-", i))
  }
  tab_grob
}

# Contour + site-level panels shared by the supplementary figures
drop_pooled <- function(res) {
  for (v in names(res$cont)) res$cont[[v]]$ALL <- NULL
  res
}
meta_grob <- function(model, tbl)
  readRDS(path_out(model, "meta_tables.RData"))[[tbl]] |> prep_display() |> get_meta_plot()

site_row <- function(label, sites, meta)
  list(row_label(label),
       wrap_plots(list(sites, wrap_elements(full = meta)), nrow = 1, widths = c(2, 1)))

add_legends <- function(legends, all_rows)
  wrap_plots(list(legends, all_rows), ncol = 1, heights = c(0.025, 1))

#### MAIN CODE ####

export <- a 

load(paste0(path_data, "2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata"))
df_wide <- df_wide |>
  filter(location != "Bilaspur",
         date >= as.Date("2015-02-14"),
         date <= as.Date("2020-02-14")) |>
  rename(rh = rh_mean)

h_cols <- c("closed", "total_patients", "location", "F", "M", "under_5", "age5_64", "over_65")
w_cols <- c("t2m", "mn2t", "mx2t", "tp", "fg10", "rh", "pm25")
p_cols <- c("elev", "popu", "hosp_n", "insurance", "decision", "meducation", "hsize", "hwealth")

if (grepl("roll", export, fixed = TRUE)) {
  results1 <- readRDS(path_out(model1, "pooled_plots.RData"))
} else {
  results1 <- readRDS(path_out(model1, "pooled_plots.RData"))
  results2 <- readRDS(path_out(model2, "pooled_plots.RData"))
}

if (export == "t2m") {
  win_models <- c(model1, model2)   # immediate, cumulative
  tot_list   <- lapply(win_models, \(m) readRDS(path_out(m, "pooled_model.rds")))
  sub_list   <- lapply(win_models, \(m) readRDS(path_out(m, "subgroup_plots.RData")))
}

#### FIGURES ####

# MF2: Time series ####

if (export == "t2m") {
  
  get_time_series <- function(df_weekly, title = NULL, month_breaks = "4 months",
                              hosp_bottom = 400, band_bottom = 800,
                              weather_min = 0, weather_max = 30) {
    
    # Weather is drawn in a band sitting on top of the hospital bars
    band_top <- hosp_top <- round(band_bottom + band_bottom / 2.667)
    to_hosp    <- \(x) band_bottom + (x - weather_min) / (weather_max - weather_min) * (band_top - band_bottom)
    from_hosp  <- \(y) (y - band_bottom) / (band_top - band_bottom) * (weather_max - weather_min) + weather_min
    thin_breaks <- function(lo, hi, step) {
      b <- seq(lo, hi, by = step)
      if (length(b) > 6) seq(lo, hi, by = 2 * step) else b
    }
    
    date_breaks <- seq(min(df_weekly$date, na.rm = TRUE), max(df_weekly$date, na.rm = TRUE), by = month_breaks)
    
    ggplot(df_weekly, aes(x = date)) +
      geom_col(aes(y = hosp_week), fill = "grey70", alpha = 0.7, width = 7) +
      geom_smooth(aes(y = hosp_week, color = "Smoothed trend"), method = "loess", span = 0.07,
                  se = FALSE, linewidth = 0.7) +
      geom_line(aes(y = to_hosp(tp_week),  color = "Total Precipitation (mm)"), linewidth = 0.7) +
      geom_line(aes(y = to_hosp(t2m_week), color = "Mean Temperature (°C)"),    linewidth = 0.7) +
      scale_x_date(breaks = date_breaks, date_labels = "%b %Y") +
      scale_y_continuous(
        name = "Total patients",
        breaks = thin_breaks(hosp_bottom, band_bottom, 100),
        sec.axis = sec_axis(~ from_hosp(.), name = "Weather variables",
                            breaks = thin_breaks(weather_min, weather_max, 10))
      ) +
      coord_cartesian(ylim = c(hosp_bottom, hosp_top)) +
      scale_color_manual(
        values = c("Mean Temperature (°C)" = "#D74B3C", "Total Precipitation (mm)" = "#1D98EB",
                   "Smoothed trend" = "black"),
        breaks = c("Mean Temperature (°C)", "Total Precipitation (mm)", "Smoothed trend")) +
      labs(title = title, color = "") +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(hjust = 0.25, size = 12),
            axis.title.y.right = element_text(hjust = 0.1, size = 12),
            panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank())
  }
  
  df_weekly <- df |>
    mutate(date = floor_date(date, "week")) |>
    group_by(date) |>
    summarise(t2m_week  = mean(t2m, na.rm = TRUE),
              tp_week   = mean(tp, na.rm = TRUE),
              hosp_week = mean(total_patients, na.rm = TRUE))
  
  save_fig("MF2.png", get_time_series(df_weekly), 8, 4)
}

# MF3/MF4/S7 (t2m): exposure-lag-response ####

if (export == "t2m") {
  
  # Daily (1) and seven-day (2) panels for one exposure ("tmp" or "tp")
  make_exp_fig <- function(v, width = c(1, 3), center = FALSE) {
    cont1 <- results1$cont[[v]]; exp1 <- results1$exp[[v]]; lag1 <- results1$lag[[v]]
    cont2 <- results2$cont[[v]]; exp2 <- results2$exp[[v]]; lag2 <- results2$lag[[v]]
    
    # contour, cumulative, then whichever lag panels exist (P5 is absent for tp)
    panels <- \(cont, exp, lag) {
      p <- list(cont$ALL, exp, lag$P5$p_blup, lag$P95$p_blup)
      lapply(Filter(Negate(is.null), p), strip_legend)
    }
    p1 <- panels(cont1, exp1, lag1)
    p2 <- panels(cont2, exp2, lag2)
    plots <- c(p1, p2)
    
    if (center) {   # common y-axis for every non-contour panel
      idx    <- setdiff(seq_along(plots), c(1, length(p1) + 1))
      y_brks <- pretty(range(sapply(plots[idx], \(p) ggplot_build(p)$layout$panel_params[[1]]$y.range)), n = 5)
      plots[idx] <- lapply(plots[idx], \(p) p +
                             scale_y_continuous(breaks = y_brks, labels = scales::label_number(accuracy = 0.01)) +
                             coord_cartesian(ylim = range(y_brks), expand = FALSE))
    }
    
    plots <- Map(\(p, lab) p + labs(title = paste0(lab, p$labels$title)),
                 plots, paste0("(", LETTERS[seq_along(plots)], ") "))
    
    legends <- wrap_plots(list(contour_legend(cont2$ALL),
                               get_legend(exp1 + theme(legend.position = "bottom"))),
                          nrow = 1, widths = width)
    wrap_plots(list(wrap_plots(plots, nrow = 2), legends), ncol = 1, heights = c(1, 0.15))
  }  
  save_fig("MF3.png", make_exp_fig("tmp", width = c(1, 3)), 14, 7)
  save_fig("MF4.png", make_exp_fig("tp",  width = c(1, 2), center = TRUE), 10, 7)
  
  # S7: site-level contours + meta tests
  leg_ref  <- results1$cont$tp$ALL
  res1     <- drop_pooled(results1)
  res2     <- drop_pooled(results2)
  
  rows <- c(
    site_row("(A)", make_panel(res1$cont$tmp), meta_grob(model1, "tbl_tmp")),
    site_row("(B)", make_panel(res2$cont$tmp), meta_grob(model2, "tbl_tmp")),
    site_row("(C)", make_panel(res1$cont$tp),  meta_grob(model1, "tbl_tp")),
    site_row("(D)", make_panel(res2$cont$tp),  meta_grob(model2, "tbl_tp"))
  )
  all_rows <- wrap_plots(rows, nrow = 4, widths = c(0.05, 1))
  legends  <- wrap_plots(list(plot_spacer(), contour_legend(leg_ref), plot_spacer()),
                         nrow = 1, widths = c(0.05, 2, 1))
  
  save_fig("supplement/S7.png", add_legends(legends, all_rows), 20, 22, dpi = 100)
  
}

# MF5: Attributable fractions ####

if (export == "t2m") {
  
  plot_df <- rbind(readRDS(path_out(model1, "attr.RData")),
                   readRDS(path_out(model2, "attr.RData")))
  
  panel_colors  <- c(cold = "#1D98EB", heat = "#D74B3C", precip = "#8968CD")
  panel_color   <- function(p) {
    hit <- names(panel_colors)[vapply(names(panel_colors), grepl, NA, x = p, ignore.case = TRUE)]
    if (length(hit)) panel_colors[[hit[1]]] else "grey50"
  }
  fill_values   <- unlist(lapply(levels(plot_df$panel), \(p)
                                 setNames(c(panel_color(p), "grey50"), paste(p, c(FALSE, TRUE), sep = "."))))
  legend_colors <- c(Cold = panel_colors[["cold"]], Heat = panel_colors[["heat"]],
                     Precipitation = panel_colors[["precip"]], Total = "grey50")
  panel_labels  <- setNames(paste0("(", LETTERS[seq_along(levels(plot_df$panel))], ") ", levels(plot_df$panel)),
                            levels(plot_df$panel))
  
  fig <- ggplot(plot_df, aes(x = location, y = af * 100,
                             fill = interaction(panel, location == "TOTAL"))) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_errorbar(aes(ymin = af_low * 100, ymax = af_high * 100), width = 0.25) +
    facet_wrap(~panel, nrow = 2, labeller = labeller(panel = panel_labels)) +
    scale_fill_manual(values = fill_values) +
    new_scale_fill() +
    # invisible points, used only to draw the legend
    geom_point(data = data.frame(legend_group = names(legend_colors)),
               aes(x = 1, y = 0, fill = factor(legend_group, levels = names(legend_colors))),
               shape = 22, size = 0, inherit.aes = FALSE) +
    scale_fill_manual(name = NULL, values = legend_colors,
                      guide = guide_legend(override.aes = list(size = 5))) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    labs(x = NULL, y = "Attributable fraction (%)") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"),
          strip.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
          legend.position = "bottom",
          legend.box.spacing = unit(0, "pt"))
  
  save_fig("MF5.png", fig, 8, 5)
}

# S6: Subgroup exposure-response ####

if (export == "t2m") {
  
  res_imm <- sub_list[[1]]$res_plots
  res_cum <- sub_list[[2]]$res_plots
  plt_imm <- lapply(res_imm, strip_legend)
  plt_cum <- lapply(res_cum, strip_legend)
  
  # One row: daily/seven-day x gender/age for a single exposure
  sub_row <- function(gender, age, xlab, name, tags) {
    plots  <- list(plt_imm[[gender]], plt_cum[[gender]], plt_imm[[age]], plt_cum[[age]])
    titles <- paste0("(", tags, ") ", rep(c("Daily", "Seven-day"), 2), " ", name, ": ",
                     rep(c("Gender", "Age"), each = 2))
    wrap_plots(Map(\(p, t) p + labs(x = xlab, title = t), plots, titles), nrow = 1)
  }
  
  fig <- wrap_plots(
    wrap_plots(list(sub_row(1, 3, "Temperature Percentile (%)",   "temp",   LETTERS[1:4]),
                    sub_row(2, 4, "Precipitation Percentile (%)", "precip", LETTERS[5:8])),
               nrow = 2),
    wrap_plots(list(get_legend(res_imm$immediate_tmp_gen), get_legend(res_imm$immediate_tmp_age)), nrow = 1),
    ncol = 1, heights = c(1, 0.1)
  ) & theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
  
  save_fig("supplement/S6.png", fig, 10, 5)
}

# S8/S9: Curve comparisons ####

# Extract the pooled RR curve from a fitted meta model
get_curve_df <- function(model, x, ref, test) {
  cp  <- meta_predict(model[[x]], city = "ALL", reference = ref, blup = TRUE)
  idx <- match(cp$at, cp$cp$predvar)
  data.frame(test = test, perc = cp$perc, tmp = cp$at,
             rr = cp$cp$allRRfit[idx], lo = cp$cp$allRRlow[idx], hi = cp$cp$allRRhigh[idx])
}

# Overlay curves (`curves$test` must be a factor)
plot_curves <- function(curves, xlab, ttl = NULL) {
  cols <- setNames(scales::hue_pal()(nlevels(curves$test)), levels(curves$test))
  ylim <- c(min(curves$lo), max(curves$hi))
  
  ggplot(curves, aes(x = perc)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = test), alpha = 0.2, show.legend = FALSE) +
    geom_line(aes(y = rr, color = test), linewidth = 1.2) +
    geom_hline(yintercept = 1, color = "gray20") +
    scale_color_manual(name = NULL, values = cols) +
    scale_fill_manual(name = NULL, values = cols) +
    scale_x_continuous(breaks = pretty(curves$perc, n = 5)) +
    scale_y_continuous(breaks = pretty(ylim, n = 5)) +
    labs(x = xlab, y = "RR", title = ttl) +
    coord_cartesian(ylim = ylim, xlim = range(curves$perc)) +
    guides(color = guide_legend(nrow = 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
          axis.line = element_line(color = "gray20"),
          axis.text = element_text(color = "gray20"),
          axis.ticks = element_line(color = "gray20"),
          legend.position = "bottom")
}

# S8: Multi-day windows
if (grepl("roll6", export)) {
  
  roll_files  <- list.files("outputs/", pattern = "roll.*\\.rds$", full.names = TRUE)
  roll_files  <- roll_files[grepl("roll[2-6](?!\\d)", basename(roll_files), perl = TRUE)]
  roll_models <- setNames(lapply(roll_files, readRDS),
                          sub("_cumulative_pooled_model\\.rds$", "", basename(roll_files)))
  
  roll_curves <- function(x, xlab, ttl = NULL, ref = "50%") {
    curves <- bind_rows(lapply(names(roll_models), \(i) get_curve_df(roll_models[[i]], x, ref, i)))
    curves$test <- gsub("^roll(\\d+)$", "\\1 days", curves$test)
    curves$test <- factor(curves$test, unique(curves$test))
    plot_curves(curves, xlab, ttl)
  }
  
  fig <- wrap_plots(
    list(roll_curves("meta_tmp", "Temperature percentiles (%)",   ttl = "(A) Consecutive temperature"),
         roll_curves("meta_tp",  "Precipitation percentiles (%)", ttl = "(B) Consecutive precipitation", ref = tp_ref)),
    nrow = 1
  ) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  save_fig("supplement/S8.png", fig, 8, 4)
}

# S9: Different temperature measures
if (grepl("mn2t", export)) {
  
  roll_files  <- list.files("outputs/", pattern = "^(mn2t|mx2t|t2m).*_pooled_model\\.rds$", full.names = TRUE)
  roll_models <- setNames(lapply(roll_files, readRDS),
                          sub("_pooled_model\\.rds$", "", basename(roll_files)))
  var_labs    <- c(mn2t = "Minimum temperature", mx2t = "Maximum temperature", t2m = "Mean temperature")
  
  measure_curves <- function(type, ttl) {
    nms    <- grep(paste0("_", type, "$"), names(roll_models), value = TRUE)
    curves <- bind_rows(lapply(nms, \(i)
                               get_curve_df(roll_models[[i]], "meta_tmp", "50%", var_labs[[sub(paste0("_", type, "$"), "", i)]])))
    curves$test <- factor(curves$test, unname(var_labs))
    plot_curves(curves, "Temperature percentiles (%)", ttl)
  }
  
  fig <- wrap_plots(list(measure_curves("immediate",  "(A) Daily effects"),
                         measure_curves("cumulative", "(B) Seven-day effects")), nrow = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  save_fig("supplement/S9.png", fig, 8, 4)
}

#### TABLES ####

# S4A/S4B: Descriptive statistics ####

if (export == "t2m") {
  
  get_desc <- function(df, h_cols, w_cols, p_cols = NULL) {
    
    h_cols   <- h_cols[h_cols != "location"]
    all_cols <- c(h_cols, w_cols, p_cols)
    locs     <- as.character(unique(df$location))
    
    # Section headings sit above the first variable of each block
    first_vars <- c("F", "under_5", "t2m", "elev")
    headings   <- c("Hospital Characteristics", "Gender", "Age",
                    "Weather Characteristics", "Socio-geographic Characteristics")
    
    summarise_col <- function(col, sub_df) {
      x   <- sub_df[[col]]
      avg <- NA_character_
      if (!all(x %in% c(0, 1, NA), na.rm = TRUE) && !is.na(sd(x, na.rm = TRUE))) {
        m   <- mean(x, na.rm = TRUE)
        avg <- if (!col %in% p_cols) {
          paste0(round(m), " (", round(min(x, na.rm = TRUE)), ", ", round(max(x, na.rm = TRUE)), ")")
        } else if (m >= 0 && m <= 1) {
          as.character(round(m, 2) * 100)
        } else {
          as.character(round(m))
        }
      }
      data.frame(Variable = sub("_b$", "", col), Average = avg)
    }
    summarise_df <- \(sub_df) bind_rows(lapply(all_cols, summarise_col, sub_df = sub_df))
    
    overall <- summarise_df(df)
    by_loc  <- lapply(setNames(locs, locs), \(l) summarise_df(df[df$location == l, ])$Average)
    out     <- tibble(Variable = overall$Variable, Total = overall$Average, !!!by_loc)
    
    # Insert heading rows
    starts  <- c(1, match(first_vars, out$Variable))
    section <- findInterval(seq_len(nrow(out)), starts)
    blank   <- \(v) as_tibble(setNames(c(list(v), rep(list(NA_character_), ncol(out) - 1)), names(out)))
    bind_rows(lapply(seq_along(headings), \(i) bind_rows(blank(headings[i]), out[section == i, ])))
  }
  
  base        <- get_desc(df_wide, h_cols[-1], w_cols, p_cols) |> replace_exact_matches(nm)
  value_cols  <- setdiff(names(base), "Variable")
  header_rows <- which(grepl("Characteristics", base$Variable))
  section_rows <- c(setdiff(which(rowSums(!is.na(base[value_cols])) == 0), header_rows),
                    2, seq(header_rows[2], nrow(base)))
  
  base_display <- base |>
    mutate(
      across(all_of(value_cols), trimws),
      across(everything(), \(x) ifelse(is.na(x), "", escape_latex_unicode(x))),
      Variable = case_when(
        row_number() %in% header_rows  ~ paste0("\\textbf{", Variable, "}"),
        row_number() %in% section_rows ~ paste0("\\hspace{1em}", Variable),
        TRUE                           ~ paste0("\\hspace{2em}\\textit{", Variable, "}")
      )
    )
  
  build_tab <- function(disp, cols, label) {
    t <- get_kbl_tab(disp[, c("Variable", cols)], colnames = c("Variable", cols), label = label) |>
      column_spec(1, width = "4.5cm", latex_column_spec = "p{4.5cm}")
    for (i in seq_along(cols) + 1) t <- t |> column_spec(i, width = "2.25cm", latex_column_spec = "p{2.25cm}")
    
    t <- gsub("%", "\\%", t, fixed = TRUE)
    # Let footnote rows wrap within the table width instead of overflowing
    gsub("\\\\multicolumn\\{([0-9]+)\\}\\{l\\}",
         "\\\\multicolumn{\\1}{p{\\\\dimexpr\\\\linewidth-2\\\\tabcolsep\\\\relax}}", t)
  }
  
  # Split locations over two tables
  split_pt   <- ceiling(length(value_cols) / 2)
  col_groups <- list(value_cols[seq_len(split_pt)], value_cols[-seq_len(split_pt)])
  for (i in 1:2) {
    writeLines(as.character(build_tab(base_display, col_groups[[i]], paste0("MT_base_desc_t", i))),
               paste0(path_report, "supplement/S4", c("A", "B")[i], ".tex"))
  }
}

# S5: Percentile lookup ####

if (export == "t2m") {
  
  datalist    <- lapply(locations, \(region) df[df$location == region, ])
  lower_pct   <- c(0.01, 0.5, 1, 2, 3, 4, 5, 10, 15, 20, 25)
  upper_pct   <- c(75, 80, 85, 90, 95, 96, 97, 98, 99, 99.5, 99.99)
  percentiles <- c(lower_pct, upper_pct)
  
  # Max 2 dp, no trailing zeros (e.g. "0.01", "25", "99.5")
  pct_label <- \(p) sub("\\.?0+$", "", formatC(p, format = "f", digits = 2))
  
  # Mean across locations of the given quantiles (back-transformed)
  mean_centile <- function(var, transform) {
    tab <- sapply(datalist, \(x) round(transform(quantile(x[[var]], percentiles / 100, na.rm = TRUE)), 2))
    setNames(rowMeans(tab), paste0(pct_label(percentiles), "%"))
  }
  
  var_specs <- tibble::tribble(
    ~variable, ~type,       ~col,           ~transform,
    "t2m",     "Daily",     "t2m",          identity,
    "t2m",     "Seven-day", "t2m_roll7",    identity,
    "tp",      "Daily",     "tp_log",       \(q) exp(q) - 2,
    "tp",      "Seven-day", "tp_roll7_log", \(q) exp(q) - 10
  )
  
  cen <- as.data.frame(pmap(var_specs, \(variable, type, col, transform) mean_centile(col, transform)),
                       col.names = paste0("V", seq_len(nrow(var_specs))))
  names(cen) <- paste(var_specs$type, replace_exact_matches(var_specs$variable, nm))
  
  tab_df <- data.frame(Percentile = rownames(cen), lapply(cen, sprintf, fmt = "%.2f"),
                       check.names = FALSE)
  
  tab <- get_kbl_tab(tab_df, colnames = names(tab_df), label = "ST_perc",
                     note = "Please note that daily and seven-day precipitation have different log transformations, which skew the zero-values.",
                     hline_after_rows = length(lower_pct)) |>
    gsub("%", "\\%", x = _, fixed = TRUE) |>
    gsub("\\midrule\\\\", "\\midrule", x = _, fixed = TRUE)
  
  writeLines(as.character(tab), paste0(path_report, "supplement/S5.tex"))
}

# RR tables (Excel) ####

if (export == "t2m") {
  
  # Percentiles to report: every point in the tails, coarser steps in the middle
  keep_pred <- \(p) p <= 5 | p >= 95 |
    (p <= 20 & p %% 5 == 0) | (p >= 80 & p %% 5 == 0) | (p > 20 & p < 80 & p %% 10 == 0)
  fmt_rr <- \(rr, lo, hi) sprintf("%.2f (%.2f, %.2f)", rr, lo, hi)
  is_sig <- \(lo, hi, rr) !(lo <= 1 & hi >= 1) & round(rr, 2) != 1.00
  
  # One entry per sheet: window (1 = daily, 2 = seven-day), exposure, reference percentile
  specs <- list(
    list(label = "Daily Temp",   win = 1, exp = "tmp", ref = "50%"),
    list(label = "7-Day Temp",   win = 2, exp = "tmp", ref = "50%"),
    list(label = "Daily Precip", win = 1, exp = "tp",  ref = "55%"),
    list(label = "7-Day Precip", win = 2, exp = "tp",  ref = "20%")
  )
  spec_names <- vapply(specs, `[[`, "", "label")
  
  add_fmt <- \(d, rr, lo, hi) mutate(d, rr_fmt = fmt_rr({{ rr }}, {{ lo }}, {{ hi }}),
                                     sig    = is_sig({{ lo }}, {{ hi }}, {{ rr }}))
  
  get_rr_table <- function(meta, city = "ALL", reference = "50%") {
    res <- city_con_res(meta = meta, city = city, reference = reference, xlab = "", title = "")$df
    mp  <- meta_predict(meta, city = city, reference = reference, blup = FALSE)
    idx <- match(mp$at, mp$cp$predvar)
    
    filtered <- res |> filter(keep_pred(pred)) |> add_fmt(RR, low, high)
    overall  <- data.frame(pred = mp$perc, RR = mp$cp$allRRfit[idx],
                           low = mp$cp$allRRlow[idx], high = mp$cp$allRRhigh[idx]) |>
      filter(pred %in% filtered$pred) |> add_fmt(RR, low, high)
    
    make_wide <- \(val) filtered |>
      dplyr::select(Percentile = pred, lag, all_of(val)) |>
      pivot_wider(names_from = lag, values_from = all_of(val), names_prefix = "Lag=") |>
      left_join(dplyr::select(overall, Percentile = pred, Cumulative = all_of(val)), by = "Percentile")
    
    list(tab = make_wide("rr_fmt"), sig = make_wide("sig"))
  }
  
  get_sub_table <- function(df) {
    level_col <- setdiff(names(df), c("perc", "exp", "fit", "lo", "hi"))
    long <- df |>
      rename(level = all_of(level_col)) |>
      mutate(Percentile = as.numeric(sub("%", "", perc))) |>
      filter(keep_pred(Percentile)) |>
      add_fmt(fit, lo, hi)
    
    make_wide <- \(val) long |>
      dplyr::select(Percentile, level, all_of(val)) |>
      pivot_wider(names_from = level, values_from = all_of(val))
    
    list(tab = make_wide("rr_fmt"), sig = make_wide("sig"))
  }
  
  site_tables <- \(s) setNames(
    lapply(c("ALL", locations), \(l)
           get_rr_table(tot_list[[s$win]][[paste0("meta_", s$exp)]], city = l, reference = s$ref)),
    c("Pooled", locations))
  
  group_tables <- \(group) setNames(
    lapply(specs, \(s) get_sub_table(sub_list[[s$win]]$res_data[[group]][[s$exp]])),
    spec_names)
  
  # Several tables stacked per sheet; significant cells in bold
  write_stacked_excel <- function(sheets, filepath, col_width = 18) {
    wb <- createWorkbook()
    title_style  <- createStyle(textDecoration = "bold", fgFill = "grey30", fontColour = "white")
    header_style <- createStyle(textDecoration = "bold", border = "TopBottom", borderStyle = "medium")
    centre_style <- createStyle(halign = "center")
    sig_style    <- createStyle(halign = "center", textDecoration = "bold")
    
    for (sheet in names(sheets)) {
      addWorksheet(wb, sheet)
      row0 <- 1
      for (tbl_name in names(sheets[[sheet]])) {
        df  <- as.data.frame(sheets[[sheet]][[tbl_name]]$tab)
        sig <- as.matrix(sheets[[sheet]][[tbl_name]]$sig[-1])   # drop Percentile column
        
        writeData(wb, sheet, x = tbl_name, startRow = row0, startCol = 1)
        addStyle(wb, sheet, title_style, rows = row0, cols = 1:ncol(df))
        row0 <- row0 + 1
        
        writeData(wb, sheet, x = df, startRow = row0, startCol = 1, headerStyle = header_style)
        addStyle(wb, sheet, centre_style, rows = row0 + seq_len(nrow(df)), cols = 2:ncol(df), gridExpand = TRUE)
        
        hits <- which(sig, arr.ind = TRUE)
        if (nrow(hits) > 0) {
          addStyle(wb, sheet, sig_style, rows = row0 + hits[, 1], cols = hits[, 2] + 1,
                   gridExpand = FALSE, stack = TRUE)
        }
        row0 <- row0 + nrow(df) + 2
      }
      setColWidths(wb, sheet, cols = 1,    widths = 12)
      setColWidths(wb, sheet, cols = 2:13, widths = col_width)
    }
    saveWorkbook(wb, filepath, overwrite = TRUE)
  }
  
  write_stacked_excel(setNames(lapply(specs, site_tables), spec_names),
                      paste0(path_report, "supplement/ST_rr_site_", export, ".xlsx"))
  write_stacked_excel(list(Gender = group_tables("gender"), Age = group_tables("age")),
                      paste0(path_report, "supplement/ST_rr_pop.xlsx"), col_width = 20)
}