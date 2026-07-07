### ------------------------------------------------------------------------ ### 
###  DNLM – Tables & Figures (HISP)
### ------------------------------------------------------------------------ ### 

rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

#### Dependencies ####

source('2_DLNM/2_Utilities.R')

library(lmtest); library(flextable); library(officer); library(purrr); library(mvmeta)
library(kableExtra); library(dplyr); library(moments); library(stringr)
library(tidyr); library(openxlsx); library(ggplot2); library(patchwork); library(cowplot)

## LaTeX table builder
get_kbl_tab <- function(dat, title, note, colnames, label,
                        spanning_header  = NULL,
                        bold_rows        = NULL,
                        italic_rows      = NULL,
                        hline_after_rows = NULL) {
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

## Escape Unicode characters to LaTeX math
escape_latex_unicode <- function(x) {
  x <- gsub("±", "$\\pm$",    x, fixed = TRUE)
  x <- gsub("°", "$^\\circ$", x, fixed = TRUE)
  x <- gsub("µ", "$\\mu$",    x, fixed = TRUE)
  x <- gsub("³", "$^{3}$",    x, fixed = TRUE)
  x
}

## Empty ggplot carrying only a bold title (patchwork spacer workaround)
get_title <- function(title, size = 16) {
  ggplot() + labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = size, face = "bold"))
}

## Formatting helpers for numeric table cells
fmt_p   <- function(x) ifelse(is.na(x), "\u2013", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_num <- function(x, d = 1) ifelse(is.na(x), "\u2013", formatC(x, format = "f", digits = d))
fmt_int <- function(x) ifelse(is.na(x), "\u2013", as.character(as.integer(x)))

#### MAIN CODE ####

# Load data
load('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData')
load(paste0(path_data, '2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata'))

df_wide <- df_wide |>
  filter(!location %in% c('Bilaspur'),
         date >= as.Date('2015-02-14'),
         date <= as.Date('2020-02-14')) |>
  rename(rh = rh_mean)

h_cols <- c('closed', 'total_patients', 'location', 'F', 'M', 'under_5',
            'age5_14', 'age15_24', 'age25_64', 'over_65', 'general', 
            'ophthalmology', 'ent', 'gynaecology')

locations <- df_wide |> distinct(location, elev) |> arrange(elev) |> pull(location) |> as.character()

#### TABLES ####
# MT1: Descriptive statistics ----

get_desc <- function(df, h_cols, w_cols) {
  
  single_loc <- length(unique(df$location)) == 1
  
  if (single_loc) {
    h_cols  <- h_cols[h_cols != "location"]
    check   <- c('F', 'under_5', 'general', 't2m')
    labels  <- c('Gender', 'Age', 'Ward', 'Weather Characteristics')
  } else {
    check   <- c(locations[1], 'F', 'under_5', 'general', 't2m')
    labels  <- c('Location', 'Gender', 'Age', 'Ward', 'Weather Characteristics')
  }
  
  summarise_col <- function(col) {
    x         <- df[[col]]
    is_binary <- all(x %in% c(0, 1, NA), na.rm = TRUE)
    
    if (is.numeric(x)) {
      data.frame(
        Variable          = sub('_b$', '', col),
        total_visits      = sum(x, na.rm = TRUE),
        perc_visits       = if (is_binary)
          round(sum(x, na.rm = TRUE) / nrow(df) * 100)
        else
          round(sum(x, na.rm = TRUE) / sum(df$total_patients, na.rm = TRUE) * 100, 1),
        mean_daily_visits = if (is_binary) NA else round(mean(x, na.rm = TRUE)),
        sd_daily_visits   = if (is_binary) NA else round(sd(x,   na.rm = TRUE)),
        min_daily_visit   = if (is_binary) NA else round(min(x,  na.rm = TRUE)),
        max_daily_visit   = if (is_binary) NA else round(max(x,  na.rm = TRUE)),
        skew_daily_visit  = if (is_binary) NA else round(skewness(x, na.rm = TRUE), 1)
      )
    } else {
      df |>
        group_by(location) |>
        summarise(
          total_visits      = sum(total_patients,            na.rm = TRUE),
          mean_daily_visits = round(mean(total_patients,     na.rm = TRUE)),
          sd_daily_visits   = round(sd(total_patients,       na.rm = TRUE)),
          min_daily_visit   = round(min(total_patients,      na.rm = TRUE)),
          max_daily_visit   = round(max(total_patients,      na.rm = TRUE)),
          skew_daily_visit  = round(skewness(total_patients, na.rm = TRUE), 1),
          .groups = 'drop'
        ) |>
        mutate(perc_visits = round(total_visits / sum(total_visits) * 100, 1)) |>
        dplyr::select(Variable = location, total_visits, perc_visits, everything())
    }
  }
  
  desc <- bind_rows(lapply(c(h_cols, w_cols), summarise_col))
  
  base_desc <- desc |>
    mutate(
      Frequency = paste0(format(round(total_visits), big.mark = ','), ' (', perc_visits, '%)'),
      Average   = ifelse(
        is.na(sd_daily_visits), NA,
        paste0(round(mean_daily_visits),
               ifelse(Variable == 'tp', ' + ', ' ± '),
               round(sd_daily_visits),
               ' (', min_daily_visit, ', ', max_daily_visit, ')')
      )
    ) |>
    dplyr::select(Variable, Frequency, Average)
  
  # Clear frequency for raw weather rows
  base_desc$Frequency[base_desc$Variable %in% gsub('(_IMD_b|_p90|_p10)$', '', w_cols)] <- NA
  
  # Sort ward rows by descending total visits
  ward_cols   <- h_cols[11:28]
  ward_vars   <- desc$Variable[desc$Variable %in% ward_cols]
  ward_idx    <- which(base_desc$Variable %in% ward_vars)
  if (length(ward_idx) > 0) {
    base_desc[ward_idx, ] <- base_desc[ward_idx, ][
      order(as.numeric(gsub(",", "", sub(" \\(.*", "", base_desc$Frequency[ward_idx]))),
            decreasing = TRUE, na.last = TRUE), ]
  }
  
  # Insert section headings
  base_desc <- bind_rows(
    tibble(Variable = 'Hospital Characteristics', Frequency = NA, Average = NA),
    base_desc
  )
  positions <- which(base_desc$Variable %in% check)
  new_rows  <- tibble(Variable = labels, Frequency = NA, Average = NA)
  for (i in rev(seq_along(positions))) {
    pos       <- positions[i]
    base_desc <- bind_rows(base_desc[1:(pos - 1), ], new_rows[i, ], base_desc[pos:nrow(base_desc), ])
  }
  
  base_desc
}

base <- get_desc(df_wide, h_cols, w_cols) |> replace_exact_matches(nm)

header_rows  <- which(grepl("Weather|Hospital", base$Variable))
section_rows <- c(setdiff(which(is.na(base$Frequency) & is.na(base$Average)), header_rows),
                  2, 3, header_rows[2]:nrow(base))

base_display <- base |>
  mutate(
    Frequency = trimws(Frequency),
    Frequency = gsub("%", "\\%", Frequency, fixed = TRUE),
    across(everything(), escape_latex_unicode),
    across(everything(), ~ifelse(is.na(.), "", .)),
    Variable  = case_when(
      row_number() %in% header_rows  ~ paste0("\\textbf{", Variable, "}"),
      row_number() %in% section_rows ~ paste0("\\hspace{1em}", Variable),
      TRUE                           ~ paste0("\\hspace{2em}\\textit{", Variable, "}")
    )
  )

tab <- get_kbl_tab(
  base_display,
  title    = "Descriptive statistics of hospital and weather characteristics",
  note     = "Binary event percentages represent the proportion of occurrence 
  over time and hospital registration percentages represent the proportion of 
  total registrations. Averages are presented as daily total registrations.
  Precipitation values are reported with positive SDs only because of its 
  inverse exponential distribution.",
  colnames = c("Variable", "Frequency n (\\%)", "Average mean $\\pm$ SD (min, max)"),
  label    = "MT_base_desc"
) |>
  column_spec(1, width = "5cm", latex_column_spec = "p{5cm}") |>
  column_spec(2, width = "4cm", latex_column_spec = "p{4cm}") |>
  column_spec(3, width = "5cm", latex_column_spec = "p{5cm}")

writeLines(as.character(tab), paste0(path_report, "II. dlnm/MT1_base_desc.tex"))

# ST: Meta-regression results ----

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

blank_row <- data.frame(
  `Meta-predictors` = "", `Stat` = "", `df` = "", `p` = "", `AIC` = "",
  `Stat ` = "", `df ` = "", `p ` = "", `I2 (%)` = "",
  check.names = FALSE, stringsAsFactors = FALSE
)

output  <- readRDS('2_DLNM/Results/meta_tables.RData')
d_tmp   <- prep_display(output$tbl_tmp)
d_tp    <- prep_display(output$tbl_tp)
combined <- rbind(blank_row, d_tmp, blank_row, d_tp)
combined[1, 1] <- 'Temperature'
combined[19, 1] <- 'Precipitation'

combined <- combined |>
  mutate(`Meta-predictors` = ifelse(
    `Meta-predictors` %in% c('Temperature', 'Precipitation', ''),
    `Meta-predictors`,
    paste0("\\hspace{1em}", `Meta-predictors`)
  ))

tab <- get_kbl_tab(
  combined,
  title           = "Meta-regression results for temperature and precipitation models",
  note            = "LR test was used to test the significance of meta-predictors with the intercept-only model as reference. Cochran Q test was used to ttest the significance of residual heterogeneity with the null hypothesis as no heterogeneity. LR = likelihood ratio; AIC = Akaike information criterion; $I^2$ = heterogeneity statistic.",
  colnames        = c("Meta-predictors", "Stat", "df", "\\textit{p}", "AIC",
                      "Stat", "df", "\\textit{p}", "(\\%)"),
  label           = "MT_meta_reg",
  spanning_header = c(" " = 1, "LR test" = 3, "Model fits" = 1, "Cochran Q test" = 3, "$I^2$" = 1),
  bold_rows       = c(1, 19),
  italic_rows     = c(2, 20),
  hline_after_rows = 18
)
writeLines(as.character(tab), paste0(path_report, "II. dlnm/ST_meta_reg.tex"))

# ST: Temporal summary statistics ----

summarise_by_group <- function(group_var, my_vars) {
  map_dfr(my_vars, function(var) {
    is_binary <- all(df_wide[[var]] %in% c(0, 1, NA), na.rm = TRUE)
    df_wide |>
      group_by(across(all_of(group_var))) |>
      summarise(
        Count = if (is_binary) sum(.data[[var]], na.rm = TRUE) else NA_integer_,
        Total = sum(!is.na(.data[[var]])),
        Mean  = if (!is_binary) round(mean(.data[[var]], na.rm = TRUE)) else NA_integer_,
        SD    = if (!is_binary) round(sd(.data[[var]],   na.rm = TRUE)) else NA_integer_,
        .groups = 'drop'
      ) |>
      mutate(
        variable = var,
        Summary  = if (is_binary)
          paste0(ifelse(Total > 0, round(100 * Count / Total), NA_integer_), '%')
        else
          paste0(Mean, ' ± ', SD),
        Time = .data[[group_var]]
      ) |>
      dplyr::select(Time, variable, Summary)
  })
}

my_vars <- c(h_cols[-3], w_cols)

df_time <- bind_rows(
  summarise_by_group('year',    my_vars) |> mutate(Group = 'Year',    Time = as.character(Time)),
  summarise_by_group('month',   my_vars) |> mutate(Group = 'Month',   Time = as.character(Time)),
  summarise_by_group('weekday', my_vars) |> mutate(Group = 'Weekday', Time = as.character(Time))
) |>
  pivot_wider(names_from = variable, values_from = Summary) |>
  mutate(
    Group = factor(Group, levels = c("Year", "Month", "Weekday")),
    Time = factor(Time, levels = c(as.character(2015:2020),
                                   levels(df_wide$month), levels(df_wide$weekday)))
  ) |>
  arrange(Group, Time)

names(df_time) <- replace_exact_matches(names(df_time), nm)


# Escape special characters in column names, then wrap in \textbf{}
names(df_time) <- names(df_time) |>
  (\(x) gsub("%", "\\%", x, fixed = TRUE))() |>
  (\(x) gsub("µ", "$\\mu$", x, fixed = TRUE))() |>
  (\(x) gsub("°", "$^\\circ$", x, fixed = TRUE))() |>
  (\(x) gsub("³", "$^3$", x, fixed = TRUE))() |>
  (\(x) paste0("\\textbf{", str_remove_all(x, "\\\\textbf\\{|\\}"), "}"))()

df_time <- df_time |>
  mutate(across(where(is.character), ~gsub("%", "\\%", ., fixed = TRUE)),
         across(where(is.character), ~gsub("±", "$\\pm$", ., fixed = TRUE)))

cols      <- df_time[, 1:2]
data_cols <- df_time[, 3:ncol(df_time)]
splits    <- list(seq(1, 11), seq(12, 18), seq(19, 22))
tables    <- lapply(splits, function(idx) cbind(cols, data_cols[, idx, drop = FALSE]))
nms       <- paste0("t", seq_along(tables))

for (i in seq_along(nms)) {
  writeLines(as.character(get_kbl_tab(
    dat      = tables[[i]],
    title    = NULL, note = if (i == length(nms)) 'Example' else '\\\\textit{(continued on next page)}',
    colnames = names(tables[[i]]), label = NULL
  )), paste0(path_report, "II. dlnm/ST_time_desc_", nms[i], ".tex"))
}

# ST: Site-specific descriptive statistics ----

base_loc <- reduce(locations, function(acc, l) {
  res           <- get_desc(df_wide |> filter(location == l), h_cols, w_cols[1:9])
  names(res)[-1] <- paste0(names(res)[-1], ' (', l, ')')
  if (is.null(acc)) res else cbind(acc, res[-1])
}, .init = NULL) |>
  replace_exact_matches(nm)

header_rows  <- which(grepl("Weather|Hospital", base_loc$Variable))
all_na_rows  <- which(apply(base_loc[-1], 1, function(x) all(is.na(x))))
section_rows <- c(setdiff(all_na_rows, header_rows), 2, 3, header_rows[2]:nrow(base_loc))

base_display <- base_loc |>
  mutate(
    across(where(is.character), ~trimws(.)),
    across(where(is.character), ~gsub("%", "\\%", ., fixed = TRUE)),
    across(everything(), escape_latex_unicode),
    across(everything(), ~ifelse(is.na(.), "", .)),
    Variable = case_when(
      row_number() %in% header_rows  ~ paste0("\\textbf{", Variable, "}"),
      row_number() %in% section_rows ~ paste0("\\hspace{1em}", Variable),
      TRUE                           ~ paste0("\\hspace{2em}\\textit{", Variable, "}")
    )
  )

col_nms  <- names(base_loc)
groups   <- gsub(".*\\((.*)\\)", "\\1", col_nms[-1])
prefixes <- gsub("\\s*\\(.*\\)", "", col_nms[-1])
cols      <- base_display[, 1, drop = FALSE]
data_cols <- base_display[, 2:ncol(base_display)]
splits    <- split(seq_len(ncol(data_cols)), ceiling(seq_len(ncol(data_cols)) / 6))
tables    <- lapply(splits, function(idx) cbind(cols, data_cols[, idx, drop = FALSE]))
nms       <- paste0("t", seq_along(tables))

for (i in seq_along(nms)) {
  idx <- splits[[i]]
  writeLines(as.character(get_kbl_tab(
    dat             = tables[[i]],
    title           = NULL,
    note            = if (i == length(nms)) 'Example' else '\\\\textit{(continued on next page)}',
    colnames        = c("Variable", prefixes[idx]),
    label           = NULL,
    spanning_header = c(" " = 1, table(groups[idx])[unique(groups[idx])])
  )), paste0(path_report, "II. dlnm/ST_loc_desc_", nms[i], ".tex"))
}

# ST: Percentile lookup ####

datalist <- lapply(locations, function(region) df_dlnm[df_dlnm$location == region, ])

cen_tmp <- sapply(datalist,function(x) {
  round(quantile(x$mx2t,c(0:100)/100,na.rm=T), 5)
})
cen_tmp <- as.data.frame(cen_tmp)
names(cen_tmp) <- locations
cen_tp <- sapply(datalist,function(x) {
  q <- quantile(x$tp,c(0:100)/100,na.rm=T)
  round(exp(q)-2,5)
})
cen_tp <- as.data.frame(cen_tp)
names(cen_tp) <- locations#

wb <- createWorkbook()
addWorksheet(wb, "Temperature")
addWorksheet(wb, "Precipitation")
writeData(wb, "Temperature", cen_tmp, rowNames = TRUE)
writeData(wb, "Precipitation", cen_tp, rowNames = TRUE)
saveWorkbook(wb, paste0(path_report, "II. dlnm/ST_percentiles.xlsx"), overwrite = TRUE)

# ST: RR tables (Excel) ----

make_stacked_excel <- function(tab_list_tmp, tab_list_tp, filepath) {
  wb <- createWorkbook()
  
  write_sheet <- function(tab_list, sheet_name) {
    addWorksheet(wb, sheet_name)
    current_row <- 1
    
    for (nm in names(tab_list)) {
      df      <- as.data.frame(tab_list[[nm]]$tab)
      sig_mat <- as.data.frame(tab_list[[nm]]$sig)
      
      writeData(wb, sheet_name, x = nm, startRow = current_row, startCol = 1)
      addStyle(wb, sheet_name,
               style = createStyle(textDecoration = "bold", fgFill = "grey30", fontColour = "white"),
               rows = current_row, cols = 1:ncol(df))
      current_row <- current_row + 1
      
      writeData(wb, sheet_name, x = df, startRow = current_row, startCol = 1,
                headerStyle = createStyle(textDecoration = "bold",
                                          border = "TopBottom", borderStyle = "medium"))
      addStyle(wb, sheet_name,
               style = createStyle(halign = "center"),
               rows  = (current_row + 1):(current_row + nrow(df)),
               cols  = 2:ncol(df), gridExpand = TRUE)
      
      for (r in seq_len(nrow(sig_mat))) {
        for (c in which(sig_mat[r, ] == TRUE)) {
          addStyle(wb, sheet_name,
                   style = createStyle(halign = "center", textDecoration = "bold"),
                   rows = current_row + r, cols = c, stack = TRUE)
        }
      }
      current_row <- current_row + nrow(df) + 1
    }
    setColWidths(wb, sheet_name, cols = 1,     widths = 10)
    setColWidths(wb, sheet_name, cols = 2:13,  widths = 18)
  }
  
  write_sheet(tab_list_tmp, "Temperature")
  write_sheet(tab_list_tp,  "Precipitation")
  saveWorkbook(wb, filepath, overwrite = TRUE)
}

get_rr_table <- function(meta, city = "ALL", reference = "50%") {
  res <- city_con_res(meta = meta, city = city, reference = reference, xlab = "", title = "")$df
  
  list2env(meta, envir = environment())
  
  if (!is.null(meta$mv)) {
    coef_use  <- coef(mv); vcov_use <- vcov(mv)
    rows_dup  <- apply(cen, 1, function(row) any(duplicated(row)))
    cen_clean <- as.data.frame(cen[!rows_dup, ])
  } else {
    city_name <- if (city == "ALL") rownames(ymat)[1] else city
    coef_use  <- as.numeric(ymat[city_name, ]); vcov_use <- Slist[[city_name]]
    cen_clean <- cen[!(duplicated(cen) | duplicated(cen, fromLast = TRUE)), , drop = FALSE]
  }
  
  perc    <- as.numeric(sub("%", "", names(rowMeans(cen_clean))))
  ref_idx <- if (reference == '20%') paste0(min(perc), '%') else reference
  cp_at   <- rowMeans(cen_clean)
  cp      <- crosspred(b_mean, coef = coef_use, vcov = vcov_use, model.link = 'log',
                       at = cp_at, cen = cp_at[ref_idx])
  
  overall_df <- data.frame(x = perc, z = cp$allRRfit, z_low = cp$allRRlow, z_high = cp$allRRhigh)
  
  keep_x <- function(x) x <= 5 | x >= 95 |
    (x >= 5  & x <= 20 & x %% 5  == 0) |
    (x >= 80 & x <= 95 & x %% 5  == 0) |
    (x >= 20 & x <= 80 & x %% 10 == 0)
  
  fmt_rr  <- function(z, lo, hi) sprintf("%.2f (%.2f, %.2f)", z, lo, hi)
  is_sig  <- function(lo, hi, z) !(lo <= 1 & hi >= 1) & round(z, 2) != 1.00
  
  filtered <- res |>
    filter(keep_x(x)) |>
    mutate(rr_fmt = fmt_rr(z, z_low, z_high), sig = is_sig(z_low, z_high, z))
  
  overall_filtered <- overall_df |>
    filter(x %in% filtered$x) |>
    mutate(rr_fmt = fmt_rr(z, z_low, z_high), sig = is_sig(z_low, z_high, z))
  
  make_wide <- function(df_long, val_col, key_prefix = "Lag=") {
    df_long |>
      dplyr::select(x, y, {{ val_col }}) |>
      pivot_wider(names_from = y, values_from = {{ val_col }}, names_prefix = key_prefix) |>
      rename(Percentile = x) |>
      left_join(overall_filtered |> dplyr::select(x, {{ val_col }}) |> rename(Percentile = x, Cumulative = {{ val_col }}),
                by = "Percentile")
  }
  
  list(tab = make_wide(filtered, rr_fmt), sig = make_wide(filtered, sig))
}

output <- readRDS('2_DLNM/Results/subgroup_models.rds')

tab_rr_site_tmp <- setNames(
  lapply(c("ALL", locations), function(l) get_rr_table(output$total_patients$meta_tmp, city = l)),
  c("Pooled", locations)
)
tab_rr_site_tp <- setNames(
  lapply(c("ALL", locations), function(l) get_rr_table(output$total_patients$meta_tp,  city = l, reference='20%')),
  c("Pooled", locations)
)
make_stacked_excel(tab_rr_site_tmp, tab_rr_site_tp, paste0(path_report, "II. dlnm/ST_rr_site.xlsx"))

outcomes <- names(output)
tab_rr_pop_tmp <- setNames(lapply(outcomes[-1], function(o) get_rr_table(output[[o]]$meta_tmp, city = 'ALL')),
                           replace_exact_matches(outcomes[-1], nm))
tab_rr_pop_tp  <- setNames(lapply(outcomes[-1], function(o) get_rr_table(output[[o]]$meta_tp,  city = 'ALL')),
                           replace_exact_matches(outcomes[-1], nm))
make_stacked_excel(tab_rr_pop_tmp, tab_rr_pop_tp, paste0(path_report, "II. dlnm/ST_rr_pop.xlsx"))

#### FIGURES ####
# MF2: Pooled exposure-response ----

results <- readRDS('2_DLNM/Results/pooled_plots.RData')
res <- results$res; cont <- results$cont; exp <- results$exp

legend_plot <- ggplot() +
  geom_line(aes(x = 0, y = 0, linetype = "City-specific (BLUP)", colour = "City-specific (BLUP)"), linewidth = 1.2) +
  geom_line(aes(x = 0, y = 0, linetype = "Pooled",               colour = "Pooled"),               linewidth = 1.2) +
  scale_linetype_manual(name = "", values = c("City-specific (BLUP)" = "dashed", "Pooled" = "solid")) +
  scale_colour_manual(  name = "", values = c("City-specific (BLUP)" = "grey70", "Pooled" = "grey30")) +
  guides(colour = guide_legend(override.aes = list(linewidth = c(1, 1.5)))) +
  theme_void() +
  theme(legend.position  = "top", legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"), legend.text = element_text(size = 10))

fig <- plot_grid(
  res$p_er_overall & theme(legend.position = "none", plot.title = element_blank()),
  legend_plot, ncol = 1, rel_heights = c(1, 0.05)
)
ggsave(paste0(path_report, 'II. dlnm/MF2_pool_exp.png'), fig, width = 8, height = 4, dpi = 300, bg = 'white')

# MF3: Exposure-lag-response figures, pooled ----

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
    theme(plot.title = element_text(size = 13, face = "bold", hjust = -0.05))
}
contour_legend  <- function(ref_plot, bar_w = 3.5, bar_h = 0.3, txt = 7, ttl = 8) {
  get_legend(ref_plot +
               guides(fill = guide_colorbar(barwidth  = unit(bar_w, "cm"), barheight = unit(bar_h, "cm"),
                                            frame.colour = "black", ticks.colour = "black")) +
               theme(legend.position = "bottom",
                     legend.text  = element_text(size = txt),
                     legend.title = element_text(size = ttl, face = 'bold')))
}

make_figure_overall <- function(contours, exposure) {
  con_tmp <- strip_title(strip_legend(contours$tmp$ALL))
  con_tp  <- strip_title(strip_legend(contours$tp$ALL))
  lag_tmp <- lapply(exposure$tmp$ALL, function(p) tag_as_title(strip_legend(p)))
  lag_tp  <- lapply(exposure$tp$ALL,  function(p) tag_as_title(strip_legend(p)))
  
  legends <- wrap_plots(
    list(plot_spacer(),
         contour_legend(contours$tmp$ALL),
         get_legend(exposure$tmp$ALL[[1]] + theme(legend.position = "bottom"))),
    nrow = 1, widths = c(0.15, 1, 2)
  )
  
  all_rows <- wrap_plots(
    list(row_label("(A)"), wrap_plots(c(list(con_tmp), lag_tmp), nrow = 1),
         row_label("(B)"), wrap_plots(c(list(con_tp),  lag_tp),  nrow = 1)),
    nrow = 2, widths = c(0.05, 1)
  )
  wrap_plots(list(all_rows, legends), ncol = 1, heights = c(1, 0.1))
}

(fig <- make_figure_overall(cont, exp))
ggsave(paste0(path_report, 'II. dlnm/MF3_pool_exp-lag.png'), fig, width = 9, height = 7, dpi = 300, bg = 'white')

# SF: Exposure-lag-response figures, site-specific ----
make_figure_sites <- function(contours, exposure) {
  
  make_exp_fig <- function(con_list, lag_list, main_title) {
    legends <- wrap_plots(
      list(plot_spacer(),
           contour_legend(con_list$ALL, bar_w = 7, bar_h = 0.5, txt = 12, ttl = 14),
           get_legend(lag_list[[2]][[1]] + theme(legend.position = "bottom")),
           plot_spacer()),
      nrow = 1, widths = c(1.5, 1, 1, 1.5)
    ) & theme(plot.margin = margin(t = 20, b = -10, unit = "pt"))
    
    rows <- lapply(locations, function(l) {
      wrap_plots(c(list(strip_title(strip_legend(con_list[[l]]))),
                   lapply(lag_list[[l]], function(p) tag_as_title(strip_legend(p)))),
                 nrow = 1)
    })
    
    # all_cells <- do.call(c, lapply(seq_along(locations), function(i) list(row_label(locations[i]), rows[[i]])))
    all_cells <- vector("list", length(locations) * 2)
    for (i in seq_along(locations)) {
      all_cells[[2*i - 1]] <- row_label(locations[i])
      all_cells[[2*i]]     <- rows[[i]]
    }
    altitude_panel <- function(cells, title){
      wrap_plots(list(get_title(title),
                      wrap_plots(cells, ncol = 2, widths = c(0.2, 1))),
                 ncol = 1, heights = c(0.02, 1))
    }

    
    all <- wrap_plots(
      list(altitude_panel(all_cells[1:10],  'Lower altitude'),
           altitude_panel(all_cells[11:20], 'Higher altitude')),
      ncol = 2
    )
    wrap_plots(list(all, legends), ncol = 1, heights = c(1, 0.08)) +
      plot_annotation(title = main_title,
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  }
  
  list(
    tmp = make_exp_fig(contours$tmp, exposure$tmp, "Temperature"),
    tp  = make_exp_fig(contours$tp,  exposure$tp,  "Precipitation")
  )
}

fig <- make_figure_sites(cont, exp)

final_plot <- wrap_plots(
  list(wrap_plots(list(get_title("(A)", 24), fig$tmp), nrow = 1,  widths=c(0.1,1)),
       wrap_plots(list(get_title("(B)", 24), fig$tp), nrow = 1, widths=c(0.1,1))),
  ncol = 1)

ggsave(paste0(path_report, 'II. dlnm/SF_site_exp-lag.png'), final_plot, width = 15, height = 25, dpi = 300, bg = 'white')

# SF: Subgroup exposure-response ----

outcomes <- h_cols[c(2,4:length(h_cols))]

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
make_subgroup_plot <- function(data, x_lab, pal, xlim = NULL, ylim = NULL) {
  p <- ggplot(data, aes(x = pred, y = RR, colour = Subgroup, fill = Subgroup, group = Subgroup)) +
    geom_hline(yintercept = 1, color = "gray40") +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1, colour = NA) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(name = NULL, values = pal) +
    scale_fill_manual(  name = NULL, values = pal) +
    labs(x = x_lab, y = "RR") +
    theme_classic() +
    theme(legend.position = "bottom")
  
  if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim, breaks = pretty(xlim, 5))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, xlim = range(data$pred, na.rm = TRUE)) +
      scale_y_continuous(breaks = pretty(ylim, 5))
  p
}
make_panel <- function(groups, pal) {
  
  temp_data   <- filter(all_df, exposure == "Temperature",   Subgroup %in% groups) |> droplevels()
  precip_data <- filter(all_df, exposure == "Precipitation", Subgroup %in% groups) |> droplevels()
  
  p_temp   <- make_subgroup_plot(temp_data,   "Temperature (%)",   pal, xlim = c(0, 20))
  p_precip <- make_subgroup_plot(precip_data, "Precipitation (%)", pal,
                                 ylim = range(c(precip_data$low, precip_data$high), na.rm = TRUE))
  
  wrap_plots(
    wrap_plots(list(p_temp   + theme(legend.position = "none"),
                    p_precip + theme(legend.position = "none")), nrow = 1),
    get_legend(p_temp),
    ncol = 1, heights = c(1, 0.1)
  )
}

all_df <- map_dfr(outcomes, function(i) {
  label <- as.factor(replace_exact_matches(i, nm))
  bind_rows(
    get_pooled_df(output[[i]]$meta_tmp, '50%') |> mutate(Subgroup = label, exposure = "Temperature"),
    get_pooled_df(output[[i]]$meta_tp,  '20%') |> mutate(Subgroup = label, exposure = "Precipitation")
  )
}) |>
  mutate(pred = as.numeric(pred)) |>
  filter((exposure == "Temperature"   & pred <= 20) |
           (exposure == "Precipitation" & pred >= 80))

gender_groups <- c("Women", "Men")
age_groups    <- c("Under 5", "Ages 5-14", "Ages 65+")
opd_groups    <- c('General', 'Ophthalmology', 'ENT', 'Gynaecology')
gender_pal <- c(scales::hue_pal()(2))
age_pal    <- c(scales::hue_pal()(3))
opd_pal    <- c(scales::hue_pal()(4))

panel_gender <- make_panel(gender_groups, gender_pal)
panel_age    <- make_panel(age_groups,    age_pal)
panel_opd    <- make_panel(opd_groups,    opd_pal)
final_plot <- wrap_plots(
  list(row_label("(A)"), wrap_plots(panel_gender, nrow = 1),
       row_label("(B)"), wrap_plots(panel_age, nrow = 1),
       row_label("(C)"), wrap_plots(panel_opd,  nrow = 1)),
  ncol = 1, heights = c(0.05, 1,0.05, 1,0.05, 1)
)

ggsave(paste0(path_report, 'II. dlnm/SF_pop_exp-lag.png'), final_plot, width = 6, height = 12, dpi = 300, bg = 'white')
