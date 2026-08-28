### ------------------------------------------------------------------------ ### 
###  DNLM – Additional Tables & Figures (HISP)
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


#### TABLES ####
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

output <- readRDS('2_DLNM/Results/cumulative_subgroup_models.rds')

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

