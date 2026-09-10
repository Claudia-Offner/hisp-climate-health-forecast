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

output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_meta_tables.RData'))
d_tmp   <- prep_display(output$tbl_tmp)
d_tp    <- prep_display(output$tbl_tp)
immediate <- rbind(blank_row, d_tmp, blank_row, d_tp)
immediate[1, 1] <- paste0(str_to_title(model1), ' Temperature')
immediate[19, 1] <- paste0(str_to_title(model1), ' Precipitation')
immediate <- immediate |>
  mutate(`Meta-predictors` = ifelse(
    `Meta-predictors` %in% c(paste0(str_to_title(model1), ' Temperature'),
                             paste0(str_to_title(model1), ' Precipitation'), ''),
    `Meta-predictors`,
    paste0("\\hspace{1em}", `Meta-predictors`)
  ))


output  <- readRDS(paste0('2_DLNM/Results/', export, '_', model1, '_meta_tables.RData'))
d_tmp   <- prep_display(output$tbl_tmp)
d_tp    <- prep_display(output$tbl_tp)
cumulative <- rbind(blank_row, d_tmp, blank_row, d_tp)
cumulative[1, 1] <- paste0(str_to_title(model2), ' Temperature')
cumulative[19, 1] <- paste0(str_to_title(model2), ' Precipitation')
cumulative <- cumulative |>
  mutate(`Meta-predictors` = ifelse(
    `Meta-predictors` %in% c(paste0(str_to_title(model2), ' Temperature'),
                             paste0(str_to_title(model2), ' Precipitation'), ''),
    `Meta-predictors`,
    paste0("\\hspace{1em}", `Meta-predictors`)
  ))

combined <- rbind(immediate, cumulative)

tab <- get_kbl_tab(
  combined,
  title           = "Meta-regression results for immediate and cumulative weather models",
  note            = "LR test was used to test the significance of meta-predictors with the intercept-only model as reference. Cochran Q test was used to test the significance of residual heterogeneity with the null hypothesis as no heterogeneity. LR = likelihood ratio; AIC = Akaike information criterion; $I^2$ = heterogeneity statistic.",
  colnames        = c("Meta-predictors", "Stat", "df", "\\textit{p}", "AIC",
                      "Stat", "df", "\\textit{p}", "(\\%)"),
  label           = "MT_meta_reg",
  spanning_header = c(" " = 1, "LR test" = 3, "Model fits" = 1, "Cochran Q test" = 3, "$I^2$" = 1),
  bold_rows       = c(1, 19, 37, 55),
  italic_rows     = c(2, 20, 38, 56),
  hline_after_rows = c(18, 36, 54)
)

tab <- sub(
  "\\\\multicolumn\\{9\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}",
  "\\\\multicolumn{9}{p{0.8\\\\textwidth}}{\\\\rule{0pt}{1em}",
  tab
)

writeLines(as.character(tab), paste0(path_report, "II. dlnm/ST_meta_reg_", export, ".tex"))

