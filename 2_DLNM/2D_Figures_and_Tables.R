### ------------------------------------------------------------------------ ### 
###  DNLM - Multiple city meta- and uni-analysis (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

####  DEPENDENCIES  ####
source('2_DLNM/2_Utilities.R')

# # Load packages
# library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
# library(ggplot2); library(dplyr); library(purrr); library(patchwork); library(openxlsx);
# library(colorspace); library(png); library(grid); library(reshape2); library(cowplot)
# library(TTR)

library(lmtest); library(flextable); library(officer); library(purrr)
library(kableExtra); library(dplyr); library(moments); library(stringr)

# Function to create basea latex table
get_kbl_tab <- function(dat, title, note, colnames, label,
                        spanning_header = NULL,
                        bold_rows = NULL,
                        italic_rows = NULL,
                        hline_after_rows = NULL){
  k <- kableExtra::kbl(
    dat,
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    linesep   = "",
    col.names = colnames,
    escape    = FALSE,
    caption   = title,
    label     = label
  ) |>
    kableExtra::kable_styling(
      latex_options = c("repeat_header", "hold_position"),
      font_size     = 9,
      full_width    = FALSE,
      position      = "center"
    )
  
  # Add spanning header row if provided
  if (!is.null(spanning_header)) {
    k <- k |> kableExtra::add_header_above(spanning_header, escape=FALSE)
  }
  
  # Bold specific rows
  if (!is.null(bold_rows)) {
    k <- k |> kableExtra::row_spec(bold_rows, bold = TRUE)
  }
  
  # Italic specific rows
  if (!is.null(italic_rows)) {
    k <- k |> kableExtra::row_spec(italic_rows, italic = TRUE)
  }
  
  # Extra hlines after specific rows
  if (!is.null(hline_after_rows)) {
    for (r in hline_after_rows) {
      k <- k |> kableExtra::row_spec(r, extra_latex_after = "\\midrule")
    }
  }
  
  k |> kableExtra::footnote(
    general           = note,
    general_title     = "",
    footnote_as_chunk = TRUE,
    escape            = FALSE
  )
}

escape_latex_unicode <- function(x) {
  x <- gsub("±", "$\\pm$",      x, fixed = TRUE)
  x <- gsub("°", "$^\\circ$",   x, fixed = TRUE)
  x <- gsub("µ", "$\\mu$",      x, fixed = TRUE)
  x <- gsub("³", "$^{3}$",      x, fixed = TRUE)
  x
}

# Function to get empty plot with title, for nested wrap_plot work around
get_title <- function(title, size=16){
  
  ggplot() +
    labs(title=title) +
    theme_void() + 
    theme(plot.title=element_text(hjust=0.5, size=size, face="bold"))
  
}

#### MAIN CODE ####

# Load data and parameters
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData'))
# load(paste0(path_project, '/ind_sens_params.RData'))
# df <- df_dlnm; rm(df_dlnm)
# 



#### TABLES ####
# MT1: Overall summary statistics (portrait) #####

# Function to calculate descriptive statistics based on select columns
get_desc <- function(df, h_cols, w_cols) {
  
  # Handling for location specific tables
  if (length(unique(df$location)) == 1){
    h_cols <- h_cols[h_cols != "location"]
    check <- c('F', 'under_5', 'general', 't2m')
    labels <- c('Gender', 'Age', 'Ward', 'Weather Characteristics')
  } else{
    check <- c(as.character(location_levels[1]), 'F', 'under_5', 'general', 't2m')
    labels <- c('Location', 'Gender', 'Age', 'Ward', 'Weather Characteristics')
  }
  
  my_vars <- c(h_cols, w_cols)
  
  # Get descriptive stats 
  (desc <- my_vars %>%
      lapply(function(col) {
        x <- df[[col]]
        
        
        if(is.numeric(x)){
          
          # Check if binary
          is_binary <- all(x %in% c(0, 1, NA), na.rm = TRUE)
          
          perc_val <- if (is_binary) {
            # For binary: percent of times '1' appears relative to number of rows
            round((sum(x, na.rm = TRUE) / nrow(df)) * 100)
          } else {
            # For non-binary: sum of x relative to total patients, rounded to 1 decimal
            round((sum(x, na.rm = TRUE) / sum(df$total_patients, na.rm = TRUE)) * 100, 1)
          }
          
          mean_val <- if (is_binary) {
            NA
          } else {
            round(mean(x, na.rm = TRUE))        # round mean to whole number
          }
          
          sd_val <- if (is_binary) {
            NA
          } else {
            round(sd(x, na.rm = TRUE))          # round SD to whole number
          }
          
          min_val <- if (is_binary) {
            NA  
          } else {
            round(min(x,na.rm = TRUE))
          }
          
          max_val <- if (is_binary) {
            NA  
          } else {
            round(max(x,na.rm = TRUE))
          }
          
          skew_val <- if (is_binary) {
            NA  
          } else {
            round(skewness(x,na.rm = TRUE), 1)
          }
          
          data.frame(
            Variable = sub('_b$', '', col),
            total_visits = sum(x, na.rm = TRUE),
            perc_visits = perc_val,
            mean_daily_visits = mean_val,
            sd_daily_visits = sd_val,
            min_daily_visit = min_val,
            max_daily_visit = max_val,
            skew_daily_visit = skew_val
          )
          
          
        } else {
          
          df %>%
            group_by(location) %>%
            summarise(
              total_visits = sum(total_patients, na.rm = TRUE),
              mean_daily_visits = round(mean(total_patients, na.rm = TRUE)),
              sd_daily_visits = round(sd(total_patients, na.rm = TRUE)),
              min_daily_visit = round(min(total_patients, na.rm = TRUE)),
              max_daily_visit = round(max(total_patients, na.rm = TRUE)),
              skew_daily_visit = round(skewness(total_patients, na.rm = TRUE))
            ) %>%
            mutate(
              perc_visits = round((total_visits / sum(total_visits)) * 100, 1)
            ) %>%
            dplyr::select(Variable=location, total_visits, perc_visits, everything())
          
        }
        
      }) %>%
      bind_rows())
  
  
  (base_desc <- desc %>%
      mutate(
        Frequency = paste0(format(round(total_visits), big.mark = ','), ' (', perc_visits, '%)'),
        Average = ifelse(
          is.na(sd_daily_visits),
          NA,
          paste0(
            round(mean_daily_visits),
            ifelse(Variable == 'tp', ' + ', ' ± '),
            round(sd_daily_visits),
            ' (', min_daily_visit, ', ', max_daily_visit, ')'
          )
          
        ),
        # Range = ifelse(
        #   is.na(skew_daily_visit),
        #   NA,
        #   # paste0(min_daily_visit, ', ',max_daily_visit, ' (', skew_daily_visit, ')')
        #   paste0(sd_daily_visits, ' (', skew_daily_visit, ')')
        # )
      ) %>%
      dplyr::select(Variable, Frequency, Average))
  
  
  # Remove frequency values for raw weather data
  base_desc$Frequency[base_desc$Variable %in%  gsub('(_IMD_b|_p90|_p10)$', '', w_cols)] <- NA
  
  # Sort ward rows by total_visits descending
  ward_cols <- h_cols[11:28]
  general_idx <- which(desc$Variable == "general")
  if (length(general_idx) > 0) {
    # Find all ward rows (from 'general' onward until next non-ward section)
    ward_vars <- desc$Variable[desc$Variable %in% ward_cols]
    ward_rows <- base_desc[base_desc$Variable %in% ward_vars, ]
    non_ward_rows <- base_desc[!base_desc$Variable %in% ward_vars, ]
    ward_idx  <- which(base_desc$Variable %in% ward_vars)
    
    base_desc[ward_idx, ] <- ward_rows[order(
      as.numeric(gsub(",", "", sub(" \\(.*", "", ward_rows$Frequency))), 
      decreasing = TRUE, na.last = TRUE), ]
  }
  
  # FORMATING
  # Insert headings at fixed positions
  base_desc <- rbind(
    tibble(Variable = 'Hospital Characteristics', Frequency = NA, Average = NA),
    base_desc
  )
  # Find positions to insert before
  positions <- which(base_desc$Variable %in% check)
  new_rows <- tibble(
    Variable = labels,
    Frequency = NA,
    Average = NA,
  )
  # Insert rows in reverse order so positions don't shift
  for (i in rev(seq_along(positions))) {
    pos <- positions[i]
    base_desc <- bind_rows(
      base_desc[1:(pos-1), ],
      new_rows[i, ],
      base_desc[pos:nrow(base_desc), ]
    )
  }
  
  return(base_desc)
}

# Load and clean data
load(paste0(path_data, '2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata'))
df_wide <- df_wide %>% filter(location != 'Bilaspur')
df_wide <- df_wide[df_wide$date >= as.Date('2015-02-14') & df_wide$date <= as.Date('2020-02-14'),]
names(df_wide) <- gsub("rh_mean", "rh", names(df_wide))

# Set cols of interest
h_cols <- c('closed', 'total_patients', 'location', 'F', 'M', 'under_5', 
            'age5_14', 'age15_24', 'age25_64', 'over_65', 'general', 'paediatrics', 
            'ophthalmology', 'gynaecology', 'orthopaedics', 'ent', 'dermatology', 
            'dentistry', 'emergency')
location_levels <- unique(df_wide$location)

# Calculate statistics
base <- get_desc(df_wide, h_cols, w_cols)
base <- replace_exact_matches(base, nm)

# Format indentation
header_rows     <- which(grepl("Weather|Hospital", base$Variable))
section_rows    <- c(setdiff( which(is.na(base$Frequency) & is.na(base$Average)), header_rows), 2,3, header_rows[2]:nrow(base))
base_display <- base %>%
  mutate(
    Frequency = trimws(Frequency),
    Frequency = gsub("%", "\\%", Frequency, fixed = TRUE),
    across(everything(), escape_latex_unicode),
    across(everything(), ~ifelse(is.na(.), "", .)),
    Variable = case_when(
      row_number() %in% header_rows    ~ paste0("\\textbf{", Variable, "}"),
      row_number() %in% section_rows ~ paste0("\\hspace{1em}", Variable),
      TRUE                              ~ paste0("\\hspace{2em}\\textit{", Variable, "}")
    )
  )

# Create latex table
tab <- get_kbl_tab(base_display, 
                  title    = "Descriptive statistics of hospital and weather characteristics",
                  note     = "Example",
                  colnames = c("Variable", "Frequency n (\\%)", "Average mean $\\pm$ SD (min, max)"),
                  label    = "MT_base_desc"
) |>
  column_spec(1, width = "5cm", latex_column_spec = "p{5cm}") |>
  column_spec(2, width = "4cm", latex_column_spec = "p{4cm}") |>
  column_spec(3, width = "5cm", latex_column_spec = "p{5cm}") 

# Save
writeLines(as.character(tab), paste0(path_report, "II. dlnm/MT1_base_desc.tex"))


# MT2: Meta-regression results (portrait) ####

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
fmt_p   <- function(x) ifelse(is.na(x), "\u2013", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_num <- function(x, d = 1) ifelse(is.na(x), "\u2013", formatC(x, format = "f", digits = d))
fmt_int <- function(x) ifelse(is.na(x), "\u2013", as.character(as.integer(x)))

# Load tables
output <- readRDS(paste0('2_DLNM/Results/meta_tables.RData'))

# Format flextable
d_tmp <- prep_display(output[[1]]) #$tbl_tmp)
d_tp  <- prep_display(output[[2]]) #$tbl_tp)
combined <- rbind(section_row(), d_tmp, section_row(), d_tp)
combined[1,1] <- 'Temperature'
combined[7,1] <- 'Precipitation'
tab <- get_kbl_tab(
  dat              = combined,
  title            = "Meta-regression results for temperature and precipitation models",
  note             = "LR = likelihood ratio; AIC = Akaike information criterion; $I^2$ = heterogeneity statistic.",
  colnames         = c("Meta-predictors", "Stat", "df", "\\textit{p}", "AIC",
                       "Stat", "df", "\\textit{p}", "(\\%)"),
  label            = "MT_meta_reg",
  spanning_header  = c(" " = 1, "LR test" = 3, "Model fits" = 1,
                       "Cochran Q test" = 3, "$I^2$" = 1),
  bold_rows        = c(1, 7),
  italic_rows      = c(2, 8),
  hline_after_rows = c(6)
)


# Write to file
writeLines(as.character(tab), paste0(path_report, "II. dlnm/MT2_meta_reg.tex"))


# ST: Temporal summary statistics (landscape) ####

# Function to summarize by a grouping variable
summarise_by_group <- function(group_var, my_vars) {
  map_dfr(my_vars, function(var) {
    x <- df_wide[[var]]
    is_binary <- all(x %in% c(0, 1, NA), na.rm = TRUE)
    
    df_wide %>%
      group_by(across(all_of(group_var))) %>%
      summarise(
        Count = if (is_binary) sum(.data[[var]], na.rm = TRUE) else NA_integer_,
        Total = sum(!is.na(.data[[var]])),
        Mean = if (!is_binary) round(mean(.data[[var]], na.rm = TRUE)) else NA_integer_,
        SD = if (!is_binary) round(sd(.data[[var]], na.rm = TRUE)) else NA_integer_,
        .groups = 'drop'
      ) %>%
      mutate(
        variable = var,
        Summary = if (is_binary) {
          pct <- ifelse(Total > 0, round(100 * Count / Total), NA_integer_)
          paste0(pct, '%')
        } else {
          paste0(Mean, ' ± ', SD)
        },
        Time = .data[[group_var]]
      ) %>%
      dplyr::select(Time, variable, Summary)
  })
}

my_vars <- c(h_cols[-3], w_cols)

# Get summaries for each time grouping
summary_weekday <- summarise_by_group('weekday', my_vars) %>% mutate(Time_Group = 'Weekday')
summary_month <- summarise_by_group('month', my_vars) %>% mutate(Time_Group = 'Month')
summary_year <- summarise_by_group('year', my_vars) %>% mutate(Time_Group = 'Year')

# Combine all summaries and pivot
df <- bind_rows(
  summary_year    %>% mutate(Time = as.character(Time)),
  summary_month   %>% mutate(Time = as.character(Time)),
  summary_weekday %>% mutate(Time = as.character(Time))
) %>%
  dplyr::select(Time_Group, Time, variable, Summary) %>%
  pivot_wider(names_from = variable, values_from = Summary) %>%
  mutate(
    Time_Group = factor(Time_Group, levels = c("Year", "Month", "Weekday")),
    Time = factor(Time, levels = c(as.character(2015:2020), levels(df_wide$month), levels(df_wide$weekday)))
  ) %>%
  arrange(Time_Group, Time)

names(df) <- replace_exact_matches(names(df), nm)
colnames(df)[1] <- 'Group'

# Column name special characters - single \ in output = \\ in R string
names(df) <- gsub("%",  "\\%",          names(df), fixed = TRUE)
names(df) <- gsub("µ",  "$\\mu$",       names(df), fixed = TRUE)
names(df) <- gsub("°",  "$^\\circ$",  names(df), fixed = TRUE)
names(df) <- gsub("³",  "$^3$",       names(df), fixed = TRUE)

# Wrap in \textbf{}
colnames(df) <- paste0("\\textbf{", str_remove_all(colnames(df), "\\\\textbf\\{|\\}"), "}")

# Cell values - % needs \% and ± needs $\pm$
df <- df %>%
  mutate(across(where(is.character), ~gsub("%",  "\\%",      ., fixed = TRUE))) %>%
  mutate(across(where(is.character), ~gsub("±",  "$\\pm$",   ., fixed = TRUE)))

# Split wide table
cols <- df[, 1:2]
data_cols <- df[, 3:ncol(df)]
splits <- list(seq(1,11), seq(12,20), seq(21,27))
# split(seq_len(ncol(data_cols)), ceiling(seq_len(ncol(data_cols)) / 10))
tables <- lapply(splits, function(idx) cbind(cols, data_cols[, idx, drop = FALSE]))
names(tables) <- paste0("t", seq_along(tables))

# Set footer
foot <- 'Example'

for (i in names(tables)) {
  
  latex_code <- get_kbl_tab(
    dat       = tables[[i]],
    title     = NULL,
    note      = if (i == names(tables)[length(tables)]) foot else '\\\\textit{(continued on next page)}',
    colnames  = names(tables[[i]]),
    label     = NULL
  )
  
  writeLines(latex_code, paste0(path_report,"II. dlnm/ST_time_desc_", i, ".tex"))
  
}

# ST: Site-specific summary descriptive (landscape) #####

# Get table
locations <- as.character(location_levels)
base_loc <- NULL

for (l in locations) {
  sub <- df_wide %>% dplyr::filter(location == l)
  res <- get_desc(sub, h_cols, w_cols[c(1:9)])
  names(res)[-1] <- paste0(names(res)[-1], ' (', l, ')')
  
  if (is.null(base_loc)) {
    base_loc <- res
  } else {
    base_loc <- cbind(base_loc, res[-1])  # drop first column if it's repeated
  }
}
base_loc <- replace_exact_matches(base_loc, nm)

# Format indentation
header_rows  <- which(grepl("Weather|Hospital", base_loc$Variable))
all_na_rows  <- which(apply(base_loc[-1], 1, function(x) all(is.na(x))))
section_rows <- c(setdiff(all_na_rows, header_rows), 2, 3, header_rows[2]:nrow(base_loc))

base_display <- base_loc %>%
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

# Extract spanning header from column names like "Frequency (Una)"
col_nms  <- names(base_loc)
groups   <- gsub(".*\\((.*)\\)", "\\1", col_nms[-1])
prefixes <- gsub("\\s*\\(.*\\)", "", col_nms[-1])

# Split wide table
cols      <- base_display[, 1, drop = FALSE]
data_cols <- base_display[, 2:ncol(base_display)]

splits <- split(seq_len(ncol(data_cols)), ceiling(seq_len(ncol(data_cols)) / 6))
tables <- lapply(splits, function(idx) cbind(cols, data_cols[, idx, drop = FALSE]))
names(splits) <- paste0("t", seq_along(splits))
names(tables) <- paste0("t", seq_along(tables))

foot <- 'Example'

for (i in names(tables)) {
  
  idx <- splits[[i]]
  
  # Subset groups and prefixes to current split
  grp_sub <- groups[idx]
  pre_sub <- prefixes[idx]
  
  # Build spanning header for this subset
  spanning_i <- c(" " = 1, table(grp_sub)[unique(grp_sub)])
  
  # Build colnames for this subset
  colnames_i <- c("Variable", pre_sub)
  
  latex_code <- get_kbl_tab(
    dat             = tables[[i]],
    title           = NULL,
    note            = if (i == names(tables)[length(tables)]) foot else '\\\\textit{(continued on next page)}',
    colnames        = colnames_i,
    label           = NULL,
    spanning_header = spanning_i
  )
  
  writeLines(latex_code, paste0(path_report, "II. dlnm/ST_loc_desc_", i, ".tex"))
}

# ST: Percentile value look up (landscape) ####

# ST: RR with CI for all exposure-lags, across locations (excel for ref) ####

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

# Load models
output <- readRDS(paste0('2_DLNM/Results/subgroup_models.rds'))

# Get site specific temperature RR tables
tab_rr_site_tmp <- lapply(c("ALL", locations), function(l) {
  get_rr_table(meta=output$total_patients$meta_tmp, city=l)
})
names(tab_rr_site_tmp) <- c("Pooled", locations)

# Get site specific precipitation RR tables
tab_rr_site_tp <- lapply(c("ALL", locations), function(l) {
  get_rr_table(output$total_patients$meta_tp, city=l)
})
names(tab_rr_site_tp) <- c("Pooled", locations)

# Export
make_stacked_excel(tab_rr_site_tmp, tab_rr_site_tp, filepath=paste0(path_report, "II. dlnm/ST_rr_site.xlsx"))

# ST: RR with CI for all exposure-lags, across population groups (excel for ref) ####

outcomes <- names(output)
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
make_stacked_excel(tab_rr_pop_tmp, tab_rr_pop_tp, filepath=paste0(path_report, "II. dlnm/ST_rr_pop.xlsx"))

#### FIGURES ####
# MF2: Pooled exposure-response ####

results <- readRDS(paste0('2_DLNM/Results/pooled_plots.RData'))
res <- results$res
con <- results$con
exp <- results$ exp

total_er_overall <- res$p_er_overall & theme(legend.position = "none", plot.title = element_blank())

legend_plot <- ggplot() +
  geom_line(aes(x=0, y=0, linetype="City-specific (BLUP)", colour="City-specific (BLUP)"), linewidth=1.2) +
  geom_line(aes(x=0, y=0, linetype="Pooled", colour="Pooled"), linewidth=1.2) +
  scale_linetype_manual(name="", values = c("City-specific (BLUP)" = "dashed", "Pooled" = "solid")) +
  scale_colour_manual(name="",   values = c("City-specific (BLUP)" = "grey70", "Pooled" = "grey30")) +
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

fig <- make_figure_overall(con, exp)
ggsave(filename=paste0(path_report, 'II. dlnm/MF3_pool_exp-lag.png'),  plot=fig, width=9, height=7,  dpi=300, bg='white')

# SF: Site specific exposure-lag-response ####

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

fig <- make_figure_sites(con, exp)
ggsave(paste0(path_report, 'II. dlnm/SF_site_exp-lag_tmp.png'), fig$tmp, width=15.5, height=11, dpi=300, bg='white')
ggsave(paste0(path_report, 'II. dlnm/SF_site_exp-lag_tp.png'),  fig$tp,  width=15.5, height=11, dpi=300, bg='white')


# SF: Subgroup exposure-lag-response ####

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

ggsave(paste0(path_report, 'II. dlnm/SF_pop_exp-lag.png'), final_plot, width=8, height=5, dpi=300, bg='white')


