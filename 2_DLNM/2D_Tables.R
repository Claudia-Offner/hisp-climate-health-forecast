### ------------------------------------------------------------------------ ### 
###  DNLM - Multiple city meta- and uni-analysis (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list=ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/offne/Documents/GitHub/hisp-climate-health-forecast/')


####  DEPENDENCIES  ####
source('2_DLNM/DLNM - Dependencies.R')

# # Load packages
# library(dlnm); library(splines); library(xtable); library(mvmeta); library(MASS);
# library(ggplot2); library(dplyr); library(purrr); library(patchwork); library(openxlsx);
# library(colorspace); library(png); library(grid); library(reshape2); library(cowplot)
# library(TTR)

library(lmtest); library(flextable); library(officer)
library(kableExtra); library(dplyr)

get_kbl_tab <- function(dat, title, note, colnames, label,
                        spanning_header = NULL,  # named list: list(text, cols)
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
    k <- k |> kableExtra::add_header_above(spanning_header)
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
    general_title     = "Note:",
    footnote_as_chunk = TRUE,
    escape            = FALSE
  )
}

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
  
  vars <- c(h_cols, w_cols)
  
  # Get descriptive stats 
  (desc <- vars %>%
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

escape_latex_unicode <- function(x) {
  x <- gsub("±", "$\\pm$",      x, fixed = TRUE)
  x <- gsub("°", "$^\\circ$",   x, fixed = TRUE)
  x <- gsub("µ", "$\\mu$",      x, fixed = TRUE)
  x <- gsub("³", "$^{3}$",      x, fixed = TRUE)
  x
}



#### MAIN CODE ####

# Load data and parameters
load(paste0(path_data, '2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata'))
load(paste0('2_DLNM/ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_dlnm.RData'))
load(paste0(path_project, '/ind_sens_params.RData'))
df <- df_dlnm; rm(df_dlnm)

# CLEAN WIDE DATA

# Remove Bilaspur (time series too short)
df_wide <- df_wide %>% filter(location != 'Bilaspur')

# Remove COVID years for now
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-14')  
df_wide <- df_wide[df_wide$date >= startdate & df_wide$date <= enddate,]

# Checks
check <- df_wide[, !grepl('lag|_c|_p|IMD', names(df_wide))]
names(check)

# Expand cols of interest
h_cols <- c('closed', 'total_patients', 'location', h_cols[-1])
w_cols <- c(w_cols[!grepl("^sh_", w_cols)],
            'mx2t_IMD_b', 'mn2t_IMD_b', 'tp_IMD_b', 'fg10_IMD_b', 
            'rh_low_IMD_b', 'rh_hig_IMD_b', 'pm25_IMD_b',
            'mx2t_p90', 'mn2t_p10', 'tp_p90', 'fg10_p90', 
            'rh_mean_p10', 'rh_mean_p90', 'pm25_p90')
location_levels <- unique(df_wide$location)


# MT1: Overall summary statistics (portrait) #####

base <- get_desc(df_wide, h_cols, w_cols[c(1:9)])[c(1:48,58),]
base <- replace_exact_matches(base, nm)

# Format
section_rows    <- which(is.na(base$Frequency) & is.na(base$Average))  # rows 1,4,13,16,22,41
base_display <- base %>%
  mutate(
    Frequency = trimws(Frequency),
    Frequency = gsub("%", "\\%", Frequency, fixed = TRUE),
    across(everything(), escape_latex_unicode),
    across(everything(), ~ifelse(is.na(.), "", .)),
    Variable  = ifelse(row_number() %in% section_rows,
                       paste0("\\textbf{", Variable, "}"),
                       paste0("\\hspace{1em}", Variable))
  )

# Create latex table
kt <- get_kbl_tab(base_display, 
                  title    = "Descriptive statistics of hospital and weather characteristics",
                  note     = "Example",
                  colnames = c("Variable", "Frequency n (\\%)", "Average mean $\\pm$ SD (min, max)"),
                  label    = "MT_base_desc"
) |>
  column_spec(1, width = "5cm", latex_column_spec = "p{5cm}") |>
  column_spec(2, width = "4cm", latex_column_spec = "p{4cm}") |>
  column_spec(3, width = "5cm", latex_column_spec = "p{5cm}") 


# Save
writeLines(
  as.character(kt),
  paste0(path_report, "II. dlnm/MT_base_desc.tex")
)


# A4: Site-specific summary descriptive (landscape) #####

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


# MT2: Meta-regression results ####

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

# Build table

output <- readRDS(paste0('2_DLNM/pooled_model.rds'))

# Define predictors
predictors_list <- list(
  elev = "Altitude",
  popu = "Population density",
  rh   = "Relative humidity",
  pm25   = "Air pollution"
  
)
(tbl_tmp <- make_meta_reg_table(meta=output$meta_tmp, predictors_list))
(tbl_tp  <- make_meta_reg_table(meta=output$meta_tp,  predictors_list))

# Format flextable
d_tmp <- prep_display(tbl_tmp)
d_tp  <- prep_display(tbl_tp)
(combined <- rbind(section_row(), d_tmp, section_row(), d_tp))
tmp_rows   <- 2:6
tp_rows    <- 8:12
label_rows <- c(1, 7)
intercept_rows <- c(2, 8)   # italic "Intercept only" rows



kt <- get_kbl_tab(
  dat     = combined,
  title   = "Meta regression",
  note    = "Example",
  colnames = c("Meta-predictors", "Stat", "df", "p", "AIC",
               "Stat", "df", "p", "I\\textsuperscript{2} (\\%)"),
  label   = "MT_meta_reg",
  spanning_header  = c(" " = 1, "LR test" = 3, " " = 1,
                       "Cochran Q test" = 3, " " = 1),
  italic_rows      = c(2, 8),   # intercept-only rows
  hline_after_rows = c(6)       # thin rule between Temperature/Precipitation
)

# Save
writeLines(
  as.character(kt),
  paste0(path_report, "II. dlnm/MT_meta_reg.tex")
)



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
  compose(i = 7, j = 1, part = "body",
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
  hline(i = 6, part = "body",   border = fp_border(width = 0.75)) |>
  
  padding(padding.top = 2, padding.bottom = 2, part = "all") |>
  set_table_properties(layout = "fixed")

ft

#      Save table ####
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


# ST: RR with CI for all exposure-lags, across locations ####

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

# ST: RR with CI for all exposure-lags, across population groups ####

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



