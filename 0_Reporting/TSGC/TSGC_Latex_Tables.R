### ------------------------------------------------------------------------ ### 
###  TSGC Results - LaTex table for-mating                                    ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/III. Tables/')
setwd(path)

####  PACKAGES & FUNCTIONS ####
library(readxl)
library(knitr)
library(kableExtra)
library(dplyr)
library(stringr)

# Helper function to insert a row at a given position
insert_row <- function(df, pos, variable_label) {
  new_row <- tibble(Variable = variable_label, Frequency = NA, Average = NA)
  bind_rows(
    df[1:pos, ],
    new_row,
    df[(pos + 1):nrow(df), ]
  )
}

# Function to replace variable names with set names
replace_exact_matches <- function(df, dict) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      col <- ifelse(col %in% names(dict), dict[col], col)
    }
    col
  })
  df
}


####  MAIN CODE ####

# Load tables ####
files <- list.files(pattern = "\\.xlsx$", ignore.case = TRUE)

for (file in files) {
  # Extract the base name without extension
  var_name <- tools::file_path_sans_ext(basename(file))
  # Read the Excel file and assign variable name
  data <- read_excel(file, .name_repair = "minimal")
  # Add \\ for % so latex does not get confused
  data[] <- lapply(data, function(x) {
    if (is.character(x)) gsub("%", "\\\\%", x) else x
  })
  assign(var_name, data)
}


# ST: ST_comp_model_Apr16  ####

df <- cbind(ST_comp_model_Apr16_noExp, ST_comp_model_Apr16_ALL[-1])
foot <- "Forecasts were conducted with and without weather exposures, for comparison. Lower AIC and BIC values indicate more suitable models and higher LogLik values indicate better model fit. MAPE: Mean Accuracy Percentage Error; AIC: Akaike Information Criterion; BIC: Bayesian Information Criterion."

# Set label & col names
colnames(df) <- paste0("\\textbf{", colnames(df), "}")

# Get latex code
latex_code <- kbl(
  df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  longtable = TRUE,
  caption = NULL
) %>%
  add_header_above(c(" " = 1, "Without weather" = 4, "With weather" = 4)) %>%
  kable_styling(position = "center", latex_options = c("repeat_header", "scale_down"), 
                repeat_header_continued = "\\textit{(continued on next page)}") %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces & repeat titles
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
latex_code <- gsub(paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{\\@\\{\\}l\\}\\{\\\\textit\\{\\(continued\\)\\}\\}\\\\\\\\"), "", latex_code)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.75\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# Export
cat(latex_code, file = "ST_comp_model_Apr16.tex")

# ST: ST_comp_model_Aug16  ####

df <- cbind(ST_comp_model_Aug16_noExp, ST_comp_model_Aug16_ALL[-1])
foot <- "Forecasts were conducted with and without weather exposures, for comparison. Lower AIC and BIC values indicate more suitable models and higher LogLik values indicate better model fit. MAPE: Mean Accuracy Percentage Error; AIC: Akaike Information Criterion; BIC: Bayesian Information Criterion."

# Set label & col names
colnames(df) <- paste0("\\textbf{", colnames(df), "}")

# Get latex code
latex_code <- kbl(
  df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  longtable = TRUE,
  caption = NULL) %>%
  add_header_above(c(" " = 1, "Without weather" = 4, "With weather" = 4)) %>%
  kable_styling(position = "center", latex_options = c("repeat_header", "scale_down"), 
                repeat_header_continued = "\\textit{(continued on next page)}") %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces & repeat titles
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
latex_code <- gsub(paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{\\@\\{\\}l\\}\\{\\\\textit\\{\\(continued\\)\\}\\}\\\\\\\\"), "", latex_code)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.75\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# latex_code <- gsub(
#   "(\\\\endfirsthead\\s*)\\\\caption\\[\\]\\{.*?\\}\\\\\\\\\\s*\\\\toprule",
#   "\\1\\\\toprule",
#   latex_code,
#   perl = TRUE
# )

# Export
cat(latex_code, file = "ST_comp_model_Aug16.tex")

# ST: ST_sens_lag ####

df <- ST_sens_lag
foot <- "Tests were run for lags 0 to 30 days for all weather exposures, across both IMD and Climatology (Clim) thresholds. Checks were conducted on two hospital outcomes to enable generalisation. Lower AIC and BIC values indicate more suitable models and higher LogLik values indicate better model fit. MAPE: Mean Accuracy Percentage Error; AIC: Akaike Information Criterion; BIC: Bayesian Information Criterion."

# Set label & col names
colnames(df) <- paste0("\\textit{", colnames(df), "}")

# Get latex code
latex_code <- kbl(
  df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  longtable = TRUE,
  caption = NULL) %>%
  add_header_above(c(" " = 2, "Outcome: Total patients" = 4, "Outcome: General medicine" = 4)) %>%
  kable_styling(position = "center", latex_options = c("repeat_header", "scale_down"),
                repeat_header_continued = "\\textit{(continued on next page)}") %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces & repeat titles
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
latex_code <- gsub(paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{\\@\\{\\}l\\}\\{\\\\textit\\{\\(continued\\)\\}\\}\\\\\\\\"), "", latex_code)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.85\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# Export
cat(latex_code, file = "ST_sens_lag.tex")

# ST: ST_sens_param  ####

df <- ST_sens_param
foot <- "Tests were run for q values (range: 0 to 0.10), estimation window (range: 1 to 12 months), forecast horizons (range: 1 to 20 days) and both IMD and Climatology (Clim) thresholds. Checks were conducted on two hospital outcomes to enable generalisation. Lower AIC and BIC values indicate more suitable models and higher LogLik values indicate better model fit. MAPE: Mean Accuracy Percentage Error; AIC: Akaike Information Criterion; BIC: Bayesian Information Criterion."

# Set label & col names
colnames(df) <- paste0("\\textit{", colnames(df), "}")

# Get latex code
latex_code <- kbl(
  df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  longtable = TRUE,
  caption = NULL) %>%
  add_header_above(c(" " = 2, "Outcome: Total patients" = 4, "Outcome: General medicine" = 4)) %>%
  kable_styling(position = "center", latex_options = c("repeat_header", "scale_down"), 
                repeat_header_continued = "\\textit{(continued on next page)}") %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces & repeat titles & continues
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
latex_code <- gsub(paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{\\@\\{\\}l\\}\\{\\\\textit\\{\\(continued\\)\\}\\}\\\\\\\\"), "", latex_code)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.85\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# Export
cat(latex_code, file = "ST_sens_param.tex")


# ST: ST_time_desc ####

df <- ST_time_desc
foot <- 'Continuous variables report mean and SD, while binary variables report percentage of occurrence, within each respective grouping. Tp values are reported with positive SDs only because of its inverse exponential distribution. Sundays were NA because Chamba hospital was closed on these days. IMD: Indian Meteorological Department thresholds; Clim: Climatology thresholds; SD: Standard deviation; Mx2t: Maximum 2-meter surface temperature; Mn2t: Minimum 2-meter surface temperature; Tp: Total precipitation; Fg10: Maximum wind gust at 10-meters.'

# Set label and format column names
colnames(df)[1] <- 'Group'
colnames(df) <- paste0("\\textbf{", str_remove_all(colnames(df), "\\\\textbf\\{|\\}"), "}")

# Annotate first row
for (col in names(df)) {
  if (any(grepl("%", df[[col]]))) {
    df[1, col] <- "\\textit{\\%}"
  } else if (any(grepl("±", df[[col]]))) {
    df[1, col] <- "\\textit{$\\mu \\pm$ SD}"
  }
}

# Split wide table
rows <- df[1:2]
t1 <- cbind(rows, df[3:10])
t2 <- cbind(rows, df[11:17])
t3 <- cbind(rows, df[18:25])
t4 <- cbind(rows, df[26:31])
t5 <- cbind(rows, df[32:35])

tables <- list(t1=t1, t2=t2, t3=t3, t4=t4, t5=t5)
for (i in names(tables)) {

    if(i =="t1"){
      # Generate LaTeX table
      latex_code <- kbl(tables[[i]], format = "latex", booktabs = TRUE, escape = FALSE, longtable = FALSE, 
                        # caption = paste0('Summary statistics of included variables (2015–2020), disaggregated by year, month and weekday \\label{App:', label, '}'),
                        align = c("l", rep("c", ncol(df) - 1))) %>%
        kable_styling(position = "left", latex_options = c("scale_down"))
      
    } else{
      
      # Generate LaTeX table
      latex_code <- kbl(tables[[i]], format = "latex", booktabs = TRUE, escape = FALSE, longtable = FALSE, 
                        caption = NULL,
                        align = c("l", rep("c", ncol(df) - 1))) %>%
        kable_styling(position = "left", latex_options = c("scale_down"))
    }

    if (i != "t5") {
      latex_code <- latex_code %>%
        footnote(
          general = "(continued on next page)",
          general_title = "",
          footnote_as_chunk = TRUE
        )
      latex_code <- gsub(
        "\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\(continued on next page\\)\\}\\\\\\\\",
        "{r}{\\\\rule{0pt}{1em}\\\\textit{(continued on next page)}}\\\\",
        latex_code
      )
    } else {
      latex_code <- latex_code %>%
        footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)
      # Add foot note formatting
      latex_code <- gsub(
        paste0("\\\\multicolumn\\{", ncol(tables[[i]]), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
        paste0("\\\\multicolumn{", ncol(tables[[i]]), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.95\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
        latex_code
      )
      latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)
    }
    
    
    # Clean up code
    latex_code <- gsub("\n\\\\addlinespace", "", latex_code)
    latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
    latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
    
    # Export
    cat(latex_code,file = paste0("ST_time_desc_", i, ".tex"))
    }

# ST: ST_samp_desc ####

df <- cbind(ST_apr_desc[,-4], ST_aug_desc[,-c(1,4)])
foot <- "* Binary event percentages represent the proportion of occurrence over time and hospital admission percentages represent the proportion of total patients. Tp values are reported with positive SDs only because of its inverse exponential distribution. IMD: Indian Meteorological Department thresholds; Clim: Climatology thresholds; SD: Standard deviation; Mx2t: Maximum 2-meter surface temperature; Mn2t: Minimum 2-meter surface temperature; Tp: Total precipitation; Fg10: Maximum wind gust at 10-meters."

# Add indentation
no_indent_1x <- c(1, 26)
no_indent_2x <- c(2:4, 7, 13, 27:38)  # These get at most one indent
df[["Variable"]] <- dplyr::case_when(
  seq_len(nrow(df)) %in% no_indent_1x ~ paste0("\\textbf{", df[["Variable"]], "}"),
  seq_len(nrow(df)) %in% no_indent_2x ~ paste0("\\hspace{1em}", df[["Variable"]]),
  TRUE ~ paste0("\\hspace{2em}\\textit{", df[["Variable"]], "}")
)

# Remove the NA values
df[is.na(df)] <- " "
# Add headings
names(df) <- c("\\textbf{Variable}",  
               "\\textbf{Frequency*} (\\textit{\\%})",
               "\\textbf{Average} (\\textit{$\\pm$SD})", 
               "\\textbf{Frequency*} (\\textit{\\%})",
               "\\textbf{Average} (\\textit{$\\pm$SD})")
# Get latex code
latex_code <- kable(df, "latex", booktabs = TRUE, escape = FALSE, longtable = TRUE, 
                    # caption = paste0("Summary statistics of all variables for the 3-month estimation period prior to each sample period, 01 April 2016 and 14 August 2016. \\label{tab:", label, "}"),
                    align = c("l", rep("c", ncol(df) - 1))) %>%
  add_header_above(c(" " = 1, "01 April 2016" = 2, "14 August 2016" = 2)) %>%
  kable_styling(position = "center", latex_options = c("repeat_header", "scale_down"),
                repeat_header_continued = "\\textit{(continued on next page)}") %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
latex_code <- gsub(paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{\\@\\{\\}l\\}\\{\\\\textit\\{\\(continued\\)\\}\\}\\\\\\\\"), "", latex_code)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.95\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# Export
cat(latex_code, file = "ST_samp_desc.tex")

# MT: MT_base_desc ####

df <- MT_base_desc
label <- 'MT_base_desc'
foot <- "* Binary event percentages represent the proportion of occurrence over time and hospital admission percentages represent the proportion of total patients. Tp values are reported with positive SDs only because of its inverse exponential distribution. IMD: Indian Meteorological Department thresholds; Clim: Climatology thresholds; SD: Standard deviation; Mx2t: Maximum 2-meter surface temperature; Mn2t: Minimum 2-meter surface temperature; Tp: Total precipitation; Fg10: Maximum wind gust at 10-meters."
# Add indentation
no_indent_1x <- c(1, 26)
no_indent_2x <- c(2:4, 7, 13, 27:38)  # These get at most one indent
df <- df %>%
  dplyr::mutate(Variable = case_when(
    row_number() %in% no_indent_1x ~ paste0("\\textbf{", Variable, "}"),  # No indent + bold
    row_number() %in% no_indent_2x ~ paste0("\\hspace{1em}", Variable),  # Single indent + italics
    TRUE ~ paste0("\\hspace{2em}\\textit{", Variable, "}")  # Double indent (normal font)
  ))
# Remove the NA values
df[is.na(df)] <- " "
# Add headings
names(df) <- c("\\textbf{Variable}",  
               "\\textbf{Frequency*} (\\textit{\\%})",
               "\\textbf{Average} (\\textit{$\\pm$SD})", 
               "\\textbf{Range} (\\textit{Skew})"
)
# Get latex code
latex_code <- kable(df, "latex", booktabs = TRUE, escape = FALSE,
                    caption = paste0("Summary statistics of variables used in analysis (2015-2020) \\label{tab:", label, "}"),
                    align = c("l", rep("c", ncol(df) - 1))) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.95\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)


cat(latex_code, file = "MT_base_desc.tex")


# MT: MT_weat_coefs_Apr16 ####

df <- MT_weat_coefs_Apr16_ALL
label <- 'MT_weat_coefs_Apr16_ALL'
foot <- "Results are reported as the percentage change in respective hospital admission growth rate. All weather exposues are provided at lag 1. Different TSGC models were run for each hospital outcome. Significant results are reported in bold. Maximum temperature and total precipitation were excluded since they did not occur over this period. CI: Confidence Interval; Mn2t: Minimum 2-meter surface temperature; Fg10: Maximum wind gust at 10-meters."

# BOLD SIGNIFICANT RESULTS
df[1:7, 2] <- paste0("\\textbf{", unlist(df[1:7, 2]), "}") # Most mn2t
df[9, 2] <- paste0("\\textbf{", unlist(df[9, 2]), "}") # general medicine
df[11:12, 2] <- paste0("\\textbf{", unlist(df[11:12, 2]), "}") # general medicine
df[3, 3] <- paste0("\\textbf{", unlist(df[3, 3]), "}") # Fg10 Men
df[7, 3] <- paste0("\\textbf{", unlist(df[7, 3]), "}") # Fg10 age 25-64

# Get latex code
latex_code <- kable(df, "latex", booktabs = TRUE, escape = FALSE,
                    caption = paste0("Weather coefficients from each hospital registration model, for the period commencing on 01 April 2016. \\label{tab:", label, "}"),
                    align = c("l", rep("c", ncol(df) - 1))) %>%
  # add_header_above(c(" " = 1, "IMD: Mn2t <0\u00B0C" = 2, "IMD: Fg10 >40km/h" = 2, " "=1)) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.95\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# Export
cat(latex_code, file = "MT_weat_coefs_Apr16.tex")


# MT: MT_weat_coefs_Aug16 ####

df <- MT_weat_coefs_Aug16_ALL
label <- 'MT_weat_coefs_Aug16_ALL'
foot <- "Results are reported as the percentage change in respective hospital admission growth rate. All weather exposues are provided at lag 1. Different TSGC models were run for each hospital outcome. Significant results are reported in bold. No weather events occurred over this period. CI: Confidence Interval; Mn2t: Minimum 2-meter surface temperature; Fg10: Maximum wind gust at 10-meters."

# BOLD SIGNIFICANT RESULTS
df[6, 2] <- paste0("\\textbf{", unlist(df[6, 2]), "}") # ages 15-24
df[6, 4] <- paste0("\\textbf{", unlist(df[6, 4]), "}") # ages 15-24
df[1:11, 3] <- paste0("\\textbf{", unlist(df[1:11, 3]), "}") # Most mn2t

# Get latex code
latex_code <- kable(df, "latex", booktabs = TRUE, escape = FALSE, #align = rep("c", ncol(df)),
                    caption = paste0("Weather coefficients from each hospital registration model, for the period commencing on 14 August 2016. \\label{tab:", label, "}"),
                    align = c("l", rep("c", ncol(df) - 1))) %>%
  # add_header_above(c(" " = 1, "IMD: Mn2t <0\u00B0C" = 2, "IMD: Fg10 >40km/h" = 2, " "=1)) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  footnote(general = paste0(foot, "}"), general_title = "Note:", footnote_as_chunk = TRUE)

# Remove periodic line spaces
latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
latex_code <- gsub(">", "$>$", latex_code, fixed = TRUE)
latex_code <- gsub("<", "$<$", latex_code, fixed = TRUE)
# Add foot note formatting
latex_code <- gsub(
  paste0("\\\\multicolumn\\{", ncol(df), "\\}\\{l\\}\\{\\\\rule\\{0pt\\}\\{1em\\}\\\\textit\\{Note:\\}"),
  paste0("\\\\multicolumn{", ncol(df), "}{l}{\\\\rule{0pt}{1em}\\\\parbox[t]{0.95\\\\linewidth}{\\\\fontsize{9}{11}\\\\selectfont\\\\textit{Note:}"),
  latex_code
)
latex_code <- gsub("\\\\\\}\\}\\\\\\\\", "}}\\\\\\\\", latex_code)

# Export
cat(latex_code, file = "MT_weat_coefs_Aug16.tex")

################################################################################

# ST: ST_ttest

# # Add headings
# names(ST_ttest)[2:3] <- c('Group 0 (\\textit{$\\mu})', 'Group 1 (\\textit{$\\mu})')
# colnames(ST_ttest) <- paste0("\\textbf{", colnames(ST_ttest), "}")
# # Get latex code
# latex_code <- kable(ST_ttest, "latex", booktabs = TRUE, escape = FALSE,
#                     caption = "A LaTeX table") %>%
#   kable_styling(latex_options = c("hold_position")) %>%
#   footnote(general = "This is a footnote for the table.",
#            general_title = "Note:", 
#            footnote_as_chunk = TRUE)
# # Remove periodic line spaces
# latex_code <- gsub("\n\\\\addlinespace", "", latex_code) 
# # Export
# cat(latex_code, file = "ST_ttest.tex")