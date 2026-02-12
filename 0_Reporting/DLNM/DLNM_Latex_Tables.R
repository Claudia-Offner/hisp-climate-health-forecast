### ------------------------------------------------------------------------ ### 
###  DLNM Results - LaTex table formatting                                    ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 2 - Causal Inference/III. Tables/')
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


# Load tables
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

# T: MT_base_desc ####

df <- MT_base_desc
label <- 'MT_base_desc'
foot <- "* Binary event percentages represent the proportion of occurrence over time and hospital admission percentages represent the proportion of total patients. Tp values are reported with positive SDs only because of its inverse exponential distribution. SD: Standard deviation; T2m: Mean 2-meter surface temperature; Mx2t: Maximum 2-meter surface temperature; Mn2t: Minimum 2-meter surface temperature; Tp: Total precipitation; Fg10: Maximum wind gust at 10-meters."
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

# Save table format
cat(latex_code, file = "MT_base_desc.tex")


# ST: ST_time_desc ####

df <- ST_time_desc
foot <- 'Continuous variables report mean and SD, while binary variables report percentage of occurrence, within each respective grouping. Tp values are reported with positive SDs only because of its inverse exponential distribution. Sundays were NA because Chamba hospital was closed on these days. SD: Standard deviation; T2m: Mean 2-meter surface temperature; Mx2t: Maximum 2-meter surface temperature; Mn2t: Minimum 2-meter surface temperature; Tp: Total precipitation; Fg10: Maximum wind gust at 10-meters.'

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
t4 <- cbind(rows, df[26:28])
tables <- list(t1=t1, t2=t2, t3=t3, t4=t4)

for (i in names(tables)) {
  
  # Generate LaTeX table
  latex_code <- kbl(tables[[i]], format = "latex", booktabs = TRUE, escape = FALSE, longtable = FALSE, caption = NULL) %>%
    kable_styling(position = "left", latex_options = c("scale_down"))
  
  if (i != "t4") {
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
  
  # Replace the \ifdim ... with \maxwidth
  latex_code_fixed <- gsub(
    "\\\\resizebox\\{\\\\ifdim\\\\width\\$>$\\\\linewidth.*?\\}\\{!\\}\\{",
    "\\\\resizebox{\\\\maxwidth}{!}{",
    latex_code
  )
  
  # Export
  cat(latex_code,file = paste0("ST_time_desc_", i, ".tex"))
}


