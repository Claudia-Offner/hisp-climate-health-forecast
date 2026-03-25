library(readr)
library(openxlsx)
library(dplyr)

# Folder
folder_path <- "C:/Users/claer14/Downloads/"
input_folder <- file.path(folder_path, "HISP")

# Output Excel file
output_excel <- "output.xlsx"

# ---- 1. DEFINE COLUMN NAME SETS FOR EACH FILE ----
# Example: adjust to your real structure
column_map <- list(
  "identifiers" = c(
    "patient_identifier_id","patient_id","identifier","identifier_type",
    "preferred","location_id","creator","date_created",
    "voided","voided_by","date_voided","void_reason","uuid"
  ),
  "hisp" = c(
    "ward_name","patient_id","age","gender","category","date","visit_type"
  ),
  "otherfile" = c("col1","col2","col3")   # add more as needed
)
# ---------------------------------------------------

# Create Excel workbook
wb <- createWorkbook()

# List CSVs
csv_files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$")

for (file_path in csv_files) {
  
  # Extract file stem for sheet name + column selection
  file_stem <- tools::file_path_sans_ext(basename(file_path))
  
  # Pick column names (if defined), else read header normally
  col_names_to_use <- column_map[[file_stem]]
  
  message("Processing: ", file_stem)
  
  # ---- 2. READ FILE WITH OPTIONAL CUSTOM COLUMN NAMES ----
  if (!is.null(col_names_to_use)) {
    df <- read_delim(
      file_path,
      delim = "tt",
      col_names = FALSE,   # because you're supplying your own
      escape_double = FALSE,
      trim_ws = TRUE
    )
    names(df) <- col_names_to_use
    
  } else {
    # fallback: use header in CSV
    df <- read_delim(
      file_path,
      delim = "tt",
      col_names = FALSE,
      escape_double = FALSE,
      trim_ws = TRUE
    )
  }
  df <- df[1]
  names(df) <- file_stem
  
  # ---- 3. ADD TO EXCEL ----
  addWorksheet(wb, file_stem)
  writeData(wb, file_stem, df)
}

# Save workbook
saveWorkbook(wb, output_excel, overwrite = TRUE)

cat("All CSVs combined into: ", output_excel, "\n")
