### ---------------------------------- ### 
### Web scrapper for regional holidays ###
### ---------------------------------- ### 


# Set data path
user <- 'claer14'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')

# Packages
library(rvest)
library(readr)
library(dplyr)
library(lubridate)

# Extract the holiday table
years <- as.character(seq(2015, 2023))
holiday_table <- data.frame()

for(y in years){
  
  url <- paste0("https://www.officeholidays.com/countries/india/himachal-pradesh/", y)
  
  # Read the HTML content
  page <- read_html(url)
  
  # Extract the holiday table
  x <- page %>%
    html_element("table.country-table") %>%
    html_table() %>%
    mutate(Date = as.Date(paste(Date, y), format = "%b %d %Y"))
  
  holiday_table <- rbind(holiday_table, x)
  
}

# Save to CSV
write_csv(holiday_table, "204_demo/ind_Adm2_OfficeHolidays_2015-2023_HP_holidays_daily.csv")

# Preview the data
print(head(holiday_table))
