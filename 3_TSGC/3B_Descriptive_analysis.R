### ------------------------------------------------------------------------ ### 
###  TSGC Multivariate Model - Descriptive analysis                          ###
### ------------------------------------------------------------------------ ### 

####  DEPENDENCIES  #####

# Packages
library(openxlsx); library(patchwork); library(moments)

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
file_import <- 'ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_wide'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/')
setwd(path1)

# Variables
vars <- c('closed', 'total_patients', 
          'F','M','under_5','age5_14','age15_24','age25_64','over_65',
          'general','ophthalmology','pediatrics','orthopedics','ent','gynecology',
          'dermatology','dentistry','surgery','emergency', 'psychiatry','oncology',
          'mx2t', 'mn2t', 'tp', 'fg10',
          'mx2t_IMD_b', 'mn2t_IMD_b', 'tp_IMD_b', 'fg10_IMD_b',
          'mx2t_p90', 'mn2t_p10', 'tp_p90', 'fg10_p90')


### Get model functions (make sure wd is in repo folder)
source(paste0(path1, '3_TSGC/TSGC_Model - Functions.R'))

### Function to summarize by a grouping variable
summarise_by_group <- function(group_var) {
  map_dfr(vars, function(var) {
    x <- df[[var]]
    is_binary <- all(x %in% c(0, 1, NA), na.rm = TRUE)
    
    df %>%
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

### Function to get corrplots
corrplotting <- function(variable, y_lab){
  
  
  p1 <- ggplot(df, aes(x = mx2t, y = {{ variable }}, color=season)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "(A) Mx2t",
      x = "°C",
      y = y_lab,
      color = "Season"
    ) +
    scale_color_manual(values = c('Winter (Jan-Feb)'=mypalette[3],
                                  'Pre-Monsoon (Mar-May)'=mypalette[4], 
                                  'Monsoon (Jun-Sep)'=mypalette[2], 
                                  'Post-Monsoon (Oct-Dec)'=mypalette[1])) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  
  p2 <- ggplot(df, aes(x = mn2t, y = {{ variable }}, color=season)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "(B) Mn2t",
      x = "°C",
      y = y_lab,
      color = "Season"
    ) +
    scale_color_manual(values = c('Winter (Jan-Feb)'=mypalette[3],
                                  'Pre-Monsoon (Mar-May)'=mypalette[4], 
                                  'Monsoon (Jun-Sep)'=mypalette[2], 
                                  'Post-Monsoon (Oct-Dec)'=mypalette[1])) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  p3 <- ggplot(df, aes(x = tp, y = {{ variable }}, color=season)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "(C) Tp",
      x = "mm",
      y = y_lab,
      color = "Season"
    ) +
    scale_color_manual(values = c('Winter (Jan-Feb)'=mypalette[3],
                                  'Pre-Monsoon (Mar-May)'=mypalette[4], 
                                  'Monsoon (Jun-Sep)'=mypalette[2], 
                                  'Post-Monsoon (Oct-Dec)'=mypalette[1])) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  p4 <- ggplot(df, aes(x = fg10, y = {{ variable }}, color=season)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "(D) Fg10",
      x = "km per hour",
      y = y_lab,
      color = "Season"
    ) +
    scale_color_manual(values = c('Winter (Jan-Feb)'=mypalette[3],
                                  'Pre-Monsoon (Mar-May)'=mypalette[4], 
                                  'Monsoon (Jun-Sep)'=mypalette[2], 
                                  'Post-Monsoon (Oct-Dec)'=mypalette[1])) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  corrplot <- (p1 | p2) / (p3 | p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  return(corrplot)
  
}

### Base statistics functions
get_desc <- function(df, vars) {
  
  # Get descriptive stats 
  (desc <- vars %>%
     lapply(function(col) {
       x <- df[[col]]
       
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
            ifelse(Variable == "tp", " + ", " ± "),
            round(sd_daily_visits)
          )
        ),
        Range = ifelse(
          is.na(skew_daily_visit),
          NA,
          paste0(min_daily_visit, ', ',max_daily_visit, ' (', skew_daily_visit, ')')
        )
      ) %>%
      dplyr::select(Variable, Frequency, Average, Range))
  
  
  # Remove frequency values for raw weather data
  base_desc$Frequency[22:25] <- NA  # or "" if you prefer blank strings
  
  # Get percentage of hospital registrations
  # perc <- as.numeric(gsub(",", "", base_desc$Frequency))/as.numeric(gsub(",", "", base_desc[1, "Frequency"])) *100
  # perc <- cbind(base_desc$Variable, perc)
  
  # FORMATING
  
  # Insert headings at fixed positions
  base_desc <- rbind(
    tibble(Variable = "Hospital Characteristics", Frequency = NA, Average = NA, Range = NA),
    base_desc
  )
  positions <- c(3, 6, 12, 25)
  labels <- c("Gender", "Age", "Ward", "Weather Characteristics")
  for (i in seq_along(positions)) {
    pos <- positions[i]
    label <- labels[i]
    new_row <- tibble(Variable = label, Frequency = NA, Average = NA, Range = NA)
    base_desc <- rbind(
      base_desc[1:pos, ],
      new_row,
      base_desc[(pos + 1):nrow(base_desc), ]
    )
  }
  
  return(base_desc)
}

####  MAIN CODE  #####

# Load data
load(paste0('3_TSGC/', file_import, '.RData'))

# 1. Get time series visual (Total patients over time period: 2015-2016) ####

# Get moving average (no NAs)
clean <- df[!is.na(df$total_patients), ]
# Select year of interest
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2017-02-13')
clean <- clean[clean$date >= startdate & clean$date <= enddate,]
variable <- xts(clean$total_patients,  order.by=clean$date)
ts_plot <- plot_moving_avg(variable, MA=7, y.lab='Total patients', colour='#8968CD',
                               plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)

# GET STACKED TS
variable <- xts(clean$total_patients,  order.by=clean$date)
ts_plot_tot <- plot_moving_avg(variable, MA=7, y.lab='Patients', colour='#8968CD',
                               plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
variable <- xts(clean$mx2t,  order.by=clean$date)
ts_plot_mx2t <- plot_moving_avg(variable, MA=7, y.lab='Celsius',colour='#D74B3C',
                                plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

variable <- xts(clean$mn2t,  order.by=clean$date)
ts_plot_mn2t <- plot_moving_avg(variable, MA=7, y.lab='Celsius', colour='#1D98EB',
                                plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

variable <- xts(clean$tp,  order.by=clean$date)
ts_plot_tp <- plot_moving_avg(variable, MA=7, y.lab='mm', colour='#A6B727',
                              plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

variable <- xts(clean$fg10,  order.by=clean$date)
ts_plot_fg10 <- plot_moving_avg(variable, MA=7, y.lab='km/hour', colour='#FEC306',
                                plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank())

# GET LEGEND ONLY
legend_data <- data.frame(
  x = 1:5, y = 1:5,
  group = c("Total Patients", "Mx2t", "Mn2t", "Tp", "Fg10"))
legend_data$group <- factor(legend_data$group, levels = c("Total Patients", "Mx2t", "Mn2t", "Tp", "Fg10"))
color_vals <- c("Total Patients" = "#8968CD",
                "Mx2t" = "#D74B3C",
                "Mn2t" = "#1D98EB",
                "Tp" = "#A6B727",
                "Fg10" = "#FEC306")
manual_legend <- ggplot(legend_data, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = color_vals, name = NULL) +
  theme_void() +
  theme(legend.position = "bottom")
# Display only the legend
legend_only <- manual_legend + guides(color = guide_legend(nrow = 1)) + theme(legend.box.margin = margin(t = -20))


# Stack the plots with relative heights and insert legend at the bottom
ts_plot_stack <- (
  ts_plot_tot / ts_plot_mx2t / ts_plot_mn2t / ts_plot_tp / ts_plot_fg10 / legend_only)+ 
  plot_layout(heights = c(1, 1, 1, 1, 1, 0.15))


# 2. Get baseline summary statistics (2015-2020) ####

base_desc <- get_desc(df, vars)

# 3. Get sample specific summary statistics ####

enddate <- as.Date('2016-04-01')
startdate <- enddate %m-% months(3)
apr <- df[df$date >= startdate & df$date <= enddate, ]

enddate <- as.Date('2016-08-14')
startdate <- enddate %m-% months(3)
aug <- df[df$date >= startdate & df$date <= enddate, ]

apr_desc <- get_desc(apr, vars)
aug_desc <- get_desc(aug, vars)

# Run t-tests ####

# # exposures <- c('mx2t_p90', 'mn2t_p10', 'tp_p90', 'fg10_p90')
# # exposures <- c('mn2t_IMD_b', 'fg10_IMD_b')
# 
# ttest <- data.frame()
# for (i in exposures){
#   x <- t.test(df$total_patients ~ df[[i]], data = df)
#   table_row <- data.frame(
#     variable = i,
#     mean_group_0 = x$estimate[1],
#     mean_group_1 = x$estimate[2],
#     mean_diff = diff(x$estimate),
#     t_stat = x$statistic,
#     df = x$parameter,
#     p_value = x$p.value,
#     ci_lower = x$conf.int[1],
#     ci_upper = x$conf.int[2]
#   )
#   ttest <- rbind(ttest, table_row)
# }
# ttest[ , -1] <- round(ttest[ , -1], 2)  # assuming 'variable' is the only non-numeric column
# rownames(ttest) <- NULL
# names(ttest) <- c('Variable', 'Group 0', 'Group 1', 'Diff.', 
#                   't', 'df', 'p-value', 'Lower CI', 'Upper CI')

#### EXPORT RESULTS ####

# ttest <- replace_exact_matches(ttest, nm)
base_desc <- replace_exact_matches(base_desc, nm)
apr_desc <- replace_exact_matches(apr_desc, nm)
aug_desc <- replace_exact_matches(aug_desc, nm)
# names(time_desc) <- gsub("_b", "", names(time_desc))
# names(time_desc) <- replace_exact_matches(names(time_desc), nm)

# Save tables
# write.xlsx(ttest, file = paste0(path3, "/III. Tables/ST_ttest.xlsx"))
write.xlsx(base_desc, file = paste0(path3, "/III. Tables/base_desc.xlsx"))
write.xlsx(apr_desc, file = paste0(path3, "/III. Tables/MT_apr_desc.xlsx"))
write.xlsx(aug_desc, file = paste0(path3, "/III. Tables/MT_aug_desc.xlsx"))

# Save figures
ggsave(paste0(path3, '/IV. Figures/MF_ts_nolab.png'), ts_plot, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/SF_ts_stack.png'), ts_plot_stack, width = 8, height = 8, dpi = 300)
