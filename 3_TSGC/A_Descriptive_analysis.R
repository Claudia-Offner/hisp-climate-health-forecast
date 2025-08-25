### ------------------------------------------------------------------------ ### 
###  TSGC Multivariate Model - Descriptive analysis                          ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'offne'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/HISP_TSGC/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/')
setwd(path2)

####  PACKAGES & FUNCTIONS ####
library(openxlsx)
library(patchwork)
library(moments)

# Get functions (make sure wd is in repo folder)
source(paste0(path1, 'TSGC_Model - Functions.R'))

# Function to summarize by a grouping variable
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

# Function to get corrplots
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

# Base statistics functions
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


####  DEPENDENCIES  #####
vars <- c('closed', 'total_patients', 
          'F','M','under_5','age5_14','age15_24','age25_64','over_65',
          'general','ophthalmology','pediatrics','orthopedics','ent','gynecology',
          'dermatology','dentistry','surgery','emergency', 'psychiatry','oncology',
          'mx2t', 'mn2t', 'tp', 'fg10',
          'mx2t_IMD_b', 'mn2t_IMD_b', 'tp_IMD_b', 'fg10_IMD_b') #,
          # 'mx2t_p90', 'mn2t_p10', 'tp_p90', 'fg10_p90')


####  MAIN CODE  #####


# 1. Import data ####

# Load data
# df <- read.csv('203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')
df <- read.csv('203_comb/ind_10km_HISP-HP_2015-2017_hosp-bigthaw_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')

# Set date type, 
df$date <- as.Date(df$date, format='%Y-%m-%d')  # Adjust format as needed

# Set factor data
df <- df %>%
  mutate(weekday=factor(weekday, levels=weekday_levels),
         month=factor(month, levels=month_levels), 
         season=factor(season, levels=season_levels),
         year=factor(year))

# Select year of interest
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2017-01-01')
# enddate <- as.Date('2020-02-13')
df <- df[df$date >= startdate & df$date <= enddate,]

# Back transform tp variables
df$tp <- exp(df$tp)

# 2. Create visual for weather-health relationships ####

# Define the inverse exponential function
exponential <- function(x, a = 1, b = 0.1, c = 0.1) {
  return(a * exp(-b * x) + c)
}


# First Graph: Normal Distribution 
x1 <- seq(-4, 4, length.out = 1000)
y1 <- dnorm(x1, mean = 0, sd = 1.25)
p5 <- qnorm(0.01, mean = 0, sd = 1)
p95 <- qnorm(0.99, mean = 0, sd = 1)
data1 <- data.frame(x1, y1)

# Add a quadratic variable to the dataset
# data1$quadratic <- (0.3 * data1$x1^2 + 0.4 * data1$x1) / 20
data1$mortality <- (0.2 * data1$x1^2) / 20  # Adjust coefficient for desired curvature

plot1 <- ggplot(data1, aes(x = x1, y = y1)) +
  geom_ribbon(data = subset(data1, x1 <= p5), aes(ymin = 0, ymax = y1), fill = "skyblue2", alpha = 0.75) +
  geom_ribbon(data = subset(data1, x1 >= p95), aes(ymin = 0, ymax = y1), fill = "indianred2", alpha = 0.75) +
  geom_line(color = "black", size = 0.75) +
  scale_x_continuous(breaks = round(c(p5, 0, p95), 2), labels = c("Cold", "Average", "Hot")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "(a) Temperature", x = 'Temperature frequency', y = "RR of mortality") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = c("skyblue3", "dimgrey", "indianred3"), size = 11),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  annotate("segment", x = p5, xend = p5, y = 0, yend = max(subset(data1, x1 <= p5)$y1), color = "black", linetype = "dashed", size = 0.5) + 
  annotate("segment", x = p95, xend = p95, y = 0, yend = max(subset(data1, x1 >= p95)$y1), color = "black", linetype = "dashed", size = 0.5)+
  # Add quadratic curve here
  geom_line(aes(y = mortality), color = "brown3", size = 1, linetype = "solid")  # Add asymmetric quadratic line


# Second Graph: Inverse Exponential Distribution
a <- 4.2 # initial steepness
b <- 1.2 # rate of decay
c <- 0.1 # plateau point
x2 <- seq(0, 4, length.out = 1000)
y2 <- exponential(x2, a=a, b=b, c=c)
p95_exp2 <- qexp(0.94, rate=b)
data2 <- data.frame(x2, y2)

data2$mortality <- exponential(x2, a=a, b=-1.2, c=c)/200

plot2 <- ggplot(data2, aes(x = x2, y = y2)) +
  geom_ribbon(data = subset(data2, x2 >= p95_exp2), aes(ymin = 0, ymax = y2), fill = "palegreen3", alpha = 1) +
  geom_line(color = "black", size = 0.75) +
  scale_x_continuous(breaks = round(c(p95_exp2-2, p95_exp2-1, p95_exp2), 2), labels = c("Light", "Average", "Heavy")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "(b) Precipitation", x = 'Rainfall frequency', y = "RR of mortality") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = c("darkgoldenrod3", "dimgrey", "forestgreen"), size = 11),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  annotate("segment", x = p95_exp2, xend = p95_exp2, y = 0, yend = max(subset(data2, x2 >= p95_exp2)$y2), color = "black", linetype = "dashed", size = 0.5)+
  geom_line(aes(y = mortality), color = "forestgreen", size = 1, linetype = "solid") +
  coord_cartesian(xlim = c(0, p95_exp2+0.5), ylim = c(0, 4))  # Set both x and y axis limits

# Second Graph: Inverse Exponential Distribution
a <- 4.2 # initial steepness
b <- 1.2 # rate of decay
c <- 0.1 # plateau point
x3 <- seq(0, 4, length.out = 1000)
y3 <- exponential(x3, a=a, b=b, c=c)
p95_exp3 <- qexp(0.94, rate=b)
data3 <- data.frame(x3, y3)

data3$mortality <- exponential(x3, a=a, b=-1.2, c=c)/200

plot3 <- ggplot(data3, aes(x = x3, y = y3)) +
  geom_ribbon(data = subset(data3, x3 >= p95_exp3), aes(ymin = 0, ymax = y3), fill = "mediumpurple1", alpha = 1) +
  geom_line(color = "black", size = 0.75) +
  scale_x_continuous(breaks = round(c(p95_exp3-2, p95_exp3-1, p95_exp3), 2), labels = c("Slow", "Average", "Fast")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "(c) Wind speed", x = 'Wind speed frequency', y = "RR of mortality") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = c("palevioletred3", "dimgrey", "darkorchid4"), size = 11),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  annotate("segment", x = p95_exp3, xend = p95_exp3, y = 0, yend = max(subset(data3, x3 >= p95_exp3)$y3), color = "black", linetype = "dashed", size = 0.5)+
  geom_line(aes(y = mortality), color = "darkorchid4", size = 1, linetype = "solid") +
  coord_cartesian(xlim = c(0, p95_exp3+0.5), ylim = c(0, 4))


# Display the plots
# grid.arrange(plot1, plot2, nrow = 1)
wh_dist <- grid.arrange(plot1, plot2, nrow = 1, widths = c(1, 1), heights = c(1))


# 3. Get time series visual (Total patients over time period: 2015-2020) ####

# Get moving average (no NAs)
clean <- df[!is.na(df$total_patients), ]
# Select year of interest
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


# 4. Get corplot visual (Total patients over time period: 2015-2020) ####
# general; ophthalmology; pediatrics; ent

corrplot_total <- corrplotting(total_patients, 'Total patients')
corrplot_general <- corrplotting(general, 'General Medicine')
corrplot_opth <- corrplotting(ophthalmology, 'Ophthalmology')
corrplot_paed <- corrplotting(pediatrics, 'Paediatrics')
corrplot_ent <- corrplotting(ent, 'ENT')

# 5. Get baseline summary statistics (2015-2020) ####

base_desc <- get_desc(df, vars)


# 6. Get temporal summary statistics ####

# Get summaries for each time grouping
summary_weekday <- summarise_by_group('weekday') %>% mutate(Time_Group = 'Weekday')
summary_month <- summarise_by_group('month') %>% mutate(Time_Group = 'Month')
summary_year <- summarise_by_group('year') %>% mutate(Time_Group = 'Year')

# Combine all summaries and pivot
time_desc <- bind_rows(summary_year, summary_month, summary_weekday) %>%
  dplyr::select(Time_Group, Time, variable, Summary) %>%
  pivot_wider(names_from = variable, values_from = Summary) %>%
  mutate(Time_Group = factor(Time_Group, levels = c("Year", "Month", "Weekday"))) %>%
  arrange(Time_Group, Time)

# 7. Get sample specific summary statistics ####

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
names(time_desc) <- gsub("_b", "", names(time_desc))
names(time_desc) <- replace_exact_matches(names(time_desc), nm)

# Save tables
# write.xlsx(ttest, file = paste0(path3, "/III. Tables/ST_ttest.xlsx"))
write.xlsx(base_desc, file = paste0(path3, "/III. Tables/T2_base_desc.xlsx"))
write.xlsx(time_desc, file = paste0(path3, "/III. Tables/ST_time_desc.xlsx"))
write.xlsx(apr_desc, file = paste0(path3, "/III. Tables/ST_apr_desc.xlsx"))
write.xlsx(aug_desc, file = paste0(path3, "/III. Tables/ST_aug_desc.xlsx"))

# Save figures
ggsave(paste0(path3, '/IV. Figures/F1_dist.png'), wh_dist, width = 7, height = 3, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/Archive/ts_nolab.png'), ts_plot, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/SF_ts_stack.png'), ts_plot_stack, width = 8, height = 8, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/F2_corrplot_total.png'), corrplot_total, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/Archive/SF_corrplot_general.png'), corrplot_general, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/Archive/SF_corrplot_opth.png'), corrplot_opth, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/Archive/SF_corrplot_paed.png'), corrplot_paed, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/Archive/SF_corrplot_ent.png'), corrplot_ent, width = 10, height = 6, dpi = 300)
