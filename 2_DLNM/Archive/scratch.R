### ------------------------------------------------------------------------ ### 
###  DNLM - Descriptive Analysis                                             ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Packages
library(dlnm); library(splines); library(dplyr); library(lubridate); library(xts);
library(ggplot2); library(ggthemes); library(patchwork); library(moments); 
library(purrr); library(tidyr); library(openxlsx)

# Turn off scientific notation
options(scipen=999)

# Path
user <- 'claer14'
file_import <- 'ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_wide'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 2 - Causal Inference/')
path4 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/2_DLNM/- Reporting/- figures/')
setwd(path1)

# Variables
weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
nm <- c(total_patients='Total patients', M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        general='General medicine', emergency='Emergency', ophthalmology='Ophthalmology', 
        pediatrics='Paediatric', orthopedics='Orthopedics', ent='ENT', 
        gynecology='Gynaecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        t2m='T2m (\u00B0C)', mx2t='Mx2t (\u00B0C)', mn2t='Mn2t (\u00B0C)', tp='Tp (mm)', fg10='Fg10 (km/h)', 
        closed='Closed', holiday='Holiday') 

### Function for plotting data with moving average
plot_moving_avg <- function(variable, MA=7, y.lab='', colour='red', plt.start.date=NULL, est.date.index, estimation.date.end){
  
  # Determine plot length
  if (is.null(plt.start.date)) {plt.start.date <- head(est.date.index, 1)}
  tmp.date <- min(estimation.date.end - 4, plt.start.date)
  s <- sprintf('%s/', format(tmp.date, '%Y-%m-%d'))
  variable <- variable[s]
  
  # Calculate the 7 day centered moving average
  ma.cent <- zoo::rollmean(variable, MA, align='center') 
  
  # Create df for plot (Data data vs. MA data)
  d <- setNames(cbind(variable, ma.cent), c('New Cases', 'Centered 7-day MA'))
  df.f <- data.frame(Date=index(d), 
                     New.Cases=coredata(d[,1]),
                     Centered.7.day.MA=coredata(d[,2]))
  
  
  # Plot Figure 1
  p1 <- ggplot(df.f, aes(x = Date)) +
    geom_line(aes(y = New.Cases, color = 'Real data'), linewidth = 0.1) +
    geom_line(aes(y = Centered.7.day.MA, color = 'Centered 7 day MA'), linewidth = 1) +
    scale_color_manual(name = '', values = c('Real data' = 'black', 'Centered 7 day MA' = colour)) +
    scale_y_continuous(n.breaks = 6) +
    scale_x_date(labels = scales::date_format('%d %b %y'), 
                 breaks = seq(min(df.f$Date), max(df.f$Date), by = '3 months')) +
    xlab('Day') + 
    ylab(y.lab) +
    theme_economist_white(gray_bg = FALSE, base_size = 12) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, vjust = 0.75),
      axis.title = element_text(size = rel(1), margin = margin(t = 10, r = 10)),
      axis.text = element_text(size = rel(1)),
      plot.title = element_text(face = 'bold', margin = margin(b = 5)),
      plot.caption = element_text(size = rel(1)),
      text = element_text(size = rel(1.1))
    ) +
    scale_linetype_manual(values = rep('solid', 3)) +
    scale_size_manual(values = rep(1, 3))
  
  return(p1) 
}

### Function summarize base statistics
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
  base_desc$Frequency[22:26] <- NA  # or "" if you prefer blank strings
  
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

### Function to summarize statistics by a grouping variable
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

### Function to replace variable names with set names
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


#### MAIN CODE ####

# Load data
load(paste0('2_DLNM/', file_import, '.RData'))

# Select year of interest
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-13')
df <- df[df$date >= startdate & df$date <= enddate,]

# Select variables of interest
vars <- c('closed', 'total_patients', 
          'F','M','under_5','age5_14','age15_24','age25_64','over_65',
          'general','ophthalmology','pediatrics','orthopedics','ent','gynecology',
          'dermatology','dentistry','surgery','emergency', 'psychiatry','oncology',
          't2m', 'mx2t', 'mn2t', 'tp', 'fg10')

# 1a. MT: Get baseline summary statistics (2015-2020) ####

base_desc <- get_desc(df, vars)

# 1b. ST: Get temporal summary statistics (2015-2020) ####

# Get summaries for each time grouping
summary_weekday <- summarise_by_group('weekday') %>% mutate(Time = as.character(Time), Time_Group = 'Weekday')
summary_month <- summarise_by_group('month') %>% mutate(Time = as.character(Time), Time_Group = 'Month')
summary_year <- summarise_by_group('year') %>% mutate(Time = as.character(Time), Time_Group = 'Year')

# Combine all summaries and pivot
time_desc <- bind_rows(summary_year, summary_month, summary_weekday) %>%
  dplyr::select(Time_Group, Time, variable, Summary) %>%
  pivot_wider(names_from = variable, values_from = Summary) %>%
  mutate(Time_Group = factor(Time_Group, levels = c("Year", "Month", "Weekday")))

# 2. MF: Get stacked time series visual (Predictors over time period: 2015-2020) ####

# Get moving average (no NAs)
clean <- df[!is.na(df$total_patients), ]
# Select year of interest
clean <- clean[clean$date >= startdate & clean$date <= enddate,]

# GET STACKED TS
variable <- xts(clean$total_patients,  order.by=clean$date)
ts_plot_tot <- plot_moving_avg(variable, MA=7, y.lab='Patients', colour='#8968CD',
                               plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

variable <- xts(clean$t2m,  order.by=clean$date)
ts_plot_t2m <- plot_moving_avg(variable, MA=7, y.lab='Celsius',colour='#D74B3C',
                               plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

variable <- xts(clean$tp,  order.by=clean$date)
ts_plot_tp <- plot_moving_avg(variable, MA=7, y.lab='mm', colour='#1D98EB',
                              plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

variable <- xts(clean$fg10,  order.by=clean$date)
ts_plot_fg10 <- plot_moving_avg(variable, MA=7, y.lab='km/hour', colour='#A6B727',
                                plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate) +
  theme(axis.title.x = element_blank())

# GET LEGEND ONLY
legend_data <- data.frame(
  x = 1:4, y = 1:4,
  group = c("Total Patients", "T2m", "Tp", "Fg10"))
legend_data$group <- factor(legend_data$group, levels = c("Total Patients", "T2m", "Tp", "Fg10"))
color_vals <- c("Total Patients" = "#8968CD",
                "T2m" = "#D74B3C",
                "Tp" = "#1D98EB",
                "Fg10" = "#A6B727")
manual_legend <- ggplot(legend_data, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = color_vals, name = NULL) +
  theme_void() +
  theme(legend.position = "bottom")
# Display only the legend
legend_only <- manual_legend + guides(color = guide_legend(nrow = 1)) + theme(legend.box.margin = margin(t = -20))


# Stack the plots with relative heights and insert legend at the bottom
ts_plot_stack <- (
  ts_plot_tot / ts_plot_t2m / ts_plot_tp / ts_plot_fg10 / legend_only)+ 
  plot_layout(heights = c(1, 1, 1, 1, 0.15))




#### EXPORT ####

# Replace names
base_desc <- replace_exact_matches(base_desc, nm)
names(time_desc) <- replace_exact_matches(names(time_desc), nm)

# Save tables
write.xlsx(base_desc, file = paste0(path3, "/III. Tables/MT_base_desc.xlsx"))
write.xlsx(time_desc, file = paste0(path3, "/III. Tables/ST_time_desc.xlsx"))

# Save figures
# ggsave(paste0(path3, '/IV. Figures/MF_ts_nolab.png'), ts_plot, width = 10, height = 6, dpi = 300)
ggsave(paste0(path3, '/IV. Figures/MF_ts_stack.png'), ts_plot_stack, width = 8, height = 8, dpi = 300)
ggsave(paste0(path4, 'MF_ts_stack.png'), ts_plot_stack, width = 9, height = 6, dpi = 300)
ggsave(paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/1_Data/HISP_Himachal_Pradesh/- Reporting/- figures/MF_ts_stack.png'), ts_plot_stack, width = 9, height = 6, dpi = 300)



#### Interaction Tests ####

# Test for interactions (p<0.05 indicates the interaction improves model)
test_interaction <- function(data, response, var1, var2, var3 = NULL, family = gaussian) {
  
  if (is.null(var3)) {
    # --- 2-way interaction ---
    f_main  <- as.formula(paste(response, "~", var1, "+", var2))
    f_inter <- as.formula(paste(response, "~", var1, "*", var2))
  } else {
    # --- 3-way interaction ---
    f_main  <- as.formula(paste(response, "~", var1, "+", var2, "+", var3, 
                                "+", var1, ":", var2, "+", var1, ":", var3, "+", var2, ":", var3))
    f_inter <- as.formula(paste(response, "~", var1, "*", var2, "*", var3))
  }
  
  # fit models
  m0 <- glm(f_main, data = data, family = family)
  m1 <- glm(f_inter, data = data, family = family)
  
  # likelihood ratio test
  result <- anova(m0, m1, test = "Chisq")
  return(result)
}

# Test for weather interactions
test_interaction(df, "total_patients", "tp", "fg10", family = gaussian) # p=0.01472 *
test_interaction(df, "total_patients", "t2m", "tp", family = gaussian) # p=0.15
test_interaction(df, "total_patients", "t2m", "fg10", family = gaussian) # p=0.0005258 ***
test_interaction(df, "total_patients", "mx2t", "tp", family = gaussian) # p=0.18
test_interaction(df, "total_patients", "mx2t", "fg10", family = gaussian) # p=0.00252 **
test_interaction(df, "total_patients", "mn2t", "tp", family = gaussian) # p=0.0756 .
test_interaction(df, "total_patients", "mn2t", "fg10", family = gaussian) # p=0.06709 .
# Check 3-way interaction (NULL)
test_interaction(df, "total_patients", "t2m", "tp", "fg10", family = gaussian) # p= 0.4012
# = Potential interactions: tp*fg10; t2m*fg10; mx2t*fg10


# Test for socio-dem interactions
test_interaction(df, "total_patients", 't2m', "under_5", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "age5_14", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "age15_24", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "age25_64", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "over_65", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "F", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "M", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "under_5", "M", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "under_5", "F", family = gaussian) # p=0.0000000000007516 ***
test_interaction(df, "total_patients", 't2m', "age5_14", "M", family = gaussian) # p=0.08904 .
test_interaction(df, "total_patients", 't2m', "age5_14", "F", family = gaussian) # p=0.00000000000643 ***
test_interaction(df, "total_patients", 't2m', "age15_24", "M", family = gaussian) # p=0.00000000000000022 ***
test_interaction(df, "total_patients", 't2m', "age15_24", "F", family = gaussian) # p=0.00000000009531 ***
test_interaction(df, "total_patients", 't2m', "age25_64", "M", family = gaussian) # p=0.00000000000003963 ***
test_interaction(df, "total_patients", 't2m', "age25_64", "F", family = gaussian) # p=0.00956 **
test_interaction(df, "total_patients", 't2m', "over_65", "M", family = gaussian) # p=0.00003949 ***
test_interaction(df, "total_patients", 't2m', "over_65", "F", family = gaussian) # p=0.0000000000000005428 ***

# Maybe interaction between age and gender across weather variables



#########################


# # Main exposure-response, by city
# (city_tmp_er <- get_city_legend(pred='mx2t', plot=main_res$p_er_tmp_city + plot_annotation(title='Exposure-response, by city: Temperature')))
# (city_tp_er <- get_city_legend(pred='tp', plot=main_res$p_er_tp_city + plot_annotation(title='Exposure-response, by city: Precipitation')))

# Main lag-response, by city
# (city_tmp_lr <- main_res$p_lr_tmp_city + plot_annotation(title=replace_exact_matches('wb', nm)))
# (city_tp_lr <- main_res$p_lr_tp_city + plot_annotation(title=replace_exact_matches('tp', nm)))

# City case study: Chamba & Dharamshala
dharam_res <- format_results(output, outcome, city='Dharamshala')


dharam_lr <- get_city_legend(wrap_plots(dharam_res$p_lr_tmp_city, dharam_res$p_lr_tp_city,
                                        nrow=2)+plot_annotation(title='Lag-response: Dharamshala', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
dharam_er <- get_city_legend(wrap_plots(dharam_res$p_er_tmp_city, dharam_res$p_er_tp_city,
                                        nrow=2)+plot_annotation(title='Exposure-response: Dharamshala', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))

cham_res <- format_results(output, outcome, city='Chamba')
cham_lr <- get_city_legend(wrap_plots(cham_res$p_lr_tmp_city, cham_res$p_lr_tp_city, 
                                      nrow=2)+plot_annotation(title='Lag-response: Chamba', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
cham_er <- get_city_legend(wrap_plots(cham_res$p_er_tmp_city, cham_res$p_er_tp_city, 
                                      nrow=2)+plot_annotation(title='Exposure-response: Chamba', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))


# Main percentile look up
cen_tables <- list()
datalist <- lapply(locations, function(region) df[df$location==region,])

for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)-10
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_tables[[v]] <- round(cen, 4)
}
names(cen_tables) <- c('Maximum temperature', 'Precipitation') # Wind speed

#### EXPORT #####

# FIGURES
exports <- list('1A_total_er_overall'=list(p=total_er_overall, w=9, h=5),
                '1B_total_lr'=list(p=total_lr, w=8, h=6),
                '1C_total_er'=list(p=total_er, w=8, h=6),
                '1D_con_tmp'=list(p=con_tmp, w=5, h=4),
                '1D_con_tp'=list(p=con_tp, w=5, h=4),
                '1D_con_tmp_sig'=list(p=con_tmp_sig, w=5, h=4),
                '1D_con_tp_sig'=list(p=con_tp_sig, w=5, h=4),
                '1D_city_tmp_er'=list(p=city_tmp_er, w=9, h=5),
                '1E_city_tp_er'=list(p=city_tp_er, w=9, h=5),
                # '1F_city_fg10_er'=list(p=city_fg10_er, w=9, h=5),
                '1G_cham_lr'=list(p=cham_lr, w=8, h=6),
                '1H_cham_er'=list(p=cham_er, w=8, h=6),
                '1I_dharam_lr'=list(p=dharam_lr, w=8, h=6),
                '1J_dharam_er'=list(p=dharam_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300, bg='white')
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"II. dlnm/tab_main_perc.xlsx"), overwrite = TRUE)

# Analysis 2 - City-specific models (city params) ####

# Get city-specific results
for (l in locations){
  sub <- df %>% filter(location==l)
  res <- run_bestfit(df=sub, locations, outcome, params=final_params[[l]], pooled=FALSE, city=l)
  results[[l]] <- res
  print(l)
}


(spec_dharam_er <- results[['Dharamshala']]$results$p_er + plot_annotation(title=paste0('Exposure-response: Dharamshala')))
(spec_cham_er <- results[['Chamba']]$results$p_er + plot_annotation(title=paste0('Exposure-response: Chamba')))

#### EXPORT #####

# FIGURES
exports <- list('2A_spec_dharam_er'=list(p=spec_dharam_er, w=8, h=6),
                '2B_spec_cham_er'=list(p=spec_cham_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300, bg='white')
}


# Analysis 3 - Seasonal differences (pooled params) ####

# Winter/pre-monsoon model (months=5)
df_win <- df %>% filter(season %in% c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)'))
win <- run_bestfit(df_win, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$winter <- win
win_res <- win$results

### Winter - Overall
win_res$p_er_overall <- wrap_plots(lapply(win_res$p_er_overall, function(p) p + theme(legend.position = "none")), nrow=1)
win_er_overall <- wrap_plots(get_title('Winter'), win_res$p_er_overall, nrow=2, heights=c(1, 40))
win_lr <- wrap_plots(get_title('Lag-response: Winter'), win_res$p_lr, nrow=2, heights=c(1, 40))
win_er <- wrap_plots(get_title('Exposure-response: Winter'), win_res$p_er, nrow=2, heights=c(1, 40))

### Winter - Percentile look up
cen_win_tables <- list()
datalist <- lapply(locations, function(region) df_win[df_win$location==region,])
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_win_tables[[v]] <- round(cen, 4)
}
names(cen_tables) <- c('Maximum temperature', 'Precipitation') #, 'Wind speed'


# Monsoon/Post-monsoon model (months=7)
df_sum <- df %>% filter(season %in% c('Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)'))
sum <- run_bestfit(df_sum, locations, outcome, params=final_params$pooled, pooled=TRUE, city=NULL)
results$summer <- sum
sum_res <- sum$results

### Summer - Overall
sum_er_overall <- wrap_plots(get_title('Summer'), sum_res$p_er_overall, nrow=2, heights=c(1, 40))
sum_lr<- wrap_plots(get_title('Lag-response: Summer'), sum_res$p_lr, nrow=2, heights=c(1, 40))
sum_er<- wrap_plots(get_title('Exposure-response: Summer'), sum_res$p_er, nrow=2, heights=c(1, 40))

### Summer - Percentile look up
cen_sum_tables <- list()
datalist <- lapply(locations, function(region) df_sum[df_sum$location==region,])
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)
  }
  colnames(cen) <- locations
  cen <- cbind(cen, Average = rowMeans(cen))
  cen_sum_tables[[v]] <- round(cen, 4)
}
names(cen_tables) <- c('Maximum temperature', 'Precipitation') #, 'Wind speed'


# MERGE PLOTS FOR COMPARISON
(seas_overall <- wrap_plots(win_er_overall, sum_er_overall, nrow=2))
(seas_lr <- wrap_plots(win_lr, plot_spacer(), sum_lr, nrow=1, widths = c(1, 0.08, 1)))
(seas_er <- wrap_plots(win_er, plot_spacer(), sum_er, nrow=1, widths = c(1, 0.08, 1)))


#### EXPORT #####
exports <- list('3A_seas_overall'=list(p=seas_overall, w=8, h=6),
                '3B_seas_lr'=list(p=seas_lr, w=15, h=6),
                '3C_seas_er'=list(p=seas_er, w=15, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300)
}

# TABLES
wb <- createWorkbook()
for (name in names(cen_win_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_win_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"II. dlnm/tab_winter_perc.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
for (name in names(cen_sum_tables)) {
  addWorksheet(wb, sheetName = name)               
  writeData(wb, sheet = name, cen_sum_tables[[name]], rowNames=TRUE)
}
saveWorkbook(wb, file=paste0(path_report,"II. dlnm/tab_summer_perc.xlsx"), overwrite = TRUE)


# Analysis 4 - Population differences (pooled params) ####

# # Get total tables and id relevant lags/exposure
# x <- meta_checks(output$meta_tmp, 1.001)
# meta_checks(output$meta_tp, 1.001)
# meta_checks(output$meta_pm25, 1.000)

# Multiple outcomes
outcomes <- h_cols[1:8]
er0_tmp <- list()
er3_tmp <- list()
er10_tp <- list()
for (i in outcomes){
  print(i)
  res <- run_bestfit(df, locations, i, params=final_params$pooled, pooled=TRUE, city=NULL)
  results[[i]] <- res
  pop_res <- res$results
  er0_tmp[[i]] <- pop_res$lags$`0`$p_er_tmp$p_blup + labs(title=replace_exact_matches(i, nm))
  er3_tmp[[i]] <- pop_res$lags$`3`$p_er_tmp$p_blup + labs(title=replace_exact_matches(i, nm))
  er10_tp[[i]] <- pop_res$lags$`10`$p_er_tp$p_blup + labs(title=replace_exact_matches(i, nm))
}

(er0_tmp <- (wrap_plots(er0_tmp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Exposure-response (Lag=0)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er3_tmp <- (wrap_plots(er3_tmp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Exposure-response (Lag=3)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))
(er10_tp <- (wrap_plots(er10_tp, ncol=4)& theme(legend.position="none"))+plot_annotation(title='Exposure-response (Lag=10)', theme=theme(plot.title=element_text(hjust=0.5, size=16, face="bold"))))

all_outcomes <- list(er0_tmp=er0_tmp, er3_tmp=er3_tmp, er10_tp=er10_tp) #, er0_fg10=er0_fg10

#### EXPORT ####

# Population plots
names(all_outcomes) <- paste0("disag_", names(all_outcomes))
for (n in names(all_outcomes)) {
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_4_', n, '.png'), plot=all_outcomes[[n]], width=8, height=5, dpi=300)
}

# Analysis 5 - Hamirpur emergency ward (emergency params) ####

# Set outcome and site
outcome <- 'emergency'
l <- 'Hamirpur'

# Main results
sub <- df %>% filter(location==l)
res <- run_bestfit(sub, locations, outcome, params=final_params$emergency, pooled=FALSE, city=l)
results$emergency <- res
emg_res <- res$results

# Overall
emg_res$p_er_overall <- wrap_plots(lapply(emg_res$p_er_overall, function(p) p + theme(legend.position = "none")), nrow=1)
(emg_er_overall <- emg_res$p_er_overall + plot_annotation(title='Overall: Hamirpur emergency ward', theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold'))))
(emg_lr <- emg_res$p_lr + plot_annotation(title='Lag-response: Hamirpur emergency ward'))
(emg_er <- emg_res$p_er + plot_annotation(title='Exposure-response: Hamirpur emergency ward'))


#### EXPORT ####

exports <- list('5A_emg_overall'=list(p=emg_er_overall, w=9, h=5),
                '5B_emg_lr'=list(p=emg_lr, w=8, h=6),
                '5C_emg_er'=list(p=emg_er, w=8, h=6))

for (e in names(exports)){
  ggsave(filename=paste0(path_report, 'II. dlnm/fig_', e, '.png'), plot=exports[[e]]$p, width=exports[[e]]$w, height=exports[[e]]$h, dpi=300)
}

#################################

# # META REGRESSION: Elevation and population density
# meta_reg <- function(meta, xlab='Exposure (%)', predictor='elev', perc=c(0.1, 0.9), 
#                      col='darkred'){
#   
#   # Get unique predictors
#   pred <- sapply(locations, function(loc) {
#     if(predictor=='popu'){
#       mean(unique(df[[predictor]][df$location==loc]))
#     } else{
#       unique(df[[predictor]][df$location==loc])
#     }
#   })
#   
#   # Remove duplicates before predictions
#   rows_with_dup <- apply(meta$cen, 1, function(row) any(duplicated(row)))
#   cen <- meta$cen[!rows_with_dup, ]
#   
#   # Get pooled regression estimates
#   (mv_reg <- mvmeta(meta$ymat~pred,meta$Slist,method="reml"))
#   pred_perc <- round(quantile(pred,perc,1))
#   predict <- predict(mv_reg,data.frame(pred=pred_perc),vcov=T)
#   cp_p1 <- crosspred(meta$b_mean,coef=predict[[1]]$fit, vcov=predict[[1]]$vcov,model.link="log", 
#                      at=rowMeans(cen), cen=rowMeans(cen)["50%"])
#   cp_p2 <- crosspred(meta$b_mean,coef=predict[[2]]$fit, vcov=predict[[2]]$vcov,model.link="log", 
#                      at=rowMeans(cen), cen=rowMeans(cen)["50%"])
#   
#   # Format tables for plotting
#   dfp1 <- data.frame(exp=as.numeric(sub("%", "", names(rowMeans(cen)))),
#                      perc=paste0(perc[1]*100, "%"), RR =cp_p1$allRRfit, 
#                      low=cp_p1$allRRlow, high=cp_p1$allRRhigh)
#   dfp2 <- data.frame(exp=as.numeric(sub("%", "", names(rowMeans(cen)))),
#                      perc=paste0(perc[2]*100, "%"), RR =cp_p2$allRRfit, 
#                      low=cp_p2$allRRlow, high=cp_p2$allRRhigh)
#   df_plot <- rbind(dfp1, dfp2)
#   df_plot$cross_one <- ifelse(df_plot$low <= 1 & df_plot$high >= 1, "cross", "no_cross")
#   
#   # Plot
#   ggplot(df_plot, aes(x = exp, y = RR, color = perc)) + #linetype = perc
#     geom_hline(yintercept = 1, color = 'gray20') +
#     geom_line(size = 1.2) +
#     geom_ribbon(aes(ymin = low, ymax = high, fill = perc), alpha = 0.2, color = NA) +
#     labs(
#       x = xlab,
#       y = "RR",
#       title = "",
#       color = "",
#       fill  = "",
#       # linetype = ""
#     ) +
#     scale_color_manual(values = c('gray20', col)) +
#     scale_fill_manual(values = c('gray20', col)) +
#     # scale_linetype_manual(values = c("dotted", "solid")) +  # assign specific line types
#     theme_classic() +
#     theme(
#       axis.line = element_line(color = 'gray20'),
#       axis.text = element_text(color = 'gray20'),
#       axis.ticks = element_line(color = 'gray20'),
#       plot.title = element_text(face = 'bold', hjust = 0.5, size = 12),
#       legend.position = 'bottom'
#     )
#   
# }
# 
# 
# meta_list <- list('meta_wb'=list(xlab='Mean temp. (%)', col='#D74B3C'),
#                   'meta_tp'=list(xlab='Rainfall (%)', col='#1D98EB'),
#                   'meta_rh'=list(xlab='Relative humidity (%)', col='#A6B727'))
#                   # 'meta_fg10'=list(xlab='Wind speed (%)', col='#8968CD'),
#                   # 'meta_pm25'=list(xlab='Air pollution (%)', col='#E569BF'))
# elev_plots <- list()
# popu_plots <- list()
# for(v in names(meta_list)){
#   elev_reg <- meta_reg(output[[v]], predictor='elev', perc=c(0.51, 0.49),
#                        xlab=meta_list[[v]]$xlab, meta_list[[v]]$col)
#   elev_plots[[v]] <- elev_reg
#   popu_reg <- meta_reg(output[[v]], predictor='popu', perc=c(0.51, 0.49),
#                        xlab=meta_list[[v]]$xlab, meta_list[[v]]$col)
#   popu_plots[[v]] <- popu_reg
# }
# 
# (elev_res <- wrap_plots(elev_plots, ncol=3)+plot_annotation(title='Exposure-resonse by elevation percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# (popu_res <- wrap_plots(popu_plots, ncol=3)+plot_annotation(title='Exposure-resonse by population density percentile', theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))))
# 







