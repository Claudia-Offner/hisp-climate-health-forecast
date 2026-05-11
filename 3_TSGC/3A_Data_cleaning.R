### ------------------------------------------------------------------------ ### 
###  TSGC - Data Cleaning (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/offne/Documents/GitHub/hisp-climate-health-forecast/')

#### DEPENDENCIES ####
source('3_TSGC/TSGC - Dependencies.R')
setwd(path_data)
file_export <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'


#### MAIN CODE ####


load('2_Active_Data/201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_SUB_clean.Rdata')

df <- df_heal_w
df$tp <- sample(0:50, nrow(df), replace = TRUE)
df$mx2t <- sample(-10:45, nrow(df), replace = TRUE)
# Create cumulative health outcomes
### NB: Set NA values to 1 in the cumulative variable - otherwise the growth curve model gets confused!!
df <- df %>%
  group_by(location) %>%
  arrange(location, date, .by_group = TRUE) %>%
  mutate(across(all_of(names(df)[6:33]), ~ cumsum(replace_na(., 1)), .names = '{.col}_c'))


# Load data
load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata')
df <- df_wide

# Remove Bilaspur
df <- df %>% filter(location != 'Bilaspur')

# Identify locations
# locations <- as.character(df %>% distinct(location, elev) %>% arrange(elev) %>% pull(location))
locations <- unique(df$location)

# Select year of interest
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-13')
df <- df[df$date >= startdate & df$date <= enddate,]

# Create 7-day rolling sum & log transform
df <- df %>%
  group_by(location) %>%
  mutate(tp_roll=rollapply(tp,   width = 7, FUN = sum,  align = "right", fill = NA)) %>%
  ungroup()
df$tp <- log(df$tp_roll+10)

# Create percentile look up
variables <- c('mx2t', 'tp')
cen_tables <- list()
datalist <- lapply(locations, function(region) df[df$location==region,])
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  # if (v=='tp'){
  #   cen <- exp(cen)-10
  # }
  colnames(cen) <- locations
  # cen <- cbind(cen, Average = rowMeans(cen))
  cen_tables[[v]] <- round(cen, 4)
}

# Set weather thresholds based on DLNM outputs
mx2t_20th <- as.data.frame(cen_tables$mx2t) %>%
  mutate(percentile = rownames(.)) %>%
  pivot_longer(cols = -percentile, names_to = "location", values_to = "thresh_mx2t") %>% 
  filter(percentile == "20%") %>% dplyr::select(location, thresh_mx2t)

tp_80th <- as.data.frame(cen_tables$tp) %>%
  mutate(percentile = rownames(.)) %>%
  pivot_longer(cols = -percentile, names_to = "location", values_to = "thresh_tp") %>% 
  filter(percentile == "80%") %>% dplyr::select(location, thresh_tp)

df <- df %>% 
  left_join(mx2t_20th , by = "location")  %>%
  left_join(tp_80th , by = "location")  %>%
  mutate(cold20_lag0 = ifelse(mx2t < thresh_mx2t, 1, 0),
         prec80_lag0 = ifelse(tp > thresh_tp, 1, 0)) %>%
  group_by(location) %>%
  mutate(
    cold20_lag3 = lag(cold20_lag0, 3),
    prec80_lag10 = lag(prec80_lag0, 10)
  ) %>%
  ungroup() %>%
  dplyr::select(-c(tp_roll, thresh_tp, thresh_mx2t), -matches("lag(2|[4-9]|1[1-4])"))

names(df)
table(df$cold20_lag0)
table(df$prec80_lag0)

#### EXPORT ####

save(df, file = paste0(path_git, '3_TSGC/', file_export, '.RData'))

