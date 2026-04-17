### ------------------------------------------------------------------------ ### 
###  Plot the results of the test period
### ------------------------------------------------------------------------ ### 

library(ggplot2); library(RColorBrewer); library(gridExtra); library(cowplot)

Sys.setlocale("LC_ALL","English")

# Prepare data  ------------------------------------------------

date_max<-max(dates_analysis)
dat2<-data_rt %>%
  filter(date_analysis==max(data_rt$date_analysis)) %>%
  filter(date>=min(dates_analysis)-remove_last_n & date<=date_max) %>%
  dplyr::select(date,location,outcome) %>% #,iICU_smooth,inHosp_smooth,inICU_smooth
  pivot_longer(cols=c("outcome"), #"iICU_smooth","inHosp_smooth","inICU_smooth"
               names_to="var",values_to="value") %>%
  mutate(var=gsub("_smooth","",var)) %>%
  filter(!is.element(location,c("COR","national")))



# Load projections from components and ensemble models -----------------------

proj<-readRDS(file=paste0(path_results,"Ensemble_ONLY.rds"))
proj$model<-factor(proj$model,levels=c(models,"Ensemble"), ordered=T)


# Compute evaluation metrics  ------------------------------------

# Error
res<-proj %>%
  left_join(dat2,by=c("date_proj"="date","location"="location","var"="var")) %>%
  mutate(in_interval95 = value >= lower_5 & value <= upper_5) %>%
  mutate(abs_error = point - value) %>%
  mutate(bias = (point - value) / pmax(value, 1))

# WIS
# from https://github.com/reichlab/covid19-forecast-evals/blob/main/code/process-zoltar-scores.R

intervals <- list(
  list(suffix="2",   range=98), list(suffix="5",   range=95),
  list(suffix="10",  range=90), list(suffix="20",  range=80),
  list(suffix="30",  range=70), list(suffix="40",  range=60),
  list(suffix="50",  range=50), list(suffix="60",  range=40),
  list(suffix="70",  range=30), list(suffix="80",  range=20),
  list(suffix="90",  range=10)
)

for (iv in intervals) {
  res[[paste0("IS_", iv$suffix)]] <- scoringutils:::interval_score(
    observed       = res$value,
    lower          = res[[paste0("lower_", iv$suffix)]],
    upper          = res[[paste0("upper_", iv$suffix)]],
    interval_range = iv$range,
    weigh          = TRUE,
    separate_results = FALSE
  )
}

res <- res %>%
  mutate(IS_100 = abs(point-value)) %>%
  mutate(wis = (IS_2+IS_5+IS_10+IS_20+IS_30+IS_40+IS_50+
                  IS_60+IS_70+IS_80+IS_90+0.5*IS_100)/11.5)  


# Plot settings  ------------------------------------------------

pal<-c("black",brewer.pal(10,"Paired"),"turquoise","darkred")
val_size<-c(rep(0.5,length(models)),1)
val_linetype<-c("dashed", rep("solid",length(models)-1),"dashed")
lab_mod<-levels(droplevels(proj$model))



# Figure 3 : Forecast trajectories of hospital admissions ------

# National level, all models
myvar<-"outcome"
y_min <- proj %>% filter(var == myvar) %>% filter(location == "metropolitan") %>%
  filter(model != "baseline") %>% filter(date_proj <= date_max) %>%
  pull(point) %>% min(na.rm = TRUE)
y_max <- proj %>% filter(var == myvar) %>% filter(location == "metropolitan") %>%
  filter(model != "baseline") %>% filter(date_proj <= date_max) %>%
  pull(point) %>% max(na.rm = TRUE)

p1 <- proj %>%
  filter(var == myvar) %>%
  filter(location == "metropolitan") %>%
  filter(model != "baseline") %>%
  filter(date_proj <= date_max) %>%
  ggplot() +
  geom_line(aes(x = date_proj, y = point, color = as.factor(prediction_horizon), group = date)) +
  geom_line(data = dat2[dat2$location == "metropolitan" & dat2$var == myvar & dat2$date >= min(dates_analysis),],
            aes(x = date, y = value), size = 0.8) +
  facet_wrap(~ model, scale = "free", ncol = 4) +
  xlab("") +
  scale_y_continuous("Hospital Admissions", limits = c(y_min, y_max)) +
  scale_color_discrete("Prediction\nhorizon") +
  coord_cartesian(xlim = c(min(dates_analysis), max(dates_analysis))) +
  theme_bw() +
  labs(tag = "A")

x <- proj %>%
  filter(var == myvar) %>%
  filter(location == "metropolitan") %>%
  filter(model != "baseline") %>%
  filter(date_proj <= date_max)

# locational level, ensemble model
lab_p <- c("3", "7", "14")
val_p<-c("powderblue","navajowhite",brewer.pal(10,"Paired")[5])
val_p2<-c("turquoise4","burlywood4",brewer.pal(10,"Paired")[6])

myvar<-"outcome"
p2<-proj %>%
  filter(var==myvar) %>%
  filter(prediction_horizon %in% c(3,7,14)) %>%
  mutate(prediction_horizon2 = factor(prediction_horizon,levels=c(14,7,3),ordered = T)) %>% 
  filter(model=="Ensemble") %>% 
  filter(location != "metropolitan") %>%
  filter(date_proj<=date_max) %>% 
  ggplot() +
  geom_ribbon(aes(x= date_proj, ymin= lower_5, ymax= upper_5, fill=as.factor(prediction_horizon2)),alpha=0.7)+
  geom_line(aes(x = date_proj, y = point, color = as.factor(prediction_horizon2))) +
  geom_line(data = dat2[dat2$location != "metropolitan" & dat2$var==myvar & dat2$date>=min(dates_analysis),], 
            aes(x = date, y = value, size = var)) +
  facet_wrap(~location, scale = "free") +
  xlab("") +
  scale_y_continuous("Hospital Admissions", limits = c(0, NA)) +
  scale_color_manual("Prediction\nhorizon", labels = rev(lab_p), values = val_p2) +
  scale_fill_manual("Prediction\nhorizon", labels = rev(lab_p), values = val_p) +
  scale_size_manual("",labels="Data",values=0.8) + 
  coord_cartesian(xlim=c(min(dates_analysis),max(dates_analysis))) +
  theme_bw()+
  labs(tag="B")


ggsave(filename = paste0(path_figures, "Figure3.png"),
       plot = plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1)),
       width = 8, height = 10, units = "in", dpi = 300)


# Figure S9 : RMSE and WIS  ----------------------------------------

# RMSE national
p1<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(location == "metropolitan") %>% 
  filter(var=="outcome") %>% 
  group_by(model,var,prediction_horizon) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, color = model , 
                size=model,linetype=model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous("RMSE") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  scale_color_manual("Model",labels=lab_mod,values=pal) +
  scale_size_manual("Model",labels=lab_mod,values=val_size) +
  scale_linetype_manual("Model",labels=lab_mod,values=val_linetype) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(4, "mm"))+
  guides(color=guide_legend(ncol=1)) +
  labs(tag = "A") + 
  ggtitle("RMSE, Himachal Pradesh")


# WIS national
p2<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(location == "metropolitan") %>% 
  filter(var=="outcome") %>% 
  group_by(model,var,prediction_horizon) %>%
  dplyr::summarize(average_bias = mean(wis,na.rm=T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, color = model , 
                size=model,linetype=model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous("WIS") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  scale_color_manual("Model",labels=lab_mod,values=pal) +
  scale_size_manual("Model",labels=lab_mod,values=val_size) +
  scale_linetype_manual("Model",labels=lab_mod,values=val_linetype) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(4, "mm"))+
  guides(color=guide_legend(ncol=1)) +
  labs(tag = "C") + 
  ggtitle("WIS, Himachal Pradesh")


# RMSE location
p1_reg<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(location %in% locations_hosp) %>% 
  filter(var=="outcome") %>% 
  group_by(model,var,prediction_horizon) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, color = model , 
                size=model,linetype=model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous("RMSE") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  scale_color_manual("Model",labels=lab_mod,values=pal) +
  scale_size_manual("Model",labels=lab_mod,values=val_size) +
  scale_linetype_manual("Model",labels=lab_mod,values=val_linetype) +
  # coord_cartesian(ylim=c(0,100)) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(4, "mm"))+
  guides(color=guide_legend(ncol=1))+
  labs(tag = "B") + 
  ggtitle("RMSE, hospital level")


# WIS location
p2_reg<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(location %in% locations_hosp) %>% 
  filter(var=="outcome") %>% 
  group_by(model,var,prediction_horizon) %>%
  dplyr::summarize(average_bias = mean(wis,na.rm=T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, color = model , 
                size=model,linetype=model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous("WIS") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  scale_color_manual("Model",labels=lab_mod,values=pal) +
  scale_size_manual("Model",labels=lab_mod,values=val_size) +
  scale_linetype_manual("Model",labels=lab_mod,values=val_linetype) +
  # coord_cartesian(ylim=c(0,50)) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(4, "mm"))+
  guides(color=guide_legend(ncol=1))+
  labs(tag = "D") + 
  ggtitle("WIS, hospital level")

ggsave(filename = paste0(path_figures, "FigureS9.png"),
       plot = plot_grid(p1, p1_reg, p2, p2_reg, ncol = 2, align = "v", axis = "lr"),
       width = 8, height = 6, units = "in", dpi = 200)

# Figure 4 : Ranks and MAPE ---------------------------------------

# Ranks
rnk<-res %>%
  filter(var=="outcome") %>%
  filter(location %in% locations) %>%
  filter(date_proj<=date_max) %>%
  group_by(location,model) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  arrange(location,average_bias)
rnk$rank<-ave(1:nrow(rnk),rnk$location,FUN=seq_along)
rnk_mean<-res %>%
  filter(var=="outcome") %>%
  filter(location %in% locations_hosp) %>%
  filter(date_proj<=date_max) %>%
  group_by(model) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  mutate(location="Reg. Ave.") %>%
  arrange(location,average_bias)
rnk_mean$rank<-1:nrow(rnk_mean)
rnk<-bind_rows(rnk,rnk_mean)

x_levels<-c(locations,"Reg. Ave.")
rnk$location<-factor(rnk$location,levels=x_levels, ordered=T)

# order models by rank
rnk$model[rnk$location=="Reg. Ave."]
rnk$model2<-factor(rnk$model,levels=rev(c("baseline","GAM1","ARIMA1", 
                                          "GAM2", "ARIMA2", 'TSGC',"Ensemble")),ordered=T)      

p1<-ggplot(rnk, 
           aes(y=model2, fill=rank, x=location)) + 
  geom_tile() + ylab(NULL) + xlab(NULL) +
  geom_text(aes(label=rank), size=5) +
  scale_fill_gradient(low = "dodgerblue4", high="white", limits=c(1,length(models)+1)) +
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 0.5, vjust = .5),
        axis.text.y = element_text(size=12))+
  labs(tag="A")


# MAPE for all targets, national
p2<-res %>%
  filter(location=="metropolitan") %>%  
  filter(model=="Ensemble") %>%
  filter(date_proj<=date_max) %>%
  group_by(var,prediction_horizon) %>%
  dplyr::summarize(average_bias = mean(abs(bias),na.rm=T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, group=var, linetype=var),col="darkred",size=1) +
  scale_linetype_manual("Target",labels=c("Hosp. adm.","ICU adm.","GW beds","ICU beds"),values=c("solid","dotdash","dashed","dotted")) +
  scale_y_continuous("MAPE") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  # coord_cartesian(ylim=c(0,0.35)) +
  theme_bw()+
  labs(tag = "B") + 
  ggtitle("MAPE, Himachal Pradesh")

# MAPE for all targets, locations
p3<-res %>%
  filter(location %in% locations_hosp) %>%
  filter(model=="Ensemble" ) %>%
  filter(date_proj<=date_max) %>% 
  group_by(var,prediction_horizon) %>%
  dplyr::summarize(average_bias = mean(abs(bias),na.rm=T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, group=var, linetype=var),col="darkred",size=1) +
  scale_linetype_manual("Target",labels=c("Hosp. adm.","ICU adm.","GW beds","ICU beds"),values=c("solid","dotdash","dashed","dotted")) +
  scale_y_continuous("MAPE") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  # coord_cartesian(ylim=c(0,0.35)) +
  theme_bw()+
  labs(tag = "C") + 
  ggtitle("MAPE, hospital level")


bottom_row <- plot_grid(p2, p3, ncol=2)


ggsave(filename    = paste0(path_figures, "Figure4.png"),
  plot = plot_grid(p1, bottom_row, ncol = 1, rel_heights = c(1.3, 1)),
  width = 8, height = 8, units = "in", dpi = 300)

