### ------------------------------------------------------------------------ ### 
###  Plot the results of the test period
### ------------------------------------------------------------------------ ### 

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

Sys.setlocale("LC_ALL","English")

# Prepare data  ------------------------------------------------

date_max<-max(dates_analysis)
dat2<-data_rt %>%
  filter(date_analysis==max(data_rt$date_analysis)) %>%
  filter(date>=min(dates_analysis)-remove_last_n & date<=date_max) %>%
  dplyr::select(date,region,iHosp_smooth) %>% #,iICU_smooth,inHosp_smooth,inICU_smooth
  pivot_longer(cols=c("iHosp_smooth"), #"iICU_smooth","inHosp_smooth","inICU_smooth"
               names_to="var",values_to="smooth_value") %>%
  mutate(var=gsub("_smooth","",var)) %>%
  filter(!is.element(region,c("COR","national")))



# Load projections from components and ensemble models -----------------------

proj<-readRDS(file=paste0(path_results,"Ensemble_ONLY.rds"))
proj$model<-factor(proj$model,levels=c(models,"Ensemble"), ordered=T)


# Compute evaluation metrics  ------------------------------------

# Error
res<-proj %>%
  left_join(dat2,by=c("date_proj"="date","region"="region","var"="var")) %>%
  mutate(in_interval95 = smooth_value >= lower_5 & smooth_value <= upper_5) %>%
  mutate(abs_error = point - smooth_value) %>%
  mutate(bias = (point - smooth_value) / pmax(smooth_value, 1))

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
    observed       = res$smooth_value,
    lower          = res[[paste0("lower_", iv$suffix)]],
    upper          = res[[paste0("upper_", iv$suffix)]],
    interval_range = iv$range,
    weigh          = TRUE,
    separate_results = FALSE
  )
}

res <- res %>%
  mutate(IS_100 = abs(point-smooth_value)) %>%
  mutate(wis = (IS_2+IS_5+IS_10+IS_20+IS_30+IS_40+IS_50+
                  IS_60+IS_70+IS_80+IS_90+0.5*IS_100)/11.5)  


# Plot settings  ------------------------------------------------

pal<-c("black",brewer.pal(10,"Paired"),"turquoise","darkred")
val_size<-c(rep(0.5,length(models)),1)
val_linetype<-c("dashed", rep("solid",length(models)-1),"dashed")
lab_mod<-levels(droplevels(proj$model))



# Figure 3 : Forecast trajectories of hospital admissions ------

# National level, all models
myvar<-"iHosp"
y_min <- proj %>% filter(var == myvar) %>% filter(region == "metropolitan") %>%
  filter(model != "baseline") %>% filter(date_proj <= date_max) %>%
  pull(point) %>% min(na.rm = TRUE)
y_max <- proj %>% filter(var == myvar) %>% filter(region == "metropolitan") %>%
  filter(model != "baseline") %>% filter(date_proj <= date_max) %>%
  pull(point) %>% max(na.rm = TRUE)

p1 <- proj %>%
  filter(var == myvar) %>%
  filter(region == "metropolitan") %>%
  filter(model != "baseline") %>%
  filter(date_proj <= date_max) %>%
  ggplot() +
  geom_line(aes(x = date_proj, y = point, color = as.factor(prediction_horizon), group = date)) +
  geom_line(data = dat2[dat2$region == "metropolitan" & dat2$var == myvar & dat2$date >= min(dates_analysis),],
            aes(x = date, y = smooth_value), size = 0.8) +
  facet_wrap(~ model, scale = "free", ncol = 4) +
  xlab("") +
  scale_y_continuous("Hospital Admissions", limits = c(y_min, y_max)) +
  scale_color_discrete("Prediction\nhorizon") +
  coord_cartesian(xlim = c(min(dates_analysis), max(dates_analysis))) +
  theme_bw() +
  labs(tag = "A")

x <- proj %>%
  filter(var == myvar) %>%
  filter(region == "metropolitan") %>%
  filter(model != "baseline") %>%
  filter(date_proj <= date_max)

# Regional level, ensemble model
lab_p <- c("3", "7", "14")
val_p<-c("powderblue","navajowhite",brewer.pal(10,"Paired")[5])
val_p2<-c("turquoise4","burlywood4",brewer.pal(10,"Paired")[6])

myvar<-"iHosp"
p2<-proj %>%
  filter(var==myvar) %>%
  filter(prediction_horizon %in% c(3,7,14)) %>%
  mutate(prediction_horizon2 = factor(prediction_horizon,levels=c(14,7,3),ordered = T)) %>% 
  filter(model=="Ensemble") %>% 
  filter(region != "metropolitan") %>%
  filter(date_proj<=date_max) %>% 
  ggplot() +
  geom_ribbon(aes(x= date_proj, ymin= lower_5, ymax= upper_5, fill=as.factor(prediction_horizon2)),alpha=0.7)+
  geom_line(aes(x = date_proj, y = point, color = as.factor(prediction_horizon2))) +
  geom_line(data = dat2[dat2$region != "metropolitan" & dat2$var==myvar & dat2$date>=min(dates_analysis),], 
            aes(x = date, y = smooth_value, size = var)) +
  facet_wrap(~region, scale = "free") +
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
       width = 10, height = 12, units = "in", dpi = 300)


# Figure S9 : RMSE and WIS  ----------------------------------------

# RMSE national
p1<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(region == "metropolitan") %>% 
  filter(var=="iHosp") %>% 
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
  ggtitle("RMSE, metropolitan France")


# WIS national
p2<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(region == "metropolitan") %>% 
  filter(var=="iHosp") %>% 
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
  ggtitle("WIS, metropolitan France")


# RMSE regional
p1_reg<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(region %in% regions_metro) %>% 
  filter(var=="iHosp") %>% 
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
  coord_cartesian(ylim=c(0,100)) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(4, "mm"))+
  guides(color=guide_legend(ncol=1))+
  labs(tag = "B") + 
  ggtitle("RMSE, regional level")


# WIS regional
p2_reg<-res %>%
  filter(date_proj<=date_max) %>% 
  filter(region %in% regions_metro) %>% 
  filter(var=="iHosp") %>% 
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
  coord_cartesian(ylim=c(0,50)) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(4, "mm"))+
  guides(color=guide_legend(ncol=1))+
  labs(tag = "D") + 
  ggtitle("WIS, regional level")

ggsave(filename = paste0(path_figures, "FigureS9.png"),
       plot = plot_grid(p1, p1_reg, p2, p2_reg, ncol = 2, align = "v", axis = "lr"),
       width = 9, height = 7, units = "in", dpi = 200)

# Figure 4 : Ranks and MAPE ---------------------------------------

# Ranks
rnk<-res %>%
  filter(var=="iHosp") %>%
  filter(region %in% regions) %>%
  filter(date_proj<=date_max) %>%
  group_by(region,model) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  arrange(region,average_bias)
rnk$rank<-ave(1:nrow(rnk),rnk$region,FUN=seq_along)
rnk_mean<-res %>%
  filter(var=="iHosp") %>%
  filter(region %in% regions_metro) %>%
  filter(date_proj<=date_max) %>%
  group_by(model) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  mutate(region="Reg. Ave.") %>%
  arrange(region,average_bias)
rnk_mean$rank<-1:nrow(rnk_mean)
rnk<-bind_rows(rnk,rnk_mean)

x_levels<-c(regions,"Reg. Ave.")
rnk$region<-factor(rnk$region,levels=x_levels, ordered=T)

# order models by rank
rnk$model[rnk$region=="Reg. Ave."]
rnk$model2<-factor(rnk$model,levels=rev(c("Ensemble", "MLR","BRT","GAM2","ARDL","RF","ARIMA2","Const2",
                                          "Const7","PL","baseline","GAM1","ARIMA1")),ordered=T)      

p1<-ggplot(rnk, 
           aes(y=model2, fill=rank, x=region)) + 
  geom_tile() + ylab(NULL) + xlab(NULL) +
  geom_text(aes(label=rank), size=5) +
  scale_fill_gradient(low = "dodgerblue4", high="white", limits=c(1,length(models)+1)) +
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 0.5, vjust = .5),
        axis.text.y = element_text(size=12))+
  labs(tag="A")


# MAPE for all targets, national
p2<-res %>%
  filter(region=="metropolitan") %>%  
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
  coord_cartesian(ylim=c(0,0.35)) +
  theme_bw()+
  labs(tag = "B") + 
  ggtitle("MAPE, metropolitan France")

# MAPE for all targets, regional
p3<-res %>%
  filter(region %in% regions_metro) %>%
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
  coord_cartesian(ylim=c(0,0.35)) +
  theme_bw()+
  labs(tag = "C") + 
  ggtitle("MAPE, regional level")


bottom_row <- plot_grid(p2, p3, ncol=2)


ggsave(filename    = paste0(path_figures, "Figure4.png"),
  plot = plot_grid(p1, bottom_row, ncol = 1, rel_heights = c(1.3, 1)),
  width = 9, height = 8, units = "in", dpi = 300)

