### ------------------------------------------------------------------------ ### 
###  Plot the results of the evaluation period
### ------------------------------------------------------------------------ ### 


# 1. Prepare data ####

# Set forecast horizon
date_max<-max(dates_analysis)+max_prediction_horizon-remove_last_n

# Format data for plotting
dat2<-data_rt %>%
  filter(date_analysis==max(data_rt$date_analysis)) %>%
  filter(date>=min(dates_analysis)-remove_last_n & date<=date_max) %>%
  dplyr::select(date,region,iHosp_smooth) %>%
  pivot_longer(cols=c("iHosp_smooth"),
               names_to="var",values_to="smooth_value") %>%
  mutate(var=gsub("_smooth","",var)) %>%
  filter(!is.element(region,c("COR","national")))


# 2. Load projections from components models ####

proj<-tibble()
for (model in models){
  proj_tmp<-read.csv2(file=paste0("3_Ensemble/Results/",model,"_period",period,".csv"),stringsAsFactors =F)
  proj<-rbind.fill(proj,proj_tmp)
}
proj$model<-factor(proj$model,levels=models, ordered=T)

proj<- proj %>%
  mutate(date=as.Date(date)) %>%
  filter(date>=min(dates_analysis))  %>%
  filter(date<=max(dates_analysis)) %>%
  mutate(date_proj=date+prediction_horizon) 

# 3. Compute evaluation metrics ####

# Error
res<-proj %>%
  left_join(dat2,by=c("date_proj"="date","region"="region","var"="var")) %>%
  mutate(in_interval95 = smooth_value >= lower_5 & smooth_value <= upper_5) %>%
  mutate(abs_error = point - smooth_value) 

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

# Plot settings ####

pal<-c("black",brewer.pal(10,"Paired"),"turquoise")
val_size<-c(rep(0.5,length(models)))
val_linetype<-c("dashed", rep("solid",length(models)-1))
lab_mod<-levels(droplevels(proj$model))

# SF4 : Forecasted trajectories of hospitalizations ####


myvar<-"iHosp"
p1<-res %>%
  filter(var==myvar) %>%
  filter(region=="metropolitan") %>% 
  filter(date_proj<=date_max) %>% 
  ggplot() +
  annotate(geom="rect",xmin=as.Date("2017-03-03"),xmax=as.Date("2017-03-18"),ymin=-Inf,ymax=Inf,fill="grey90",alpha=0.5)+
  geom_line(aes(x = date_proj, y = point, color = as.factor(prediction_horizon), group=date),size=0.8) +
  geom_line(data = dat2[dat2$region=="metropolitan" & dat2$var==myvar & dat2$date>=min(dates_analysis),], 
            aes(x = date, y = smooth_value), size = 0.8) +
  xlab("") +
  facet_wrap(~model, scale = "free", ncol=3) +
  scale_y_continuous("Hospital Admissions", limits = c(0, NA)) +
  scale_color_discrete("Prediction\nhorizon") +
  scale_alpha(range=c(0.1,1)) +
  coord_cartesian(xlim=c(min(dates_analysis),max(dates_analysis)+14)) +
  theme_bw() + 
  guides(alpha = FALSE)    


ggsave(filename=paste0(path_figures, 'FigureS4.png'), plot=p1, 
       width=12, height=12, units = "in", dpi=300, bg='white')


# MF1 : RMSE and WIS ####

# RMSE national 
p2<-res %>%
  filter(region == "metropolitan") %>% 
  filter(var=="iHosp") %>% 
  group_by(model,var,prediction_horizon) %>%
  dplyr::summarize(average_bias = sqrt(mean(abs_error^2,na.rm=T))) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = prediction_horizon, y = average_bias, color = model, size=model,linetype=model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous("RMSE") +
  scale_x_continuous("Days ahead", breaks = unique(res$prediction_horizon)) + 
  scale_color_manual("Model",labels=lab_mod,values=pal) +
  scale_size_manual("Model",labels=lab_mod,values=val_size) +
  scale_linetype_manual("Model",labels=lab_mod,values=val_linetype) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(2, "mm"))+
  guides(color=guide_legend(ncol=1)) +
  labs(tag = "A") + 
  ggtitle("RMSE, metropolitan France")


# WIS national
p3<-res %>%
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
  #coord_cartesian(ylim=c(0,0.45)) +
  theme_bw()+
  theme( legend.text=element_text(size=10),
         legend.key.size = unit(2, "mm"))+
  guides(color=guide_legend(ncol=1)) +
  labs(tag = "C") + 
  ggtitle("WIS, metropolitan France")



# RMSE regional
p2_reg<-res %>%
  filter(region %in% regions_metro) %>% 
  filter(region != "COR") %>%
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
         legend.key.size = unit(2, "mm"))+
  guides(color=guide_legend(ncol=1))+
  labs(tag = "B") + 
  ggtitle("RMSE, regional level")


# WIS regional
p3_reg<-res %>%
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
         legend.key.size = unit(2, "mm"))+
  guides(color=guide_legend(ncol=1))+
  labs(tag = "D") + 
  ggtitle("WIS, regional level")

ggsave(filename = paste0(path_figures, "Figure1.png"),
       plot = plot_grid(p2, p2_reg, p3, p3_reg, ncol = 2, align = "v", axis = "lr"),
       width = 9, height = 7, units = "in", dpi = 300)


