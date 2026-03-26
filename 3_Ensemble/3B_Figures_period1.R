### ------------------------------------------------------------------------ ### 
###  Plot the results of the evaluation period
### ------------------------------------------------------------------------ ### 


# Prepare data  ------------------------------------------------

date_max<-max(dates_analysis)+max_prediction_horizon-remove_last_n
dat2<-data_rt %>%
  filter(date_analysis==max(data_rt$date_analysis)) %>%
  filter(date>=min(dates_analysis)-remove_last_n & date<=date_max) %>%
  dplyr::select(date,region,iHosp_smooth) %>%
  pivot_longer(cols=c("iHosp_smooth"),
               names_to="var",values_to="smooth_value") %>%
  mutate(var=gsub("_smooth","",var)) %>%
  filter(!is.element(region,c("COR","national")))



# Load projections from components models -----------------------

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

# exclusion of lockdown
proj$inclusion<-1
proj$inclusion[which(proj$date>=as.Date("2020-10-20") & proj$date<=as.Date("2020-11-04") & 
                       proj$date_proj>=as.Date("2020-11-03"))]<-0

# Compute evaluation metrics  ------------------------------------

res<-proj %>%
  left_join(dat2,by=c("date_proj"="date","region"="region","var"="var")) %>%
  mutate(in_interval95 = smooth_value >= lower_5 & smooth_value <= upper_5) %>%
  mutate(abs_error = point - smooth_value) 

# WIS 
res$IS_2<-scoringutils:::interval_score(observed = res$smooth_value,
                                        lower=res$lower_2,
                                        upper=res$upper_2,
                                        interval_range=98,
                                        weigh = TRUE,
                                        separate_results = FALSE)

res$IS_5<-scoringutils:::interval_score(observed = res$smooth_value,
                                        lower=res$lower_5,
                                        upper=res$upper_5,
                                        interval_range=95,
                                        weigh = TRUE,
                                        separate_results = FALSE)

res$IS_10<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_10,
                                         upper=res$upper_10,
                                         interval_range=90,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_20<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_20,
                                         upper=res$upper_20,
                                         interval_range=80,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_30<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_30,
                                         upper=res$upper_30,
                                         interval_range=70,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_40<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_40,
                                         upper=res$upper_40,
                                         interval_range=60,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_50<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_50,
                                         upper=res$upper_50,
                                         interval_range=50,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_60<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_60,
                                         upper=res$upper_60,
                                         interval_range=40,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_70<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_70,
                                         upper=res$upper_70,
                                         interval_range=30,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_80<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_80,
                                         upper=res$upper_80,
                                         interval_range=20,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_90<-scoringutils:::interval_score(observed = res$smooth_value,
                                         lower=res$lower_90,
                                         upper=res$upper_90,
                                         interval_range=10,
                                         weigh = TRUE,
                                         separate_results = FALSE)

res$IS_100<-abs(res$point-res$smooth_value)

res <- res %>%
  mutate(wis = (IS_2+IS_5+IS_10+IS_20+IS_30+IS_40+IS_50+
                  IS_60+IS_70+IS_80+IS_90+0.5*IS_100)/11.5)  
# from https://github.com/reichlab/covid19-forecast-evals/blob/main/code/process-zoltar-scores.R


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
  geom_line(aes(x = date_proj, y = point, color = as.factor(prediction_horizon), alpha=inclusion,group=date),size=0.8) +
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
  filter(inclusion==1) %>% 
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
         legend.key.size = unit(2, "mm"))+
  guides(color=guide_legend(ncol=1)) +
  labs(tag = "A") + 
  ggtitle("RMSE, metropolitan France")


# WIS national
p3<-res %>%
  filter(inclusion==1) %>% 
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
  filter(inclusion==1) %>% 
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
  filter(inclusion==1) %>% 
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


