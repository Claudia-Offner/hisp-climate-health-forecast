### ------------------------------------------------------------------------ ### 
###  TSGC - Time series selection (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/offne/Documents/GitHub/hisp-climate-health-forecast/')

#### DEPENDENCIES ####
source('3_TSGC/TSGC - Dependencies.R')
setwd(path_git)
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'

# Function to get a stacked time series with legend for specific start/end dates 
get_ts_stack <- function(data, startdate, enddate, plot_specs, title) {
  
  # Reusable plotting wrappers
  make_ts_plot <- function(data, column, ylab, color, startdate, enddate) {
    
    variable <- xts(data[[column]], order.by = data$date)
    
    p <- plot_moving_avg(variable, MA = 14, y.lab = ylab, colour = color, est.date.index=NULL,
                    plt.start.date = startdate, estimation.date.end = enddate)
    
    if(column =='total_patients'){
     p <- p +  scale_y_continuous(limits = c(0, 1250))
     
    }
      
    if(column =='tp'){
      res <- p + theme(axis.title.x = element_blank())
    } else {
      res <- p + theme(axis.title.x = element_blank(), axis.text.x  = element_blank(), axis.ticks.x = element_blank())
    }
    
    return(res)
    
  }
  
  make_manual_legend <- function(plot_specs, nrow = 1, top_margin = -20) {
    
    legend_labels <- sapply(plot_specs, `[[`, "label")
    legend_colors <- setNames(
      sapply(plot_specs, `[[`, "color"),
      legend_labels
    )
    
    legend_data <- data.frame(
      x = seq_along(legend_labels),
      y = seq_along(legend_labels),
      group = factor(legend_labels, levels = legend_labels)
    )
    
    ggplot(legend_data, aes(x = x, y = y, color = group)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = legend_colors, name = NULL) +
      theme_void() +
      guides(color = guide_legend(nrow = nrow)) +
      theme(
        legend.position = "bottom",
        legend.box.margin = margin(t = top_margin)
      )
  }
  
  # Generate plots
  ts_plots <- lapply(plot_specs, function(spec) {
    make_ts_plot(data = data, column = spec$col, ylab = spec$ylab, color = spec$color, 
                 startdate = startdate, enddate = enddate) })
  
  legend_only <- make_manual_legend(plot_specs)
  
  # Stack plots + legend
  ts_plot_stack <- wrap_plots(c(ts_plots, list(legend_only)), ncol = 1, heights = c(rep(1, length(ts_plots)), 0.15)) +
    plot_annotation(title=title, theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  return(ts_plot_stack)
  
}

# Function to get stacked series by location
get_ts_loc <- function (loc, df) {
  
  sub <- df %>% filter(location==loc)
  clean <- sub[!is.na(sub$total_patients), ]
  
  # Get full TS
  startdate <- min(df$date)
  enddate <- max(df$date)-30
  # clean <- clean[clean$date >= startdate & clean$date <= enddate,]
  variable <- xts(clean$total_patients,  order.by=clean$date)
  (ts_plot <- plot_moving_avg(variable, MA=14, y.lab='Total patients', colour='#8968CD', plt.start.date=startdate, estimation.date.end=enddate)+
      ggtitle(paste0(loc, ' ')) + theme(plot.title = element_text(size=16, hjust = 0.5)))
  
  # # Get param select TS (2015-2017)
  # startdate <- as.Date('2015-02-14')
  # enddate <- as.Date('2017-02-13')
  # clean_2017 <- clean[clean$date >= startdate & clean$date <= enddate,]
  # (ts_stack_2017 <- get_ts_stack(data=clean_2017, startdate=startdate, enddate=enddate, plot_specs=plot_specs,
  #                                title=paste0(loc, ' (Feb 2015 to Feb 2017)')))
  # 
  # # Get param select TS (2017-2019)
  # startdate <- as.Date('2015-02-14')
  # enddate <- as.Date('2019-02-13')
  # clean_2019 <- clean[clean$date >= startdate & clean$date <= enddate,]
  # (ts_stack_2019 <- get_ts_stack(data=clean_2019, startdate=startdate, enddate=enddate, plot_specs=plot_specs,
  #                                title=paste0(loc, ' (Feb 2015 to Feb 2019)')))
  
  return(list(ts_plot=ts_plot)) #, ts_stack_2017=ts_stack_2017, ts_stack_2019=ts_stack_2019))
  
}

####  MAIN CODE  #####

# Load data
load(paste0('3_TSGC/', file_import, '.RData'))

# NB: We are going to use Feb 2015-2017 for parameter tuning. We will use Feb 
#     2017-2019 for prediction. Chamba and Dharamshala will be used as case studies.

# Sec variables of interest with respective specifications
plot_specs <- list(
  list(col = "total_patients", label = "Total Patients", ylab = "Patients", color = "#8968CD"),
  list(col = "mx2t", label = "Max Temperature (\u00B0C)", ylab = "Celsius",  color = "#D74B3C"),
  list(col = "tp", label = "Total precipitaiton (mm)", ylab = "Millimeters",  color = "#1D98EB")
)

# 1. Time series visuals ####

ts_solan <- get_ts_loc(loc='Solan', df)
ts_nala <- get_ts_loc(loc='Nalagarh', df)
ts_shimla <- get_ts_loc(loc='Shimla', df)

#### EXPORT ####

# Figures
figs <- list('2A_ts_solan'=list(p=ts_solan$ts_plot, w=8, h=5),
             # '2A_ts_solan_stack'=list(p=ts_solan$ts_stack_2019, w=6, h=7)#,
             '2B_ts_nala'=list(p=ts_nala$ts_plot, w=8, h=5),
             '2C_ts_shimla'=list(p=ts_shimla$ts_plot, w=8, h=5)
             )

for (f in names(figs)){
  ggsave(paste0(path_report, 'III. tsgc/fig_desc_', f, '.png'), figs[[f]]$p, width = figs[[f]]$w, height = figs[[f]]$h, dpi = 300, bg='white')
}


