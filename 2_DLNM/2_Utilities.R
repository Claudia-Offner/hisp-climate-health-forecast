### ------------------------------------------------------------------------ ### 
### DLNM - Utilities
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/


#### FUNCTIONS ####

# Function to replace exact matches in df from a dictionary
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

# Comparative metrics
comp_metrics <- function(model, cb=NULL, cen=NULL) {
  qaic <- function(model) {
    c_hat <- sum(residuals(model, type = "pearson")^2) / df.residual(model)
    c_hat <- max(1, c_hat)
    k  <- length(coef(model)) + 1
    
    if (family(model)$family == "quasipoisson") {
      # logLik() is NA for quasi-families; use the equivalent Poisson log-likelihood instead
      mu <- fitted(model)
      y  <- model$y
      ll <- sum(dpois(y, mu, log = TRUE))
    } else {
      ll <- as.numeric(logLik(model))
    }
    
    -2 * ll / c_hat + 2 * k
  }
  res <- list()
  res[['qaic']]  <- qaic(model)
  res[['acf']]  <- acf(residuals(model, type='pearson'), plot=FALSE, lag.max=30)$acf
  res[['pacf']] <- pacf(residuals(model, type='pearson'), plot=FALSE, lag.max=30)$acf
  if (!is.null(cb)) {
    pred <- crosspred(cb, model, cen)
    res[['esti']] <- pred$coef[4]
  }
  return(res)
}

# Function to run model with static formula
run_model <- function(df, groups, outcome, fam_type='binom',
                      longterm_df=8, autocorr_lag=14, lag_window=c(0, 10),
                      arglag_tmp=3, argvar_tmp=1, arglag_tp=3, argvar_tp=1,
                      group_meta=NULL, subgroup_var=NULL) {
  
  # Meta-analysis
  meta_analysis <- function(grps, coefs, vcovs, arg_var, arg_lag, ranges, cen, lag_window,
                            group_meta=NULL, subgroup_var=NULL) {
    
    bound  <- colMeans(ranges)
    b_mean <- crossbasis(seq(bound[1], bound[2], length=30),
                         lag=lag_window, argvar=arg_var, arglag=arg_lag)
    
    mv <- mixmeta(coefs ~ 1, S=vcovs, method="ml")
    
    if (nrow(coefs) > 10 && !is.null(subgroup_var)) {
      
      # Set group data
      gm <- group_meta[match(rownames(coefs), group_meta$groups), , drop=FALSE]
      gm$sg <- gm[[subgroup_var]]
      
      # # Pool by group, then by location
      # mv_sub0 <- mixmeta(coefs ~ 0+sg, S=vcovs, data=gm, method="ml")
      # mv_sub1 <- update(mv_sub0, random=~1|location)
      # lrt <- lrtest(mv, mv_sub1)
      mv_sub <- mixmeta(coefs ~ 0+sg, S=vcovs, data=gm, random = ~1 | location, method="ml")
      lrt <- lrtest(mv, mv_sub)
      
      # Store
      return(list(coefs=coefs, vcovs=vcovs, arg_var=arg_var, arg_lag=arg_lag, 
                  ranges=ranges, bound=bound, b_mean=b_mean, cen=cen, 
                  mv=mv_sub, lrt=lrt, var=subgroup_var,levels=unique(gm[[subgroup_var]])))
    } else{
      
      return(list(coefs=coefs, vcovs=vcovs, arg_var=arg_var, arg_lag=arg_lag, 
                  ranges=ranges, bound=bound, b_mean=b_mean, cen=cen, 
                  mv=mv, lrt=NULL, var=NULL, levels=NULL))
    }
    
  }
  
  # 1. Data and dependency prep ----
  datalist <- lapply(groups, function(g) df[df$groups == g, ])
  names(datalist) <- groups
  datalist <- Filter(function(x) nrow(x) > 0, datalist)
  
  grps_used <- names(datalist)
  m     <- length(datalist)
  years <- length(unique(df$year))
  
  # TMP
  ranges_tmp <- t(sapply(datalist, function(x) range(x$tmp, na.rm=TRUE)))
  cen_tmp    <- sapply(datalist, function(x) quantile(x$tmp, c(0:100)/100, na.rm=TRUE))
  arglag_tmp <- list(fun='ns', df=arglag_tmp)
  argvar_tmp <- if (argvar_tmp == 1) list(fun="lin") else list(fun="ns", df=argvar_tmp)
  df_tmp     <- prod(arglag_tmp$df, argvar_tmp$df)
  
  # TP
  ranges_tp <- t(sapply(datalist, function(x) range(x$tp, na.rm=TRUE)))
  cen_tp    <- sapply(datalist, function(x) quantile(x$tp, c(0:100)/100, na.rm=TRUE))
  arglag_tp <- list(fun='ns', df=arglag_tp)
  argvar_tp <- if (argvar_tp == 1) list(fun="lin") else list(fun="ns", df=argvar_tp)
  df_tp     <- prod(arglag_tp$df, argvar_tp$df)
  
  # 2. Stage 1: Fit the DLNM ----
  Slist_model  <- setNames(vector('list', m), grps_used)
  metric_comp  <- setNames(vector('list', m), grps_used)

  coefs_tmp  <- matrix(NA, m, df_tmp, dimnames=list(grps_used, paste0('spl', seq(df_tmp))))
  vcovs_tmp <- setNames(vector('list', m), grps_used)
  
  coefs_tp   <- matrix(NA, m, df_tp, dimnames=list(grps_used, paste0('spl', seq(df_tp))))
  vcovs_tp  <- setNames(vector('list', m), grps_used)
  
  # Suppress warnings for prediction beyond boundaries; reset on exit
  old_warn <- options(warn=-1)
  on.exit(options(old_warn), add=TRUE)
  
  for (i in seq(m)) {
    cat(i, '')
    sub <- datalist[[i]]
    
    # Autocorrelation lags
    if (autocorr_lag == 0) {
      lags <- ''
    } else {
      lag_terms <- paste0(outcome, '_lag', 1:autocorr_lag)
      lags <- paste0('+ ', paste(lag_terms, collapse=' + '))
    }
    
    # Cross-bases
    cb.tmp <- crossbasis(sub$tmp, lag=lag_window, argvar=argvar_tmp, arglag=arglag_tmp)
    cb.tp  <- crossbasis(sub$tp,  lag=lag_window, argvar=argvar_tp,  arglag=arglag_tp)
    controls <- '+ weekday + holiday + rh + fg10 + pm25'
    
    form <- as.formula(
      paste0(outcome, " ~ cb.tmp + cb.tp + ns(time, df=", longterm_df, "*", years, ')',
             lags, controls)
    )
    
    # Fit
    if (fam_type == 'quasi') {
      model <- glm(formula=form, family=quasipoisson(), data=sub)
    } else if (fam_type == 'binom') {
      model <- glm.nb(formula=form, data=sub)
    } else {
      stop('Unknown model type: ', fam_type)
    }
    
    # Extract coef and vcov
    predtmp <- crosspred(cb.tmp, model, cen=median(sub$tmp, na.rm=TRUE))
    coefs_tmp[i, ]   <- predtmp$coef
    vcovs_tmp[[i]]  <- predtmp$vcov
    
    predtp <- crosspred(cb.tp, model, cen=median(sub$tp, na.rm=TRUE))
    coefs_tp[i, ]    <- predtp$coef
    vcovs_tp[[i]]   <- predtp$vcov
    
    # Comparative metrics
    Slist_model[[i]] <- model
    metric_comp[[i]]  <- comp_metrics(model)
    metric_comp[[i]]$esti_tmp <- unname(predtmp$coefficients[1])
    metric_comp[[i]]$esti_tp  <- unname(predtp$coefficients[1])
  }
  
  # 3. Stage 2: Multivariate meta-analysis ----
  meta_comp <- list()
  meta_comp[['avg_tmp']]  <- mean(sapply(metric_comp, function(x) x$esti_tmp))
  meta_comp[['avg_tp']]   <- mean(sapply(metric_comp, function(x) x$esti_tp))
  meta_comp[['sum_qaic']]  <- sum(sapply(metric_comp, function(x) x$qaic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(metric_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']]  <- rowMeans(sapply(metric_comp, function(x) as.numeric(x$acf)))
  
  meta_tmp <- meta_analysis(grps=grps_used, coefs=coefs_tmp, vcovs=vcovs_tmp,
                            arg_var=argvar_tmp, arg_lag=arglag_tmp,
                            ranges=ranges_tmp, cen=cen_tmp, lag_window=lag_window,
                            group_meta=group_meta, subgroup_var=subgroup_var)
  
  meta_tp <- meta_analysis(grps=grps_used, coefs=coefs_tp, vcovs=vcovs_tp,
                           arg_var=argvar_tp, arg_lag=arglag_tp,
                           ranges=ranges_tp, cen=cen_tp, lag_window=lag_window,
                           group_meta=group_meta, subgroup_var=subgroup_var)
  
  list(models=Slist_model, meta_comp=meta_comp, meta_tmp=meta_tmp, meta_tp=meta_tp)
}


# Function to extract predictions
meta_predict <- function(meta, city='ALL', reference='50%', blup=FALSE){
  list2env(meta, envir=environment())
  perc <- as.numeric(sub("%", "", names(rowMeans(cen))))
  reference <- if(reference=='0%') paste0(min(perc), '%') else reference
  
  if (city=='ALL') {
    at <- rowMeans(cen)
    cp <- crosspred(b_mean, coef=coef(mv), vcov=vcov(mv), model.link='log', 
                    at=at, cen=at[reference])
  } else {
    i <- which(rownames(coefs) == city)
    if(blup==TRUE){
      blup <- blup(mv, vcov=TRUE)
      coef=blup[[i]]$blup; vcov=blup[[i]]$vcov
    } else{
      coef=coefs[i, ]; vcov=vcovs[[i]]
    }
    at <- cen[, i]
    cp <- crosspred(b_mean, coef=coef, vcov=vcov, model.link='log', from=ranges[i, 1], to=ranges[i, 2], 
                    at=at, cen=at[reference])
  }
  
  return(list(cp=cp, perc=perc, at=at))
}

# Function to plot pooled & city specific exposure-response relationships, cumulatively and by lag
city_exp_res <- function(city='ALL', meta, col='lightblue', xlab=NULL, lag=NULL, title=NULL, reference='50%'){
  
  # FUNCTION: Formatting for multi-city plots
  plot_rr_curve <- function(city=city, df_overall, df_specific1=NULL, df_specific2=NULL, xlab, col, title=NULL, legend='Hospital-specific') {
    
    if (is.null(df_specific1) && is.null(df_specific2)) {
      
      ggplot() +
        geom_hline(yintercept=1, color='gray20') +
        # geom_line(data=df_specific1, aes(x=pred, y=RR, group=city, color=legend, linetype=legend), size=0.7) +
        geom_ribbon(data=df_overall, aes(x=pred, ymin=low, ymax=high), fill=col, alpha=0.2) +
        geom_line(data=df_overall, aes(x=pred, y=RR, color='Hospital-specific', linetype='Hospital-specific'), size=1.2) +
        scale_color_manual(name=NULL, values=c('Hospital-specific'=col)) +
        scale_linetype_manual(name=NULL, values=c('Hospital-specific'='solid')) +
        scale_x_continuous(breaks=pretty(df_overall$pred, n=5)) +
        scale_y_continuous(breaks=pretty(c(min(df_overall$low), max(df_overall$high)), n=5)) +
        labs(x=xlab, y='RR', title=title) +
        coord_cartesian(ylim=c(min(df_overall$low), max(df_overall$high)), xlim=range(df_overall$pred)) +
        theme_classic() +
        theme(
          axis.line=element_line(color='gray20'),
          axis.text=element_text(color='gray20'),
          axis.ticks=element_line(color='gray20'),
          plot.title=element_text(face='bold', hjust=0.5, size=12),
          legend.position='bottom'
        )
      
    } else {
      
      if(city=='ALL'){
        
        ggplot() +
          geom_hline(yintercept=1, color='gray20') +
          geom_line(data=df_specific1, aes(x=pred, y=RR, group=city, color=legend, linetype=legend), size=0.7) +
          geom_ribbon(data=df_overall, aes(x=pred, ymin=low, ymax=high), fill=col, alpha=0.2) +
          geom_line(data=df_overall, aes(x=pred, y=RR, color='Pooled', linetype='Pooled'), size=1.2) +
          scale_color_manual(name=NULL, values=setNames(c('grey70', col), c(legend, 'Pooled'))) +
          scale_linetype_manual(name=NULL, values=setNames(c('dashed', 'solid'), c(legend, 'Pooled'))) +
          scale_x_continuous(breaks=pretty(df_overall$pred, n=5)) +
          scale_y_continuous(breaks=pretty(c(min(df_overall$low), max(df_overall$high)), n=5)) +
          labs(x=xlab, y='RR', title=title) +
          coord_cartesian(ylim=c(min(df_overall$low), max(df_overall$high)), xlim=range(df_overall$pred)) +
          theme_classic() +
          theme(
            axis.line=element_line(color='gray20'),
            axis.text=element_text(color='gray20'),
            axis.ticks=element_line(color='gray20'),
            plot.title=element_text(face='bold', hjust=0.5, size=12),
            legend.position='bottom'
          )
        
      } else{
        
        # Combine data depending on whether df_specific2 exists
        if (is.null(df_specific2)) {
          df_all <- rbind(
            data.frame(pred=df_overall$pred, RR=df_overall$RR, low=df_overall$low, high=df_overall$high, type='Pooled'),
            data.frame(pred=df_specific1$pred, RR=df_specific1$RR, low=df_specific1$low, high=df_specific1$high, type='BLUP')
          )
          
          color_vals <- c('Pooled'='black', 'BLUP'=col)
          linetype_vals <- c('Pooled'='dashed', 'BLUP'='solid')
          size_vals <- c('Pooled'=0.7, 'BLUP'=1.2)
          
        } else {
          df_all <- rbind(
            data.frame(pred=df_overall$pred, RR=df_overall$RR, low=df_overall$low, high=df_overall$high, type='Pooled'),
            data.frame(pred=df_specific1$pred, RR=df_specific1$RR, low=df_specific1$low, high=df_specific1$high, type='Raw'),
            data.frame(pred=df_specific2$pred, RR=df_specific2$RR, low=df_specific2$low, high=df_specific2$high, type='BLUP')
          )
          
          color_vals <- c('Pooled'='black', 'Raw'=col, 'BLUP'='black')
          linetype_vals <- c('Pooled'='solid', 'Raw'='solid', 'BLUP'='dashed')
          size_vals <- c('Pooled'=0.7, 'Raw'=1.2, 'BLUP'=0.7)
        }
        
        # Create plot
        ggplot() +
          geom_hline(yintercept=1, color='gray20') +
          geom_ribbon(data=df_specific1, aes(x=pred, ymin=low, ymax=high), fill=col, alpha=0.2) +
          geom_line(data=df_all, aes(x=pred, y=RR, color=type, linetype=type, size=type)) +
          scale_color_manual(name=NULL, values=color_vals) +
          scale_linetype_manual(name=NULL, values=linetype_vals) +
          scale_size_manual(name=NULL, values=size_vals) +
          scale_x_continuous(breaks=scales::breaks_pretty(n=6)) +
          scale_y_continuous(breaks=pretty(c(min(df_all$low), max(df_all$high)), n=6)) +
          labs(x=xlab, y='RR', title=title) +
          coord_cartesian(ylim=c(min(df_all$low), max(df_all$high)), xlim=range(df_overall$pred)) +
          theme_classic() +
          theme(
            axis.line=element_line(color='gray20'),
            axis.text=element_text(color='gray20'),
            axis.ticks=element_line(color='gray20'),
            plot.title=element_text(face='bold', hjust=0.5, size=12),
            legend.position='bottom'
          )
      }
      
    }
    
  }
  
  # city_exp_res — index matRRfit/allRRfit via match(at, cp$predvar), not positionally
  all       <- meta_predict(meta, city='ALL', reference=reference)
  city_raw  <- map(locations, ~ meta_predict(meta, city=.x, reference=reference)) %>% setNames(locations)
  city_blup <- map(locations, ~ meta_predict(meta, city=.x, reference=reference, blup=TRUE)) %>% setNames(locations)
  
  idx_all <- match(all$at, all$cp$predvar)
  
  if (!is.null(lag)) {
    df_overall <- data.frame(pred=all$perc, RR=all$cp$matRRfit[idx_all, lag + 1],
                             low=all$cp$matRRlow[idx_all, lag + 1], high=all$cp$matRRhigh[idx_all, lag + 1], row.names=NULL)
    df_city_raw <- map_dfr(city_raw, function(x) {
      idx <- match(x$at, x$cp$predvar)
      data.frame(pred=x$perc, RR=x$cp$matRRfit[idx, lag + 1],
                 low=x$cp$matRRlow[idx, lag + 1], high=x$cp$matRRhigh[idx, lag + 1], row.names=NULL)
    }, .id="city")
    df_city_blup <- map_dfr(city_blup, function(x) {
      idx <- match(x$at, x$cp$predvar)
      data.frame(pred=x$perc, RR=x$cp$matRRfit[idx, lag + 1],
                 low=x$cp$matRRlow[idx, lag + 1], high=x$cp$matRRhigh[idx, lag + 1], row.names=NULL)
    }, .id="city")
  } else {
    df_overall <- data.frame(pred=all$perc, RR=all$cp$allRRfit[idx_all],
                             low=all$cp$allRRlow[idx_all], high=all$cp$allRRhigh[idx_all], row.names=NULL)
    df_city_raw <- map_dfr(city_raw, function(x) {
      idx <- match(x$at, x$cp$predvar)
      data.frame(pred=x$perc, RR=x$cp$allRRfit[idx], low=x$cp$allRRlow[idx], high=x$cp$allRRhigh[idx], row.names=NULL)
    }, .id="city")
    df_city_blup <- map_dfr(city_blup, function(x) {
      idx <- match(x$at, x$cp$predvar)
      data.frame(pred=x$perc, RR=x$cp$allRRfit[idx], low=x$cp$allRRlow[idx], high=x$cp$allRRhigh[idx], row.names=NULL)
    }, .id="city")
  }
  
  p_raw  <- plot_rr_curve(city=city, df_specific1=df_city_raw, df_specific2=df_city_blup, df_overall=df_overall, xlab=xlab, col=col, title=title)
  p_blup <- plot_rr_curve(city=city, df_specific1=df_city_blup, df_overall=df_overall, xlab=xlab, col=col, title=title, legend='Hospital-specific (BLUP)')
  
  return(list(df_overall=df_overall, df_city_raw=df_city_raw, df_city_blup=df_city_blup,
              p_raw=p_raw, p_blup=p_blup))
  
}

# Function to plot overall & city specific 2D contour plots
city_con_res <- function(city='ALL', meta, xlab=NULL, title=NULL, reference='50%', 
                         sig=TRUE, shared_lims=NULL, shared_brks=NULL) {
  
  plot_contour <- function(df, xlab="", title="", ncol=300, nbreaks=15, contour_size=0.35, 
                             sig=TRUE, shared_lims=NULL, shared_brks=NULL) {
    
    pal <- colorRampPalette(c('#1D98EB', 'white', '#D74B3C'))(ncol)
    brks <- pretty(range(df$RR, na.rm=TRUE), n=nbreaks)
    
    if (!is.null(shared_lims) && !is.null(shared_brks)) {
      lims <- shared_lims
      brks <- shared_brks
    } else {
      max_dev <- max(abs(brks - 1))
      lims <- c(1 - max_dev, 1 + max_dev)
      brks <- pretty(lims, n=nbreaks)
    }
    
    p <- ggplot(df, aes(x=pred, y=lag, fill=RR)) +
      geom_raster(interpolate=TRUE) +
      scale_fill_gradient2(
        low='#1D98EB', mid='white', high='#D74B3C',
        midpoint=1, limits=lims, breaks=brks, name="RR",
        na.value="grey85",
        guide=guide_colorbar(
          barwidth=unit(12, "lines"),
          barheight=unit(0.6, "lines"),
          ticks=TRUE,
          frame.colour="black",
          ticks.colour="black",
          title.position="top",
          title.hjust=0.5,
          direction="horizontal"
        )
      ) +
      scale_x_continuous(breaks=pretty(df$pred, n=5)) +
      scale_y_continuous(breaks=pretty(df$lag, n=6)) +
      labs(x=xlab, y="Lag (days)", title=title) +
      coord_cartesian(expand=FALSE) +
      theme_classic(base_size=13) +
      theme(
        panel.grid=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=0.8),
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        axis.ticks=element_line(colour="black"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10, angle=45, hjust=1),
        legend.position='bottom',
        plot.title=element_text(size=12, face="bold", hjust=0.5)
      )
    
    
    # Overlay significance boundary if requested
    if (sig) {
      
      p <- p +
        geom_tile_pattern(
          data=subset(df, sig_dir == 'neg'),
          aes(x=pred, y=lag), pattern="crosshatch", pattern_fill="#12669e",
          pattern_colour='#12669e', fill=NA, color=NA, 
          pattern_density=0.05, pattern_spacing=0.05,
          pattern_angle=45, pattern_alpha=0.2, pattern_res=200,
          width=1, height=1
        ) +
        geom_tile_pattern(
          data=subset(df, sig_dir == 'pos'),
          aes(x=pred, y=lag), pattern="crosshatch", pattern_fill="#a0392e",
          pattern_colour="#a0392e", fill=NA, color=NA, 
          pattern_density=0.05, pattern_spacing=0.05,
          pattern_angle=45, pattern_alpha=0.2, pattern_res=200,
          width=1, height=1 
        )
    }
    p
  }
  
  get_df <- function(cp, perc, at) {
    idx <- match(at, cp$predvar)
    df <- expand.grid(pred=perc, lag=seq(cp$lag[1], cp$lag[2]))
    df$RR      <- as.vector(cp$matRRfit[idx, ])
    df$low     <- as.vector(cp$matRRlow[idx, ])
    df$high    <- as.vector(cp$matRRhigh[idx, ])
    df$sig_dir <- with(df, ifelse((low > 1 | high < 1) & abs(RR - 1) > 0.01, ifelse(RR > 1, "pos", "neg"), "none"))
    df
  }

  if (city == 'ALL') {
    pred <- meta_predict(meta, city=city, reference=reference, blup=FALSE)
    # plot_df <- get_df(pred$cp, pred$perc)
    plot_df <- get_df(pred$cp, pred$perc, pred$at)
    p <- plot_contour(df=plot_df, xlab=xlab, title=title, sig=sig, shared_lims=shared_lims, shared_brks=shared_brks)
  } else {
    pred <- meta_predict(meta, city=city, reference=reference, blup=FALSE)
    # plot_df <- get_df(pred$cp, pred$perc)
    plot_df <- get_df(pred$cp, pred$perc, pred$at)
    p <- plot_contour(df=plot_df, xlab=xlab, title=city, sig=sig, shared_lims=shared_lims, shared_brks=shared_brks)
  }
  
  return(list(df=plot_df, p=p))
  
}


#### DEPENDENCIES ####

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD', '#0D2C6E', '#910030FF', '#127804FF', '#52C5B0FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9', '#4462A2', '#91003080', '#12780480', '#52C5B080') 

w_cols <- c('t2m', 'mn2t','mx2t', 'dt', 'wb', 'tp', 'fg10', 'rh', 'pm25')
h_cols <- c('total_patients',
            'F','M','under_5','age5_64','over_65',
            'general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
            'gynaecology', 'orthopaedics', 'dermatology', 'surgery',
            'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
            'ncd', 'antiretroviral', 'substance_counselling', 'physiotherapy')
p_cols <- c('elev', 'popu', 'hosp_n', 'insurance', 'decision', 'meducation', 'hsize', 'hwealth')

location_levels <- c("Una", 'Kangra', "Nahan", "Hamirpur", 'Mandi', "Chamba", "Dharamshala", "Solan", "Kullu",  "Reckong Peo")
season_levels   <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
month_levels    <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
weekday_levels  <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
age_levels      <- c('under_5', 'age5_64',  'over_65','Unknown')
gender_levels   <- c('F', 'M', 'O')
ward_levels     <- c('general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
                     'gynaecology', 'orthopaedics', 'dermatology', 'surgery',
                     'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
                     'ncd', 'mch', 'antiretroviral', 'substance_counselling', 'physiotherapy')

# Dictionary of names for figures and tables
nm <- c(season='Season', closed='Closed', holiday='Holiday', total_patients='Total patients', 
        location='Location', weekday='Weekday', time='Time', 
        patients_lag1='Total Patients (Lag=1)', patients_lag7='Total Patients (Lag=7)',
        # POPULATION GROUPS
        M='Men', F='Women',
        under_5='Under 5', age5_64='Ages 5-64', over_65='Ages 65+', 
        # WARDS
        general='General', emergency='Emergency', ophthalmology='Ophthalmology', 
        paediatrics='Paediatric', orthopaedics='Orthopaedics', ent='ENT', 
        gynaecology='Gynaecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        pulmonary='Pulmonary', disability='Disability', ncd='NCD', mch='MCH',
        antiretroviral='HIV Treatment', substance_counselling='Substance counselling',
        physiotherapy='Physiotherapy', 
        # CONTINUOUS WEATHER
        tmp='Temperature (\u00B0C)', t2m='Mean temp. (\u00B0C)', mx2t='Max temp. (\u00B0C)', mn2t='Min temp. (\u00B0C)',
        rh='Rel. humidity (%)', tp='Precipitation (mm)', fg10='Windspeed (km/hr)', pm25='PM 2.5 (µg/m³)',
        dt="Diurnal range (\u00B0C)", wb="Wet-bulb temp. (\u00B0C)",
        # Controls
        lat= 'Latitude', lon= 'Longitude',elev= "Altitude (m)",
        popu= "Population count (per 1km)",hosp_n= "Number of hospitals", #hosp_teach= "Teaching hospital",
        insurance= "No health insurance (%)",decision= "No female decision-making (%)",
        meducation= "No maternal education (%)",hsize= "Large household size (%)",
        hwealth= "Poorest households (%)", mcv1= "MCV1 vaccination coverage",
        dtp1= "DTP1 vaccination coverage",dtp3= "DTP3 vaccination coverage",
        anc= "Antenatal care coverage",malnutrition="Malnutrition prevalence",
        pnc= "Postnatal care coverage",sba= "Skilled birth attendance",
        tt= "Tetanus toxoid coverage"
        )
        # # IMD DEFINITIONS
        # mx2t_IMD='IMD: Mx temp. >30\u00B0C', mn2t_IMD='IMD: Min temp. <0\u00B0C', 
        # tp_IMD='IMD: Precipitation >15mm', fg10_IMD='IMD: Wind speed >40km/h', pm25_IMD='IMD: PM2.5 >60 µg/m³', 
        # rh_low_IMD='IMD: Humidity <30%', rh_hig_IMD='IMD: Humidity >60%',
        # # CLIM. DEFINITIONS
        # mx2t_p90='Clim: Max temp. >90th', mn2t_p10='Clim: Min temp. <10th', 
        # tp_p90='Clim: Precipitation >90th', fg10_p90='Clim: Wind speed >90th', pm25_p90='Clim: PM2.5 >90th', 
        # rh_mean_p10 ='Clim: Humidity <10th', rh_mean_p90='Clim: Humidity >90th') 




# # Function to get usable tables and figures
# get_results <- function(output, city='ALL', lags=NULL, probs='50%', outcome){
#   
#   # Function to extract model comparison criteria
#   ftab <- function(m, mref=NULL) {
#     
#     # RETURN NOTHING IF m IS NULL
#     if (is.null(m)) return(rep(NA_real_, 12))
#     
#     # HETEROGENEITY AND IC STATS
#     q <- qtest(m)
#     het <- c(q$Q[1],q$df[1],q$pvalue[1],(q$Q[1]-q$df[1])/q$Q[1]*100,AIC(m),BIC(m))
#     # LR TEST (ONLY FOR META-REGRESSION)
#     if(!is.null(mref)) {
#       lrstat <- -2*(logLik(mref)-logLik(m))[1]
#       df <- attr(logLik(m),'df')-attr(logLik(mref),'df')
#       pvalue <- 1-pchisq(lrstat,df)
#       lr <- c(lrstat,df,pvalue)
#     }
#     # WALD TEST (ONLY FOR META-REGRESSION)
#     if(!is.null(mref)) {
#       coef <- coef(m)[-grep('Int',names(coef(m)))]
#       vcov <- vcov(m)[-grep('Int',names(coef(m))),-grep('Int',names(coef(m)))]
#       waldstat <- coef%*%solve(vcov)%*%coef
#       df <- length(coef)
#       pvalue <- 1-pchisq(waldstat,df)
#       wald <- c(waldstat,df,pvalue)
#     }
#     # outputS
#     if(!is.null(mref)) {
#       return(c(het,lr,wald))
#     } else return(c(het,rep(NA,6))) 
#   }
#   
#   meta_tmp <- output$meta_tmp
#   meta_tp <- output$meta_tp
#   
#   # Model selection criteria
#   tab <- matrix(NA,2,12)
#   colnames(tab) <- c('Q','df','p','I-square','AIC','BIC','stat','df','p','stat','df','p')
#   rownames(tab) <- c('tmp','tp') #, 'fg10', 'rh', 'fg10', 'pm25')
#   tab[1,] <- ftab(m=meta_tmp$mv)
#   tab[2,] <- ftab(m=meta_tp$mv)
#   
#   # Exposure response (with city specs) by lag
#   p_er_tmp <- city_exp_res(city=city, lag=lags, meta=meta_tmp, col='#D74B3C', xlab='Max. Temperature (%)', title=NULL)
#   p_er_tp <- city_exp_res(city=city, lag=lags, meta=meta_tp, col='#1D98EB', xlab='Precipitation (%)', title=NULL, reference='0%')
#   
#   # Lag response (with city specs) by percentile
#   q_tmp <- rowMeans(meta_tmp$cen)[probs]
#   q_tp <- rowMeans(meta_tp$cen)[probs]
#   p_lr_tmp <- city_lag_res(city=city, quant_v=q_tmp, meta=meta_tmp, col='#D74B3C', xlab='Max. Temperature (%)', title=NULL)
#   p_lr_tp <- city_lag_res(city=city, quant_v=q_tp, meta=meta_tp, col='#1D98EB', xlab='Precipitation (%)', title=NULL, reference='0%')
#   
#   return(list(tab=tab, 
#               p_er_tmp=p_er_tmp, p_er_tp=p_er_tp, #p_er_fg10=p_er_fg10, p_er_rh=p_er_rh, p_er_pm25=p_er_pm25, 
#               p_lr_tmp=p_lr_tmp, p_lr_tp=p_lr_tp)) #, p_lr_fg10=p_lr_fg10, p_lr_rh=p_lr_rh, p_lr_pm25=p_lr_pm25)) 
#   
# }
# 
# 
# # Function to extract and format relevant results and figures
# format_results <- function(output, outcome, city=NULL){
#   
#   # OVERALL
#   overall <- get_results(output, city='ALL', lags=NULL, outcome=outcome)
#   # LAG
#   quant_test <- c('10%', '25%', '75%', '90%')
#   quant <- list()
#   for (q in quant_test){
#     quant[[q]] <- get_results(output, city='ALL', lags=NULL, probs=q, outcome=outcome)
#   }
#   # EXPOSURE
#   # lag_test <- c(0, 3, 7, 10)
#   # lags <- list()
#   # for (l in lag_test){
#   #   lags[[as.character(l)]] <- get_results(output, city='ALL', lags=l, outcome=outcome)
#   # }
#   # 
#   if(!is.null(output$meta_tmp$mv)){
#     
#     # Format overall plots
#     all <- overall[!names(overall) %in% 'tab']
#     p_blup_all <- lapply(all, function(x) x[['p_blup']])
#     p_blup_er <- p_blup_all[grepl('er', names(p_blup_all))] #|> head(-2)
#     p_blup_lr <- p_blup_all[grepl('lr', names(p_blup_all))] #|> head(-2)
#     p_er_overall <- wrap_plots(p_blup_er, nrow=1)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#     
#     # Format lag plots
#     p_lr_tmp <- wrap_plots(lapply(quant, function(x) x$p_lr_tmp$p_blup), ncol=4) & theme(legend.position='none')
#     p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     p_lr <- wrap_plots(p_lr_tmp, p_lr_tp, nrow=2)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#     
#     # Format exposure plots
#     # p_er_tmp <- wrap_plots(lapply(lags, function(x) x$p_er_tmp$p_blup), ncol=4) & theme(legend.position='none')
#     # p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     # p_er <- wrap_plots(p_er_tmp, p_er_tp, nrow=2)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#     
#     if(!is.null(city)){
#       
#       locations <- city
#       city <- list()
#       city_er <- list()
#       city_lr <- list()
#       for (i in locations) {
#         city_er[[i]] <- list()
#         city_lr[[i]] <- list()
#         for (l in lag_test) {
#           city_er[[i]][[as.character(l)]] <- get_results(output, city=i, lags=l, outcome=outcome)
#           print(paste0(i, ' ', l))
#         }
#         for (q in quant_test) {
#           city_lr[[i]][[as.character(q)]] <- get_results(output, city=i, probs=q, outcome=outcome)
#           print(paste0(i, ' ', q))
#         }
#         
#         p_er_tmp_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_tmp$p_city + theme(legend.position="none")), ncol=4)+
#           plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#         p_er_tp_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_tp$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#           plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#         # p_er_fg10_city <- wrap_plots(lapply(city_er[[i]], function(x) x$p_er_fg10$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#         #   plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#  
#         p_lr_tmp_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_tmp$p_city + theme(legend.position="none")), ncol=4)+
#           plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#         p_lr_tp_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_tp$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#           plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#         # p_lr_fg10_city <- wrap_plots(lapply(city_lr[[i]], function(x) x$p_lr_fg10$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#         #   plot_annotation(title=city, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#   
#       }
#       
#     } else {
#       
#       city <- list()
#       city_er <- list()
#       city_lr <- list()
#       for (i in locations){
#         city[[i]] <- get_results(output, city=i, lag=NULL, probs='50%', outcome=outcome)
#       }
#       p_er_tmp_city <- wrap_plots(lapply(city, function(x) x$p_er_tmp$p_city + theme(legend.position="none")), ncol=4)+
#         plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       p_er_tp_city  <- wrap_plots(lapply(city, function(x) x$p_er_tp$p_city + theme(legend.position="none")), ncol=4)+
#         plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       # p_er_fg10_city  <- wrap_plots(lapply(city, function(x) x$p_er_fg10$p_city + theme(legend.position="none")), ncol=4)+
#       #   plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       # p_er_rh_city <- wrap_plots(lapply(city, function(x) x$p_er_rh$p_city + theme(legend.position="none")), ncol=4)+
#       #   plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       # p_er_pm25_city  <- wrap_plots(lapply(city, function(x) x$p_er_pm25$p_city + theme(legend.position="none")), ncol=4)+
#       #   plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       
#       p_lr_tmp_city  <- wrap_plots(lapply(city, function(x) x$p_lr_tmp$p_city + theme(legend.position="none")), ncol=4)+
#         plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       p_lr_tp_city <- wrap_plots(lapply(city, function(x) x$p_lr_tp$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#         plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       # p_lr_fg10_city <- wrap_plots(lapply(city, function(x) x$p_lr_fg10$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#       #   plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       # p_lr_rh_city  <- wrap_plots(lapply(city, function(x) x$p_lr_rh$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#       #   plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       # p_lr_pm25_city <- wrap_plots(lapply(city, function(x) x$p_lr_pm25$p_city + theme(legend.position="none", plot.title=element_blank())), ncol=4)+
#       #   plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#       
#     }
#   
#     return(list(overall=overall, p_er_overall=p_er_overall, 
#                 quant=quant, p_lr=p_lr, lags=lags, p_er=p_er))
#                 # city=city, city_er=city_er, city_lr=city_lr,
#                 # p_er_tmp_city=p_er_tmp_city, p_er_tp_city=p_er_tp_city, # p_er_fg10_city=p_er_fg10_city, #p_er_rh_city=p_er_rh_city, p_er_pm25_city=p_er_pm25_city,
#                 # p_lr_tmp_city=p_lr_tmp_city, p_lr_tp_city=p_lr_tp_city)) #p_lr_fg10_city=p_lr_fg10_city)) #p_lr_rh_city=p_lr_rh_city, p_lr_pm25_city=p_lr_pm25_city))
#                 # 
#     } else {
#   
#     # Format overall plots
#     all <- overall[!names(overall) %in% 'tab']
#     p_raw_all <- lapply(all, function(x) x[['p_raw']])
#     p_raw_er <- p_raw_all[grepl('er', names(p_raw_all))] #|> head(-2)
#     p_raw_lr <- p_raw_all[grepl('lr', names(p_raw_all))] #|> head(-2)
#     p_er_overall <- wrap_plots(p_raw_er, nrow=1)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#     
#     # Format lag plots
#     p_lr_tmp <- wrap_plots(lapply(quant, function(x) x$p_lr_tmp$p_raw), nrow=1) & theme(legend.position='none')
#     p_lr_tp <- wrap_plots(lapply(quant, function(x) x$p_lr_tp$p_raw), nrow=1) & theme(legend.position='none', plot.title=element_blank())
#     # p_lr_fg10 <- wrap_plots(lapply(quant, function(x) x$p_lr_fg10$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     # p_lr_rh <- wrap_plots(lapply(quant, function(x) x$p_lr_rh$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     # p_lr_pm25 <- wrap_plots(lapply(quant, function(x) x$p_lr_pm25$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     p_lr <- wrap_plots(p_lr_tmp, p_lr_tp, nrow=2)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#     
#     # Format exposure plots
#     p_er_tmp <- wrap_plots(lapply(lags, function(x) x$p_er_tmp$p_raw), nrow=1) & theme(legend.position='none')
#     p_er_tp <- wrap_plots(lapply(lags, function(x) x$p_er_tp$p_raw), nrow=1) & theme(legend.position='none', plot.title=element_blank())
#     # p_er_fg10 <- wrap_plots(lapply(lags, function(x) x$p_er_fg10$p_raw), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     # p_er_rh <- wrap_plots(lapply(lags, function(x) x$p_er_rh$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     # p_er_pm25 <- wrap_plots(lapply(lags, function(x) x$p_er_pm25$p_blup), ncol=4) & theme(legend.position='none', plot.title=element_blank())
#     p_er <- wrap_plots(p_er_tmp, p_er_tp, nrow=2)+plot_annotation(title=outcome, theme=theme(plot.title=element_text(hjust=0.5, size=16, face='bold')))
#   
#     return(list(overall=overall, quant=quant, lags=lags, 
#                 p_er_overall=p_er_overall, p_lr=p_lr, p_er=p_er))
#   }
# }