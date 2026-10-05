### ------------------------------------------------------------------------ ### 
### DLNM - Utilities
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/


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



# Function to run model with static formula
run_model <- function(df, groups, outcome, fam_type='binom',
                      longterm_df=8, autocorr_lag=14, lag_window=c(0,10),
                      arglag_tmp=3, argvar_tmp=3, arglag_tp=4, argvar_tp=1){
  
  # Function to calculate comparative metrics for meta-DLNM
  comp_metrics <- function(model, cb=NULL, cen=NULL) {
    
    qaic <- function(model) {
      c_hat <- sum(residuals(model, type = "pearson")^2) / df.residual(model)
      c_hat <- max(1, c_hat)
      k <- length(coef(model)) + 1 
      ll <- as.numeric(logLik(model))
      -2 * ll / c_hat + 2 * k
    }
    res <- list()
    # Get sum qaic/ACF/PACF
    # res[['aic']] <- summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion
    res[['aic']] <- qaic(model)
    res[['acf']] <- acf(residuals(model, type='pearson'), plot=FALSE, lag.max=30)$acf
    res[['pacf']] <- pacf(residuals(model, type='pearson'), plot=FALSE, lag.max=30)$acf
    
    if(!is.null(cb)){
      # Get average of estimates (lag 4)
      pred <- crosspred(cb, model, cen)
      res[['esti']] <- pred$coef[4]  
    } 
    
    return(res) 
  }
  
  # Function to get reference meta-analysis
  meta_analysis <- function(groups, ymat, Slist, arg_var, arg_lag, ranges, cen, lag_window){
    
    if(length(groups) > 10){
      m_age <- mixmeta(ymat ~ 0 + gender, S = Slist, random = ~ 1 | location,
                       data = meta, method = "reml")
    } else{
      mv <- mvmeta(ymat, Slist, method='ml')
    }
    
    # Meta-analysis w/ prediction basis
    bound <- c(colMeans(ranges))
    mean <- seq(bound[1],bound[2],length=30)
    b_mean <- crossbasis(mean, lag=lag_window, argvar=arg_var, arglag=arg_lag)
    
    return(list(ymat=ymat, Slist=Slist, arg_var=arg_var, arg_lag=arg_lag, ranges=ranges, mv=mv, 
                bound=bound, b_mean=b_mean, cen=cen))
    
  }
  
  fam_type='binom'; groups=unique(df$groups); longterm_df=8; autocorr_lag=0; lag_window=c(0,10); arglag_tmp=3; arglag_tp=4; argvar_tmp=3; argvar_tp=1
  
  # 1. Data and dependency prep ----
  
  # LIST OF DATAFRAMES FOR ALL GROUPS
  datalist <- lapply(groups, function(g) df[df$groups == g, ])
  names(datalist) <- groups
  
  # Drop empty datasets
  datalist <- Filter(function(x) nrow(x) > 0, datalist)
  
  # Use names that survived filtering
  grps_used <- names(datalist)
  m <- length(datalist)
  years <- length(unique(df$year))
  
  # TMP
  ranges_tmp <- t(sapply(datalist,function(x) range(x$tmp,na.rm=T)))
  cen_tmp <- sapply(datalist,function(x) {
    quantile(x$tmp,c(0:100)/100,na.rm=T)
  })
  arglag_tmp <- list(fun='ns', df=arglag_tmp)
  argvar_tmp <- if (argvar_tmp == 1) list(fun="lin") else list(fun="ns", df=argvar_tmp)
  df_tmp <- prod(arglag_tmp$df, argvar_tmp$df)
  
  # TP
  ranges_tp <- t(sapply(datalist,function(x) range(x$tp,na.rm=T)))
  cen_tp <- sapply(datalist,function(x) {
    quantile(x$tp,c(0:100)/100,na.rm=T)
  })
  arglag_tp <- list(fun='ns', df=arglag_tp)
  argvar_tp <- if (argvar_tp == 1) list(fun="lin") else list(fun="ns", df=argvar_tp)
  df_tp <- prod(arglag_tp$df, argvar_tp$df)
  
  # 2. Stage 1: Fit the DLNM ----
  
  # BUILT OBJECTS WHERE outputs WILL BE STORED:
  
  Slist_comp <- vector('list',m)
  Slist_model <- vector('list', m)
  names(Slist_comp) <- grps_used
  names(Slist_model) <- grps_used
  
  ymat_tmp <- matrix(NA,m,df_tmp,dimnames=list(grps_used,paste('spl',seq(df_tmp),sep='')))
  Slist_tmp <- vector('list',m)
  names(Slist_tmp) <- grps_used
  
  ymat_tp <- matrix(NA,m,df_tp,dimnames=list(grps_used,paste('spl',seq(df_tp),sep='')))
  Slist_tp <- vector('list',m)
  names(Slist_tp) <- grps_used
  
  # mmi_tmp <- data.frame(location = locations, mmi_tmp = NA_real_, mmi_tmp_pctl = NA_character_)
  
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  options(warn=-1)
  
  system.time(
    for(i in seq(m)) {
      
      # PRINT ITERATION
      cat(i,'')
      
      # LOAD
      sub <- datalist[[i]]
      
      # SELECT AUTOCORR LAGS
      if(autocorr_lag==0){
        lags <- ''
      } else {
        lag_terms <- paste0(outcome, '_lag', 1:autocorr_lag)
        lags <- paste0('+ ', paste(lag_terms, collapse=' + '))
      }
      
      # CREATE THE SPLINES - CUMULATIVE
      cb.tmp <- crossbasis(sub$tmp, lag=lag_window, argvar=argvar_tmp, arglag=arglag_tmp)
      cb.tp <- crossbasis(sub$tp, lag=lag_window, argvar=argvar_tp, arglag=arglag_tp)
      controls <- '+ weekday + holiday + rh + fg10 + pm25' 
      
      # Set formula
      form <- as.formula(
        paste0(outcome,
               " ~ cb.tmp + cb.tp + ns(time, df=", longterm_df, "*", years, ')',
               lags, controls)
      )
      
      # RUN THE MODEL
      if(fam_type=='quasi'){
        model <- glm(formula=form, family=quasipoisson(), data=sub)
      } else if (fam_type=='binom'){
        model <- glm.nb(formula=form,data=sub)
      } else {
        stop('Unknown model type: ', fam_type)
      }
      
      # EXTRACT AND SAVE THE RELATED COEF AND VCOV
      
      # TMP
      predtmp <- crosspred(cb.tmp, model, cen=median(sub$tmp, na.rm=TRUE))
      ymat_tmp[i,] <- predtmp$coef
      Slist_tmp[[i]] <- predtmp$vcov
      # in_range <- predtmp$predvar >= cen_tmp['25%', i] & predtmp$predvar <= cen_tmp['75%', i]
      # mmi <- predtmp$predvar[in_range][which.min(predtmp$allRRfit[in_range])]
      # predtmp <- crosspred(cb.tmp, model, cen=mmi)
      # mmi_tmp$mmi_tmp[i] <- mmi
      # mmi_tmp$mmi_tmp_pctl[i] <- rownames(cen_tmp)[which.min(abs(cen_tmp[, i] - mmi))]
      
      # TP
      predtp <- crosspred(cb.tp, model, cen=median(sub$tp, na.rm=TRUE))
      ymat_tp[i,] <- predtp$coef
      Slist_tp[[i]] <- predtp$vcov
      
      # EXTRACT COMPARATIVE METRICS
      Slist_comp[[i]] <- comp_metrics(model)
      Slist_model[[i]] <- model
      Slist_comp[[i]]$esti_tmp <- unname(predtmp$coefficients[1])
      Slist_comp[[i]]$esti_tp <- unname(predtp$coefficients[1])
      
    })
  
  # 3. Stage 2: Multivariate Meta-analysis ----
  
  # Meta metrics
  meta_comp <- list()
  meta_comp[['avg_tmp']] <- mean(sapply(Slist_comp, function(x) x$esti_tmp))
  meta_comp[['avg_tp']] <- mean(sapply(Slist_comp, function(x) x$esti_tp))
  meta_comp[['sum_aic']] <- sum(sapply(Slist_comp, function(x) x$aic))
  meta_comp[['sum_pacf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$pacf)))
  meta_comp[['sum_acf']] <- rowMeans(sapply(Slist_comp, function(x) as.numeric(x$acf)))
  
  # Meta coefs
  meta_tmp <- meta_analysis(groups=groups, ymat=ymat_tmp, Slist=Slist_tmp, arg_var=argvar_tmp, arg_lag=arglag_tmp, ranges=ranges_tmp, cen=cen_tmp, lag_window=lag_window)
  meta_tp <- meta_analysis(groups=groups, ymat=ymat_tp, Slist=Slist_tp, arg_var=argvar_tp, arg_lag=arglag_tp, ranges=ranges_tp, cen=cen_tp, lag_window=lag_window)
  
  # RESET WARNING
  options(warn=0)
  
  return(list(meta_comp=meta_comp, meta_tmp=meta_tmp, meta_tp=meta_tp)) #models=Slist_model, 
  
}

# Function to extract predictions
meta_predict <- function(meta, city='ALL', reference='50%', blup=FALSE){
  
  # Unlist meta-analysis w/ prediction basis
  list2env(meta, envir=environment())
  
  cen  <- cen[!apply(cen, 1, function(row) any(duplicated(row))), ] # remove duplicates
  perc <- as.numeric(sub("%", "", names(rowMeans(cen))))
  reference <- if(reference=='0%') paste0(min(perc), '%') else reference # set reference
  
  if (city=='ALL') {
    cp <- crosspred(b_mean, coef=coef(mv), vcov=vcov(mv), model.link='log', 
                    at=rowMeans(cen), cen=rowMeans(cen)[reference])
  } else {
    i <- which(rownames(ymat) == city)
    if(blup==TRUE){
      blup <- blup(mv, vcov=TRUE)
      coef=blup[[i]]$blup; vcov=blup[[i]]$vcov
    } else{
      coef=ymat[i, ]; vcov=Slist[[i]]
    }
    
    cp <- crosspred(b_mean, coef=coef, vcov=vcov, model.link='log', from=ranges[i, 1], to=ranges[i, 2], 
                    at=cen[, i], cen=cen[reference, i])
  }
  
  return(list(cp=cp, perc=perc))
  
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
          scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
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
  
  # Make all predictions (pooled, city-raw, city-blup)
  all       <- meta_predict(meta, city = 'ALL', reference = reference)
  city_raw  <- map(locations, ~ meta_predict(meta, city = .x, reference = reference)) %>% setNames(locations)
  city_blup <- map(locations, ~ meta_predict(meta, city = .x, reference = reference, blup = TRUE)) %>% setNames(locations)
  
  # Get appropriate data frame
  if (!is.null(lag)) {
    df_overall <- data.frame(pred = all$perc, RR = all$cp$matRRfit[, lag + 1],low = all$cp$matRRlow[, lag + 1], high = all$cp$matRRhigh[, lag + 1], row.names = NULL)
    df_city_raw <- map_dfr(city_raw, function(x) {data.frame(pred = x$perc, RR = x$cp$matRRfit[, lag + 1],low = x$cp$matRRlow[, lag + 1], high = x$cp$matRRhigh[, lag + 1], row.names = NULL)}, .id = "city")
    df_city_blup <- map_dfr(city_blup, function(x) {data.frame(pred = x$perc, RR = x$cp$matRRfit[, lag + 1],low = x$cp$matRRlow[, lag + 1], high = x$cp$matRRhigh[, lag + 1], row.names = NULL)}, .id = "city")
  } else {
    df_overall <- data.frame(pred = all$perc, RR = all$cp$allRRfit,low = all$cp$allRRlow, high = all$cp$allRRhigh, row.names = NULL)
    df_city_raw <- map_dfr(city_raw, function(x) {data.frame(pred = x$perc, RR = x$cp$allRRfit, low = x$cp$allRRlow, high = x$cp$allRRhigh, row.names = NULL)}, .id = "city")
    df_city_blup <- map_dfr(city_blup, function(x) {data.frame(pred = x$perc, RR = x$cp$allRRfit, low = x$cp$allRRlow, high = x$cp$allRRhigh, row.names = NULL)}, .id = "city")
  }
  
  # Plot pooled and city specific (raw + blup) estimates
  if(city=='ALL'){
    
    p_raw  <- plot_rr_curve(city=city, df_specific1=df_city_raw, df_overall=df_overall, xlab=xlab, col=col, title=title, legend='Hospital-specific (Raw)')
    p_blup <- plot_rr_curve(city=city, df_specific1=df_city_blup, df_overall=df_overall, xlab=xlab, col=col, title=title, legend='Hospital-specific (BLUP)')
    
  } else{
    title <- if (!is.null(lag)) paste0('Lag=', lag) else city
    df_city_raw <- df_city_raw[df_city_raw$city==city,]
    df_city_blup <- df_city_blup[df_city_blup$city==city,]
    p_raw  <- plot_rr_curve(city=city, df_specific1=df_city_raw, df_specific2=df_city_blup, df_overall=df_overall, xlab=xlab, col=col, title=title)
    p_blup <- plot_rr_curve(city=city, df_specific1=df_city_blup, df_overall=df_overall, xlab=xlab, col=col, title=title, legend='Hospital-specific (BLUP)')
    
  }
  
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
          data = subset(df, sig_dir == 'neg'),
          aes(x = pred, y = lag), pattern = "crosshatch", pattern_fill = "#12669e",
          pattern_colour = '#12669e', fill = NA, color = NA, 
          pattern_density = 0.05, pattern_spacing = 0.05,
          pattern_angle = 45, pattern_alpha = 0.2, pattern_res = 200,
          width = 1, height = 1
        ) +
        geom_tile_pattern(
          data = subset(df, sig_dir == 'pos'),
          aes(x = pred, y = lag), pattern = "crosshatch", pattern_fill = "#a0392e",
          pattern_colour = "#a0392e", fill = NA, color = NA, 
          pattern_density = 0.05, pattern_spacing = 0.05,
          pattern_angle = 45, pattern_alpha = 0.2, pattern_res = 200,
          width = 1, height = 1 
        )
    }
    p
  }
  
  get_df <- function(cp, perc) {
    df <- expand.grid(pred=perc, lag=seq(cp$lag[1], cp$lag[2]))
    df$RR      <- as.vector(cp$matRRfit)
    df$low  <- as.vector(cp$matRRlow)
    df$high <- as.vector(cp$matRRhigh)
    df$sig_dir <- with(df, ifelse((low > 1 | high < 1) & abs(RR - 1) > 0.01, ifelse(RR > 1, "pos", "neg"), "none"))
    df
  }
  
  
  
  if (city == 'ALL') {
    pred <- meta_predict(meta, city=city, reference=reference)
    plot_df <- get_df(pred$cp, pred$perc)
    p <- plot_contour(df=plot_df, xlab=xlab, title=title, sig=sig, shared_lims=shared_lims, shared_brks=shared_brks)
  } else {
    pred <- meta_predict(meta, city=city, reference=reference, blup=TRUE)
    plot_df <- get_df(pred$cp, pred$perc)
    p <- plot_contour(df=plot_df, xlab=xlab, title=city, sig=sig, shared_lims=shared_lims, shared_brks=shared_brks)
  }
  
  return(list(df=plot_df, p=p))
  
}




