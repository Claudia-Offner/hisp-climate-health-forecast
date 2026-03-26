### ------------------------------------------------------------------------ ### 
### Ensemble - Dependencies
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation and set seed
options(scipen=999)
set.seed(123)


#### PATHS ####
user <- 'claer14'
path_data <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/', user, '/Documents/GitHub/hisp-climate-health-forecast/')
path_project <- paste0(path_git, '3_Ensemble/')
path_report <- paste0(path_git, '0_Reporting/')
path_results <- paste0(path_project, "Results/")
path_figures <- paste0(path_report, "III. ensemble/")

# FUNCTIONS: Paireau - Real time smoother ####
# ── Basis matrix: polynomial design matrix of degree d over window [-h, h] ──
Xhd = function(h, d) {
  res = rep(1, 2*h+1)
  for (i in 1:d) { res = cbind(res, seq(-h, h, by=1)**i) }
  as.matrix(res)
}

# ── Kernel weights: bi-square kernel over window, downweighting observations
#    far from centre. Returns diagonal weight matrix of size (2h+1) ──
K = function(h) { diag((15/16 * (1-seq(-h/(h+1), h/(h+1), 1/(h+1))^2)^2)) / (h+1) }


# ── Precompute smoother weight matrices for all three window positions ──
# h     : bandwidth (half-window size)
# d     : polynomial degree
# dZ    : degree used for asymmetric correction (defaults to d)
# d.sig : tuning parameter controlling strength of asymmetric correction
# Returns:
#   W        — weight matrix for a full symmetric window (interior points)
#   W.before — list of h weight matrices for the start of the series
#              (window is right-truncated: only future observations available)
#   W.after  — list of h weight matrices for the end of the series
#              (window is left-truncated: only past observations available)
comp.weights <- function(h, d, dZ, d.sig) {
  if (is.null(dZ)) dZ = d
  dU = 1:(dZ-1)
  require(MASS)
  
  K = K(h)
  X = Xhd(h, d)
  W = ginv(t(X) %*% K %*% X) %*% t(X) %*% K   # standard LPR weights
  
  W.before <- list()
  for (k in 1:h) {
    # k steps before enough data exists; shrink window by dropping first p rows
    p  = h - k + 1
    Xp = X[-(1:p), ];  Kp = K[-(1:p), -(1:p)]
    Dp = diag(rep(1, 2*h+1-p))
    Zp = Xp[, dZ, drop=F];  Up = Xp[, dU, drop=F]
    # Ridge-like penalty on the dZ-th polynomial component (asymmetric correction)
    Qp = Dp + (d.sig)^2 * Zp %*% t(Zp)
    Qp.inv = ginv(Qp)
    Up.Qpinv.Up = ginv(t(Up) %*% Qp.inv %*% Up)
    Lp = Qp.inv %*% Up %*% Up.Qpinv.Up
    Mp = Qp.inv - Lp %*% t(Up) %*% Qp.inv
    # Combine truncated weights with a correction from the missing observations
    Vp = t(W[1:3, -(1:p)]) + (Lp %*% t(X[(1:p), dU, drop=F]) +
                                (d.sig)^2 * Mp %*% Zp %*% t(X[(1:p), dZ, drop=F])) %*% t(W[1:3, (1:p)])
    W.before[[k]] = Vp
  }
  
  W.after <- list()
  for (k in 1:h) {
    # k steps into the tail; shrink window by dropping last rows
    p  = 2*h + 1 - k
    Xp = X[1:p, ];  Kp = K[1:p, 1:p]
    Dp = diag(rep(1, p))
    Zp = Xp[, dZ, drop=F];  Up = Xp[, dU, drop=F]
    Qp = Dp + (d.sig)^2 * Zp %*% t(Zp)
    Qp.inv = ginv(Qp)
    Up.Qpinv.Up = ginv(t(Up) %*% Qp.inv %*% Up)
    Lp = Qp.inv %*% Up %*% Up.Qpinv.Up
    Mp = Qp.inv - Lp %*% t(Up) %*% Qp.inv
    Vp = t(W[1:3, 1:p]) + (Lp %*% t(X[-(1:p), dU, drop=F]) +
                             (d.sig)^2 * Mp %*% Zp %*% t(X[-(1:p), dZ, drop=F])) %*% t(W[1:3, -(1:p)])
    W.after[[k]] = Vp
  }
  
  return(list(W=W, W.before=W.before, W.after=W.after))
}


# ── Core Local Polynomial Regression (LPR) smoother ──
# Fits a degree-d polynomial in a rolling kernel-weighted window of half-width h.
# At the boundaries (first/last h points) uses asymmetric weights to avoid
# requiring future data (real-time safe).
#
# data      : data frame with a column named by log.incid
# h         : bandwidth — half the window width in days
# d         : polynomial degree for the fit
# dZ        : degree used for the asymmetric boundary correction (default = d)
# d.sig     : strength of asymmetric correction (0 = none, larger = stronger)
# log.incid : name of the column containing log(incidence)
#
# Returns a list of smoothed trajectories:
#   y.h — standard symmetric LPR (requires future data, not real-time safe)
#   y.a — truncated window LPR at boundaries (naive boundary handling)
#   y.t — leave-one-out CV residuals (used to select bandwidth)
#   y.v — asymmetrically corrected smoother (real-time safe)
#   y.c — final corrected smoother (= y.v; used downstream)
#   CV.bandwidth  — CV score for bandwidth selection
#   CV.asymetric  — CV score for d.sig selection
smoothLPR <- function(data, h, d, dZ=NULL, d.sig=0, log.incid="log.incid") {
  y = data[, log.incid]
  if (is.null(dZ)) dZ = d
  dU = 1:(dZ-1)
  require(MASS)
  
  res.y  = matrix(rep(NA, 3*length(y)), nrow=3)  # symmetric LPR
  res.ya = matrix(rep(NA, 3*length(y)), nrow=3)  # truncated window
  res.yt = matrix(rep(NA, 3*length(y)), nrow=3)  # LOO-CV values
  res.yv = matrix(rep(NA, 3*length(y)), nrow=3)  # asymmetric correction
  res.yc = matrix(rep(NA, 3*length(y)), nrow=3)  # final output
  
  K = K(h);  X = Xhd(h, d)
  W = ginv(t(X) %*% K %*% X) %*% t(X) %*% K
  
  # Asymmetric weight vector V for the real-time boundary (current edge of data)
  D = diag(rep(1, h+1))
  Z = X[1:(h+1), dZ, drop=F];  U = X[1:(h+1), dU, drop=F]
  Q = D + (d.sig)^2 * Z %*% t(Z)
  Q.inv = ginv(Q)
  U.Qinv.U = ginv(t(U) %*% Q.inv %*% U)
  L = Q.inv %*% U %*% U.Qinv.U
  M = Q.inv - L %*% t(U) %*% Q.inv
  V = W[1, 1:(h+1)] + (L %*% t(X[-(1:(h+1)), dU, drop=F]) +
                         (d.sig)^2 * M %*% Z %*% t(X[-(1:(h+1)), dZ, drop=F])) %*% W[1, -(1:(h+1))]
  
  # Start of series: only k + h observations available (k = 1..h)
  for (k in 1:h) {
    p  = h - k + 1
    Xp = X[-(1:p), ];  Kp = K[-(1:p), -(1:p)]
    Wp = ginv(t(Xp) %*% Kp %*% Xp) %*% t(Xp) %*% Kp
    res.ya[, k] = (Wp %*% y[seq(1, h+k, 1)])[1:3]    # naive truncated fit
    Dp = diag(rep(1, 2*h+1-p))
    Zp = Xp[, dZ, drop=F];  Up = Xp[, dU, drop=F]
    Qp = Dp + (d.sig)^2 * Zp %*% t(Zp)
    Qp.inv = ginv(Qp)
    Up.Qpinv.Up = ginv(t(Up) %*% Qp.inv %*% Up)
    Lp = Qp.inv %*% Up %*% Up.Qpinv.Up
    Mp = Qp.inv - Lp %*% t(Up) %*% Qp.inv
    Vp = t(W[1:3, -(1:p)]) + (Lp %*% t(X[(1:p), dU, drop=F]) +
                                (d.sig)^2 * Mp %*% Zp %*% t(X[(1:p), dZ, drop=F])) %*% t(W[1:3, (1:p)])
    res.yv[, k] = t(y[seq(1, k+h, 1)]) %*% Vp      # corrected fit
    res.yc[, k] = t(y[seq(1, k+h, 1)]) %*% Vp
  }
  
  # Interior: full symmetric window available
  for (k in (h+1):(length(y)-h)) {
    res.y[, k]    = (W %*% y[(k-h):(k+h)])[1:3]
    res.ya[, k]   = (W %*% y[(k-h):(k+h)])[1:3]
    res.yt[, k]   = 1/(1-W[1,h+1]) * res.y[1,k] - W[1,h+1]/(1-W[1,h+1]) * y[k]  # LOO-CV
    res.yv[1, k]  = t(y[(k-h):(k)]) %*% V          # real-time estimate (no future data)
    res.yc[, k]   = (W %*% y[(k-h):(k+h)])[1:3]
  }
  
  # End of series: symmetric window would require future data
  for (k in (length(y)-h+1):length(y)) {
    p  = 2*h + 1 - (k - length(y) + h)
    Xp = X[1:p, ];  Kp = K[1:p, 1:p]
    Wp = ginv(t(Xp) %*% Kp %*% Xp) %*% t(Xp) %*% Kp
    res.ya[, k] = (Wp %*% y[seq(k-h, length(y), 1)])[1:3]
    Dp = diag(rep(1, p))
    Zp = Xp[, dZ, drop=F];  Up = Xp[, dU, drop=F]
    Qp = Dp + (d.sig)^2 * Zp %*% t(Zp)
    Qp.inv = ginv(Qp)
    Up.Qpinv.Up = ginv(t(Up) %*% Qp.inv %*% Up)
    Lp = Qp.inv %*% Up %*% Up.Qpinv.Up
    Mp = Qp.inv - Lp %*% t(Up) %*% Qp.inv
    Vp = t(W[1:3, 1:p]) + (Lp %*% t(X[-(1:p), dU, drop=F]) +
                             (d.sig)^2 * Mp %*% Zp %*% t(X[-(1:p), dZ, drop=F])) %*% t(W[1:3, -(1:p)])
    res.yv[, k] = t(y[seq(k-h, length(y), 1)]) %*% Vp
    res.yc[, k] = t(y[seq(k-h, length(y), 1)]) %*% Vp
  }
  
  list(y.h=res.y, y.a=res.ya, y.t=res.yt, y.v=res.yv, y.c=res.yc,
       CV.bandwidth  = sum((y - res.yt[1,])^2 / sum(!is.na(y - res.yt[1,])), na.rm=T),
       CV.asymetric  = sum((res.yv[1,] - res.yt[1,])^2, na.rm=T) /
         sum(!is.na(res.yv[1,] - res.yt[1,]), na.rm=T))
}


# ── Estimate and remove day-of-week reporting effects ──
# Computes a baseline smooth, finds the mean residual per weekday,
# then subtracts those weekday offsets from the log-incidence series.
# Returns the corrected data and the weekday effect table.
remove.day.effect <- function(data, bwd, d, dZ, d.sig, log.incid) {
  data$wd   = strftime(data$day, "%a")   # weekday label (Mon, Tue, ...)
  data$wd.n = strftime(data$day, "%u")   # weekday number (1 = Mon)
  res.wd    = smoothLPR(data, h=bwd, d=d, dZ=dZ, d.sig=0.0, log.incid=log.incid)
  data$delta.days = data[, log.incid] - res.wd$y.h[1, ]   # residual from smooth
  wd.effect = aggregate(data$delta.days, list(wd=data$wd, wd.n=data$wd.n), mean, na.rm=T)
  names(wd.effect) = c("wd", "day.num", "day.effect")
  wd.effect = wd.effect[order(wd.effect$day.num), ]
  wd.effect$wd = factor(wd.effect$wd, levels=wd.effect$wd)
  data = merge(data, wd.effect, by="wd")
  data = data[order(data$day), ]
  data$log.incid.cor = data[, log.incid] - data$day.effect   # corrected series
  return(list(data=data, wd.effect=wd.effect))
}


# ── Top-level wrapper: real-time smoother with automatic tuning ──
# Selects the best bandwidth (h) and asymmetric correction strength (d.sig)
# by cross-validation, then returns the final smoothed series.
#
# data     : data frame with columns named by day and log.incid
# day      : name of the date column
# log.incid: name of the log(incidence) column
# n.weeks  : how many weeks of recent data to use
# d / dZ   : polynomial degree / asymmetric degree (passed to smoothLPR)
# weekday  : if TRUE, remove day-of-week effects before smoothing
# bwd      : bandwidth used for the weekday-effect baseline smooth
# grid.d.sig / grid.bw : grids searched during CV (sensible defaults provided)
#
# Returns: smoothed series (y), full smoother object, CV tables,
#          weekday effects, raw data, and selected BW / DS parameters.
real.time.smoother <- function(data, day="day", log.incid="log.incid",
                               n.weeks=6, d=3, dZ=3, weekday=F, bwd=7,
                               grid.d.sig=NULL, grid.bw=NULL) {
  if (is.null(grid.d.sig)) { grid.d.sig = sort(c(0, 10^(-seq(1,6)), 2.5*10^(-seq(1,6)),
                                                 5*10^(-seq(1,6)), 7.5*10^(-seq(1,6)))) }
  if (is.null(grid.bw))    { grid.bw = seq(5, 15) }
  
  data = data[data$day > max(data$day) - n.weeks*7, ]   # restrict to recent n.weeks
  data = data[order(data$day), ]
  if (length(data$day) < n.weeks*7) warning("number of weeks is less than requested")
  
  if (weekday) {
    tmp.data   = remove.day.effect(data, bwd, d, dZ, d.sig=0, log.incid)
    wd.effect  = tmp.data$wd.effect
    data       = tmp.data$data
  } else {
    wd.effect          = NULL
    data$log.incid.cor = data[, log.incid]
  }
  
  # CV over bandwidth grid
  res.CV.bandwidth = sapply(grid.bw, function(h)
    smoothLPR(data, h=h, d=d, dZ=dZ, d.sig=0.1, log.incid="log.incid.cor")$CV.bandwidth)
  CV.bandwidth = data.frame(h=grid.bw, CV=res.CV.bandwidth)
  best.BW = CV.bandwidth$h[which.min(CV.bandwidth$CV)]
  
  # CV over d.sig grid (using best bandwidth)
  res.CV.asymetric = sapply(grid.d.sig, function(ds)
    smoothLPR(data, h=best.BW, d=d, dZ=dZ, d.sig=ds, log.incid="log.incid.cor")$CV.asymetric)
  CV.asymetric = data.frame(d.sig=grid.d.sig, CV=res.CV.asymetric)
  best.DS = CV.asymetric$d.sig[which.min(CV.asymetric$CV)]
  
  res.fin = smoothLPR(data, h=best.BW, d=d, dZ=dZ, d.sig=best.DS, log.incid="log.incid.cor")
  
  return(list(y=res.fin$yc[1,], smooth=res.fin,
              CV.asymetric=CV.asymetric, CV.bandwidth=CV.bandwidth,
              wd.effect=wd.effect, data=data, BW=best.BW, DS=best.DS))
}


# ── Bootstrap confidence intervals for the smoothed series ──
# Resamples residuals (smooth.data$data$log.incid.cor - fitted) R times,
# re-applies the smoother, and returns 2.5 / 50 / 97.5 percentiles
# (or min/max if minmax = TRUE) for the level, first, and second derivative.
#
# smooth.data : output from real.time.smoother()
# d / dZ      : must match values used in real.time.smoother()
# R           : number of bootstrap replicates
confint.LPR <- function(smooth.data, d, dZ, y="log.incid", R=1000, minmax=F) {
  n  = length(smooth.data$data$day)
  DS = smooth.data$DS;  h = smooth.data$BW
  W  = comp.weights(h=h, d=d, dZ=dZ, d.sig=DS)
  res.all = NULL
  
  for (i in 1:R) {
    # Residual bootstrap: resample residuals and add back to fitted values
    y.bs   = smooth.data$smooth$y.c[1,] +
      sample(smooth.data$data$log.incid.cor - smooth.data$smooth$y.c[1,])
    res.yc = matrix(rep(NA, 3*n), nrow=3)
    for (k in 1:h)           res.yc[, k]  = t(y.bs[seq(1, k+h, 1)]) %*% W$W.before[[k]]
    for (k in (h+1):(n-h))   res.yc[, k]  = (W$W %*% y.bs[(k-h):(k+h)])[1:3]
    for (k in (n-h+1):n)     res.yc[, k]  = t(y.bs[seq(k-h, n, 1)]) %*% W$W.after[[k-(n-h)]]
    res.all = rbind(res.all, cbind(rep(i, n), t(res.yc)))
  }
  
  res.all = data.frame(res.all)
  names(res.all) = c("i", "y", "y.prime", "y.second")
  res.all$day = smooth.data$data$day
  
  quant.fn = if (minmax) list(lo=min, hi=max) else
    list(lo=function(x) quantile(x, 0.025, na.rm=T),
         hi=function(x) quantile(x, 0.975, na.rm=T))
  
  aa.lo  = aggregate(res.all[, c("y","y.prime","y.second")], list(day=res.all$day), quant.fn$lo)
  aa.up  = aggregate(res.all[, c("y","y.prime","y.second")], list(day=res.all$day), quant.fn$hi)
  aa.med = aggregate(res.all[, c("y","y.prime","y.second")], list(day=res.all$day),
                     function(x) quantile(x, 0.5, na.rm=T))
  names(aa.lo)  = c("day","y.lo","y.prime.lo","y.second.lo")
  names(aa.up)  = c("day","y.up","y.prime.up","y.second.up")
  names(aa.med) = c("day","y.med","y.prime.med","y.second.med")
  
  aa = merge(merge(aa.lo, aa.med, by="day"), aa.up, by="day")
  return(aa)
}


# ── Return all R bootstrap trajectories (rather than just quantiles) ──
# Useful for plotting individual resampled paths or custom uncertainty summaries.
# Arguments and bootstrap procedure identical to confint.LPR().
boot.traj.LPR <- function(smooth.data, d, dZ, y="log.incid", R=1000) {
  n  = length(smooth.data$data$day)
  DS = smooth.data$DS;  h = smooth.data$BW
  W  = comp.weights(h=h, d=d, dZ=dZ, d.sig=DS)
  res.all = NULL
  
  for (i in 1:R) {
    y.bs   = smooth.data$smooth$y.c[1,] +
      sample(smooth.data$data$log.incid.cor - smooth.data$smooth$y.c[1,])
    res.yc = matrix(rep(NA, 3*n), nrow=3)
    for (k in 1:h)           res.yc[, k] = t(y.bs[seq(1, k+h, 1)]) %*% W$W.before[[k]]
    for (k in (h+1):(n-h))   res.yc[, k] = (W$W %*% y.bs[(k-h):(k+h)])[1:3]
    for (k in (n-h+1):n)     res.yc[, k] = t(y.bs[seq(k-h, n, 1)]) %*% W$W.after[[k-(n-h)]]
    res.all = rbind(res.all, cbind(rep(i, n), t(res.yc)))
  }
  
  res.all = data.frame(res.all)
  names(res.all) = c("i", "y", "y.prime", "y.second")
  res.all$day = smooth.data$data$day
  return(res.all)
}


# FUNCTIONS: Paireau - Simulate beds utilities ####

require(Rcpp)
require(dfoptim)
require(forecast)
sourceCpp("3_Ensemble/postprocess.cpp")

# Build the function used to estimate the mean and the coefficient of
#  variation for the duration of stay in hospital
build_fHosp <- function(data, sims, pICUs, tauHosp_ICU, start_fit_at,
                        fixed_cv = -1) {
  n_tot <- nrow(data)
  if (start_fit_at <= min(data$date)) {
    ind <- 1
  } else {
    ind <- which(data$date == start_fit_at)
  }
  if (fixed_cv < 0) { # Coefficient of variation must be estimated
    fHosp <- function(x) {
      sim_inHosp <- sim_beds_Hosp(sims$iHosp, pICUs, tauHosp_ICU, x[1], x[2],
                                  FALSE)
      mean((sim_inHosp[ind:n_tot] - data$inHosp[ind:n_tot])^2)
    }
  } else { # Fixed coefficient of variation
    fHosp <- function(x) {
      sim_inHosp <- sim_beds_Hosp(sims$iHosp, pICUs, tauHosp_ICU, x[1],
                                  fixed_cv, FALSE)
      mean((sim_inHosp[ind:n_tot] - data$inHosp[ind:n_tot])^2)
    }
  }
  return(fHosp)
}

# Build the function used to estimate the mean and the coefficient of
#  variation for the duration of stay in ICU
build_fICU <- function(data, sims, start_fit_at, fixed_cv = -1) {
  n_tot <- nrow(data)
  if (start_fit_at <= min(data$date)) {
    ind <- 1
  } else {
    ind <- which(data$date == start_fit_at)
  }
  if (fixed_cv < 0) { # Coefficient of variation must be estimated
    fICU <- function(x) {
      inICU_sim <- sim_beds_ICU(sims$iICU, x[1], x[2], FALSE)
      mean((inICU_sim[ind:n_tot] - data$inICU[ind:n_tot])^2)
    }
  } else { # Fixed coefficient of variation
    fICU <- function(x) {
      inICU_sim <- sim_beds_ICU(sims$iICU, x[1], fixed_cv, FALSE)
      mean((inICU_sim[ind:n_tot] - data$inICU[ind:n_tot])^2)
    }
  }
  return(fICU)
}

# Compute stay in hosp and ICU (mean and coefficient of variation)
compute_durations <- function(data, sims, pICUs, tauHosp_ICU, start_fit_at,
                              fixed_cv_hosp = -1, fixed_cv_ICU = -1) {
  if (nrow(data) != nrow(sims)) {
    stop("Data and simulations do not have the same length!")
  }
  if (nrow(data) != length(pICUs)) {
    stop("Data and pICUs do not have the same length!")
  }
  if (start_fit_at > max(data$date)) {
    stop("Fit must start within data window!")
  }
  # TODO: check that data and sims have all required columns
  
  # Optimize fit for hospital beds
  fHosp <- build_fHosp(data, sims, pICUs, tauHosp_ICU, start_fit_at,
                       fixed_cv_hosp)
  if (fixed_cv_hosp < 0) { # Estimate coefficient of variation
    fit_Hosp <- nmkb(c(10., 1.), fHosp, lower = c(0.5, 1e-2), upper = c(60, 3))
    mu_tauHosp_R <- fit_Hosp$par[1]
    cv_tauHosp_R <- fit_Hosp$par[2]
  } else {
    fit_Hosp <- optim(c(10.), fHosp, method = "Brent",
                      lower = c(0.5), upper = c(60))
    mu_tauHosp_R <- fit_Hosp$par[1]
    cv_tauHosp_R <- fixed_cv_hosp
  }
  
  # Optimize fit for ICU beds
  fICU <- build_fICU(data, sims, start_fit_at, fixed_cv_ICU)
  if (fixed_cv_ICU < 0) { # Estimate coefficient of variation
    fit_ICU <- nmkb(c(10., 1.), fICU, lower = c(0.5, 1e-2), upper = c(60, 3))
    mu_tauICU_R <- fit_ICU$par[1]
    cv_tauICU_R <- fit_ICU$par[2]
  } else {
    fit_ICU <- optim(c(10.), fICU, method = "Brent",
                     lower = c(0.5), upper = c(60))
    mu_tauICU_R <- fit_ICU$par[1]
    cv_tauICU_R <- fixed_cv_ICU
  }
  
  return(list(mu_tauHosp_R = mu_tauHosp_R,
              cv_tauHosp_R = cv_tauHosp_R,
              mu_tauICU_R = mu_tauICU_R,
              cv_tauICU_R = cv_tauICU_R))
}


# Return the estimated probability of ICU given hosp and the estimated durations
#  of stay in hosp and ICU
calibrateV2 <- function(obs, model_sims, start_fit_at, fixed_cv_hosp = -1,
                        fixed_cv_ICU = -1, tauHosp_ICU = 1.5,
                        n_to_forecast = 1, project_pICU = "const") {
  # NOTE
  # 1) dataframe 'obs' must contain all data since March 1st
  # 2) start_fit_at = start date for calibration
  # 3) tauHosp_ICU = Average delay hosp -> ICU admission
  if (!(all(c("date", "iICU") %in% colnames(obs)))) {
    stop("Dataframe containing data must include 'date' and 'iICU' columns!")
  }
  if (!(all(c("date", "iHosp") %in% colnames(model_sims)))) {
    stop("Model dataframe must include 'date' and 'iHosp' columns!")
  }
  if (any(obs$date != model_sims$date)) {
    stop("Data and model have different dates!")
  }
  if (start_fit_at > max(obs$date)) {
    stop("Fit must start within data window!")
  }
  if (n_to_forecast < 1) {
    stop("'n_to_forecast' must be a positive number!")
  }
  if (!(project_pICU %in% c("const", "ets"))) {
    warning("'project_pICU' can be either 'const' or 'ets', using 'const'!")
    project_pICU <- "const"
  }
  
  # Pull relevant columns
  obs_iICU <- obs$iICU
  sims_iHosp <- model_sims$iHosp
  
  # Estimate pICU (vector)
  window <- 10
  raw_pICUs <- estimate_pICUs(sims_iHosp, obs_iICU, tauHosp_ICU, window)
  # Remove negative numbers (replaced with average)
  raw_pICUs[raw_pICUs < 0] <- mean(raw_pICUs[raw_pICUs > 0])
  fit_pICUs <- ets(raw_pICUs)
  pICUs <- pmin(pmax(fit_pICUs$fitted, 0), 1) # Enforce in [0, 1]
  
  # Reconstruct ICU admissions
  estim_iICU <- simulate_ICU_admissions(sims_iHosp, pICUs, tauHosp_ICU)
  
  trajectories <- model_sims %>%
    mutate(iICU = estim_iICU)
  
  # Estimate durations of stay in hosp and ICU
  estimates <- compute_durations(obs, trajectories, pICUs, tauHosp_ICU,
                                 start_fit_at, fixed_cv_hosp, fixed_cv_ICU)
  
  # Project future pICU
  if (project_pICU == "const") {
    last_pICU <- pICUs[length(pICUs)]
    pICUs <- c(pICUs, rep(last_pICU, n_to_forecast))
  } else {
    proj_pICUs <- forecast(fit_pICUs, h = n_to_forecast)
    pICUs <- c(pICUs, pmin(pmax(proj_pICUs$mean, 0), 1)) # Enforce in [0, 1]
  }
  # Add pICUs to list to return
  estimates$pICUs <- pICUs
  
  return(estimates)
}

# Project ICU admissions, hosp and ICU beds
project <- function(model_sims, pICUs, mu_tauHosp_R, cv_tauHosp_R,
                    mu_tauICU_R, cv_tauICU_R, tauHosp_ICU = 1.5) {
  # NOTE
  # 1) tauHosp_ICU = Average delay hosp -> ICU admission
  if (!(all(c("date", "iHosp") %in% colnames(model_sims)))) {
    stop("Model dataframe must include 'date' and 'iHosp' columns!")
  }
  if (nrow(model_sims) != length(pICUs)) {
    stop("Model dataframe and vector pICU must have same length!")
  }
  if (any(pICUs < 0.) | any(pICUs > 1.)) {
    stop("pICU must be in [0, 1]!")
  }
  
  estim_iICU <- simulate_ICU_admissions(model_sims$iHosp, pICUs, tauHosp_ICU)
  
  trajectories <- model_sims %>%
    mutate(iICU = estim_iICU,
           inHosp = sim_beds_Hosp(iHosp, pICUs, tauHosp_ICU,
                                  mu_tauHosp_R, cv_tauHosp_R, FALSE),
           inICU = sim_beds_ICU(iICU, mu_tauICU_R, cv_tauICU_R, FALSE))
  
  return(trajectories)
}

# FUNCTIONS: Paireau - Change tools ####

poly.smooth <- function(y, deg=1, deg0=deg, K=NULL, tau=NULL, beta=NULL, x=NULL, w=1, gamma=0) {
  ina <- which(!is.na(y))
  n <- length(y)
  if (is.null(x))
    x <- 1:n
  x <- x[ina]
  n <- length(x)
  
  if (is.null(K))
    K <- round(n/10)
  
  if (is.null(tau))
    tau <- round(seq((3*x[1] - x[2])/2 , (3*x[n] - x[n-1])/2, length.out = (K+1)) - 0.5) + 0.5
  else {
    if (tau[1] > x[1])
      tau <- c((3*x[1] - x[2])/2 , tau)
    if (tau[length(tau)] < x[n])
      tau <- c(tau, (3*x[n] - x[n-1])/2)
  }
  K <- length(tau)-1
  
  X <-   cbind(1,x)
  if (deg>1) 
    for (j in 2:deg)
      X <- cbind(X, x^j)
  
  if (K>1) {
    for (k in (2:K)) {
      xk <- pmax(0, x-tau[k])
      for (m in deg0:deg)
        X <- cbind(X, xk^m)
      #    X <- cbind(X, xk^deg)
    }
  }
  
  if (is.null(beta)) {
    Xw <- X*w
    A <- diag(c(rep(0, deg+1), rep(gamma, K-1)))
    beta <- solve(t(Xw)%*%X + A)%*%t(Xw)%*%y[ina]
  } 
  yp <- y
  yp[ina] <- X%*%beta
  return(list(tau=tau, yp=yp, x=x, X=X, beta=c(beta)))
}

estim.tau <- function(y, Kmax=12, Kmin=1, lmin=7, gamma=0) {
  Kj <- Kmax
  r.opt <- vector(mode="list", length=Kmax)
  r.opt[[Kmax]] <- estim.tau1(y, K=Kmax, lmin=lmin, gamma=gamma)
  for (Kj in ((Kmax-1):Kmin)) {
    tau <- r.opt[[Kj+1]]$tau
    r.min <- Inf
    for (k in (2:(Kj+1))) {
      rk <- estim.tau1(y, tau=tau[-k], lmin=lmin, gamma=gamma)
      if (rk$value<r.min) {
        r.min <- rk$value
        r.opt[[Kj]] <- rk
      }
    }
  }
  U <- NULL
  T <- matrix(NA, nrow=Kmax+1, ncol=(Kmax-Kmin+1))
  for (k in Kmin:Kmax) {
    U <- c(U, r.opt[[k]]$value)
    T[1:(k+1), k-Kmin+1] <- r.opt[[k]]$tau
  }
  U <- data.frame(K=Kmin:Kmax, U=U)
  T <- as.data.frame(T)
  names(T) <- Kmin:Kmax
  return(list(U=U, r.opt=r.opt, tau=T))
} 

estim.tau1 <- function(y, K=NULL, tau=NULL, lmin=NULL, gamma=0) {
  y0 <- y
  r.ini <- poly.smooth(y, K=K, tau=tau, gamma=gamma)
  yp <- r.ini$yp
  v.min <- sum((y-yp)^2)
  tau <- r.ini$tau
  J <- length(tau) 
  K <- J-1
  X <- r.ini$X
  x <- r.ini$x
  A <- diag(c(0, 0, rep(gamma, K-1)))
  if (J>2) {
    nc <- 2
    while (nc>0) {
      nc <- 0
      for (j in seq(2,(J-1))) {
        ind.j <- j+1
        tau.c <- tau
        tau.c[j] <- tau[j]-1
        X.c <- X
        X.c[,ind.j] <-  pmax(0, x-tau.c[j])
        yp <- X.c%*%solve(t(X.c)%*%X.c + A)%*%t(X.c)%*%y
        v.c <- sum((y-yp)^2)
        #       if (is.na(v.c))  browser()
        dm <- F
        while (v.c < v.min & tau.c[j]-tau.c[j-1]>=lmin) {
          nc <- nc+1
          dm <- T
          v.min <- v.c
          tau <- tau.c
          X <- X.c
          tau.c[j] <- tau[j]-1
          X.c[,ind.j] <-  pmax(0, x-tau.c[j])
          yp <- X.c%*%solve(t(X.c)%*%X.c + A)%*%t(X.c)%*%y
          v.c <- sum((y-yp)^2)
        }
        if (dm==F) {
          tau.c <- tau
          tau.c[j] <- tau[j]+1
          X.c[,ind.j] <-  pmax(0, x-tau.c[j])
          yp <- X.c%*%solve(t(X.c)%*%X.c + A)%*%t(X.c)%*%y
          v.c <- sum((y-yp)^2)
          dp <- F
          while (v.c < v.min  & tau.c[j+1]-tau.c[j]>=lmin) {
            nc <- nc+1
            dp <- T
            v.min <- v.c
            tau <- tau.c
            X <- X.c
            tau.c[j] <- tau[j]+1
            X.c[,ind.j] <-  pmax(0, x-tau.c[j])
            yp <- X.c%*%solve(t(X.c)%*%X.c + A)%*%t(X.c)%*%y
            v.c <- sum((y-yp)^2)
          }
        }
      } 
    }
  }
  r.opt <- poly.smooth(y, tau=tau, gamma=gamma)
  return(list(tau=tau, yp=r.opt$yp, beta=r.opt$beta, value=v.min))
}

interval <- function(r=NULL, y=NULL, n.pred=30, n.prev=7, alpha=0.05, gamma=0, deg=1) {
  tau <- r$tau
  d <- data.frame(x=1:length(y), y=y, yp=r$yp)
  
  s <- sd(d$y-d$yp)
  p0 <- poly.smooth(y, tau=tau, gamma=gamma, deg=deg)
  X <- p0$X
  
  K <- length(tau) -1
  A <- diag(c(rep(0,deg+1), rep(gamma, K-1)))
  var.beta <- solve(t(X)%*%X + A)%*%t(X)%*%X%*%solve(t(X)%*%X + A)*s^2
  
  y1 <- c(rep(0, n.prev), y, rep(0, n.pred))
  x1 <- (-n.prev+1):(length(y)+n.pred)
  tau1 <- tau
  tau1[1] <- tau1[1] - n.prev
  tau1[length(tau)] <- tau1[length(tau)] + n.pred
  X1 <- poly.smooth(y1, tau=tau1, x=x1, gamma=gamma, deg=deg)$X
  
  beta <- r$beta
  
  qa <- qt(c(alpha, 1-alpha), length(y)-ncol(X))
  d1 <- data.frame(x = x1)
  
  yp1 <- X1%*%beta
  var.f <- diag(X1%*%var.beta%*%t(X1))
  sigma.f1 <- sqrt(var.f)
  sigma.f2 <- sqrt(var.f + s^2)
  conf.dist <- data.frame(fit=yp1, lwr.conf=yp1+qa[1]*sigma.f1, upr.conf=yp1+qa[2]*sigma.f1)
  pred.dist <- data.frame(fit=yp1, lwr.conf=yp1+qa[1]*sigma.f2, upr.conf=yp1+qa[2]*sigma.f2)
  d1[c("fit","lwr.conf", "upr.conf")] <- conf.dist
  d1[c("lwr.pred", "upr.pred")] <- pred.dist[,2:3]
  
  return(list(data=d1, var.beta=var.beta))
}

interval2 <- function(r=NULL, y=NULL, n.pred=30, lag=0, n.prev=0, alpha=0.05, gamma=0, deg=1) {
  tau <- r$tau
  d <- data.frame(x=1:length(y), y=y, yp=r$yp)
  
  s.par <- optim(c(var(d$y-d$yp), 0), sigma.ll, d=d)$par
  s2 <- (s.par[1] + s.par[2]*d$yp)^1
  p0 <- poly.smooth(y, tau=tau, gamma=gamma, deg=deg)
  X <- p0$X
  
  K <- length(tau) -1
  A <- diag(c(rep(0,deg+1), rep(gamma, K-1)))
  var.beta <- solve(t(X/s2)%*%X + A)%*%t(X)%*%X%*%solve(t(X)%*%X + A)
  
  y1 <- c(rep(0, n.prev), y, rep(0, n.pred-lag))
  x1 <- (-n.prev+1):(length(y)+n.pred-lag)
  tau1 <- tau
  tau1[1] <- tau1[1] - n.prev
  tau1[length(tau)] <- tau1[length(tau)] + n.pred - lag
  X1 <- poly.smooth(y1, tau=tau1, x=x1, gamma=gamma, deg=deg)$X
  
  beta <- r$beta
  
  n=length(y)-lag
  a1 <- rep(1*alpha, nrow(X1))
  ia <- (n-1):(n+n.pred)
  a1[ia] <- 1.6*a1[ia]*(0.07^seq(0,0.5*(length(ia)-1),by=0.5))
  
  qa1 <- qt(a1/2, length(y)-ncol(X))
  qa2 <- qt(1-a1/2, length(y)-ncol(X))
  
  d1 <- data.frame(x = x1)
  
  yp1 <- X1%*%beta
  var.f <- diag(X1%*%var.beta%*%t(X1))
  var.f[n:(n+n.pred)] <- 2*var.f[n:(n+n.pred)]
  sigma.f1 <- sqrt(var.f)
  sigma.f2 <- sqrt(var.f + pmax(0,s.par[1] + s.par[2]*yp1))
  conf.dist <- data.frame(fit=yp1, lwr.conf=yp1+qa1*sigma.f1, upr.conf=yp1+qa2*sigma.f1)
  pred.dist <- data.frame(fit=yp1, lwr.conf=yp1+qa1*sigma.f2, upr.conf=yp1+qa2*sigma.f2)
  d1[c("fit","lwr.conf", "upr.conf")] <- conf.dist
  d1[c("lwr.pred", "upr.pred")] <- pred.dist[,2:3]
  
  return(list(data=d1, var.beta=var.beta))
}

seq.detect <- function(y, S=0.0075, lmin=7) {
  n <- length(y)
  k <- 0
  tau.last <- 0.5
  tau.est <- tau.detect <- res <- NULL
  i2 <- 0
  while (i2 < n) {
    i1 <- tau.last + 0.5
    i2 <- i1 + 2*lmin -1
    test <- F
    k <- k+1
    while (!test & i2 < n) {
      i2 <- i2 + 1
      yj <- y[i1:i2]
      nj <- length(yj)
      X0 <- cbind(1, 1:nj)
      yp0 <- X0%*%solve(t(X0)%*%X0)%*%t(X0)%*%yj
      e0 <- sum((yj-yp0)^2)
      es <- NULL
      vj <- (lmin+0.5):(nj-lmin+0.5)
      for (s in vj) {
        Xs <- cbind(X0, pmax((1:nj)-s, 0))
        yps <- Xs%*%solve(t(Xs)%*%Xs)%*%t(Xs)%*%yj
        es <- c(es, e0 - sum((yj-yps)^2))
      }
      res <- rbind(res, data.frame(n=i2, tau=vj[min(which.max(es))]+tau.last-0.5, e=max(es), k=as.factor(k)))
      if (max(es) > S) {
        tau.last <- vj[min(which.max(es))] + tau.last - 0.5
        tau.est <- c(tau.est, tau.last)
        tau.detect <- c(tau.detect, i2)
        test <- T
      }
    }
  }
  res.tau <- data.frame(tau.est=tau.est, tau.detect=tau.detect, k=as.factor(1:length(tau.est)))
  return(list(res.tau=res.tau, res.e=res))
}

pred.fhosp <- function(A, lag=8) {
  A <- A[which(!is.na(A$fhosp)),]
  lmod <- lm(log(fhosp)~log(fR), data=A)
  a <- lmod$coefficients[2]
  b <- lmod$coefficients[1]
  A <- A %>% mutate(f2=exp(a*log(fR) + b))
  n <- nrow(A)
  P <- rbind(A[,c("date","fhosp")], data.frame(date=A$date[n]+(1:lag), fhosp=A$fhosp[n]/A$f2[n-lag]*A$f2[(n-lag+1):n]))
  P <- P %>% filter(!is.na(fhosp))
  P <- P %>% mutate(khosp=c(NA, diff(fhosp))/c(NA,fhosp[-nrow(P)])) 
  P$khosp[!is.na(P$khosp) & P$khosp< -0.12] <- -0.12
  return(P)
}

sigma.ll <- function(p, d) {
  s2 <- (p[1] + p[2]*d$yp)^1
  if (min(s2) <=0) 
    obj <- Inf
  else
    obj <- sum(((d$y-d$yp)^2)/s2 + log(s2))
  return(obj)
}


#### DEPENDENCIES ####

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD', '#0D2C6E', '#910030FF', '#127804FF', '#52C5B0FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9', '#4462A2', '#91003080', '#12780480', '#52C5B080') 

nm <- c(total_patients='Total patients', M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        general='General medicine', emergency='Emergency', ophthalmology='Ophthalmology', 
        paediatrics='Paediatric', orthopedics='Orthopedics', ent='ENT', 
        gynecology='Gynecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        mx2t_IMD='IMD: Mx2t >30\u00B0C', mn2t_IMD='IMD: Mn2t <0\u00B0C', tp_IMD='IMD: Tp >15mm', fg10_IMD='IMD: Fg10 >40km/h',
        mx2t_p90='Clim: Mx2t >90th', mn2t_p10='Clim: Mn2t <10th', tp_p90='Clim: Tp >90th', fg10_p90='Clim: Fg10 >90th',
        mx2t='Mx2t (\u00B0C)', mn2t='Mn2t (\u00B0C)', tp='Tp (mm)', fg10='Fg10 (km/h)', 
        closed='Closed', holiday='Holiday') 

weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
h_cols <- c('total_patients', 'F', 'M', 'y1_l', 'y1_17', 'y18_60', 'y60_g', 
            'general_medicine', 'ophthalmology', 'pediatrics', 'orthopedics', 
            'ent', 'gynecology', 'dermatology', 'surgery', 'dentistry', 
            'emergency', 'psychiatry', 'oncology')





