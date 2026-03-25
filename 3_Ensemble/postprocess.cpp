#include <Rcpp.h>
using namespace Rcpp;


// -----------------------------------------------------------------------------
// Helper function to determine time period of input time step
// -----------------------------------------------------------------------------
int get_index(int it, IntegerVector spoints) {
  // No change points
  if (spoints.size() == 0) return 0;
  // Before first change point
  if (it < spoints[0]) return 0;
  // Intermediate change points
  for (int curr = 1; curr < spoints.size(); ++curr) {
    if ( (it >= spoints[curr - 1]) && (it < spoints[curr]) ) {
      return curr;
    }
  }
  // After last change point
  return spoints.size();
}


// [[Rcpp::export]]
NumericVector estimate_pICUs(NumericVector iHosp, NumericVector iICU,
                             double tauHosp_ICU, int window) {
  if (iHosp.size() != iICU.size()) {
    stop("ERROR: check length of iHosp and iICU!");
  }
  int n_steps = iHosp.size();
  
  // Proba hosp -> ICU
  NumericVector ws(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    double x1 = step == 0 ? 0.0 : step - 0.5;
    double x2 = step + 0.5;
    ws[step] = R::pexp(x2, tauHosp_ICU, true, false) -
      R::pexp(x1, tauHosp_ICU, true, false);
  }
  
  NumericVector the_sum(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    for (int curr = step; curr < n_steps; ++curr) {
      the_sum[curr] += iHosp[step] * ws[curr - step];
    }
  }
  
  NumericVector pICUs(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    double num = 0.0;
    double den = 0.0;
    for (int curr = 0; curr < window; ++curr) {
      if (step - curr < 0) break;
      num += iICU[step - curr];
      den += the_sum[step - curr];
    }
    if ((den > 0.0) && (den >= num)) {
      pICUs[step] = num / den;
    } else {
      pICUs[step] = -1.;
    }
  }
  
  return pICUs;
}


// [[Rcpp::export]]
double estimate_pICU(NumericVector iHosp, NumericVector iICU,
                     double tauHosp_ICU, int window) {
  NumericVector pICUs = estimate_pICUs(iHosp, iICU, tauHosp_ICU, window);
  return pICUs[pICUs.size() - 1];
}


// [[Rcpp::export]]
NumericVector simulate_ICU_admissions(NumericVector iHosp, NumericVector pICUs,
                                      double tauHosp_ICU) {
  if (iHosp.size() != pICUs.size()) {
    stop("ERROR: check length of iHosp and pICUs!");
  }
  int n_steps = iHosp.size();

  // Proba hosp -> ICU
  NumericVector ws(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    double x1 = step == 0 ? 0.0 : step - 0.5;
    double x2 = step + 0.5;
    ws[step] = R::pexp(x2, tauHosp_ICU, true, false) -
      R::pexp(x1, tauHosp_ICU, true, false);
  }

  NumericVector iICU(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    for (int curr = step; curr < n_steps; ++curr) {
      iICU[curr] += pICUs[step] * iHosp[step] * ws[curr - step];
    }
  }
  
  return iICU;
}


// -----------------------------------------------------------------------------
// Simulate beds: delay distribution = gamma parametrized by mean and cv
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector sim_beds_Hosp(NumericVector iHosp, NumericVector pICUs,
                            double tauHosp_ICU, double mu_tauHosp_R,
                            double cv_tauHosp_R, bool draw = true) {
  if (iHosp.size() != pICUs.size()) {
    stop("ERROR: check length of iHosp and pICU_Hosp!");
  }
  int n_steps = iHosp.size();
  NumericVector inHosp_ave(n_steps);
  IntegerVector inHosp(n_steps);
  
  // Survival function for hosp -> ICU
  NumericVector surv_ICU(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    double x = step == 0 ? 0.5 : step;
    surv_ICU[step] = 1. - R::pexp(x, tauHosp_ICU, true, false);
  }
  // Survival function for hosp -> recovery
  NumericMatrix surv_R(n_steps);
  double curr_shape = 1.0 / pow(cv_tauHosp_R, 2);
  double curr_scale = pow(cv_tauHosp_R, 2) * mu_tauHosp_R;
  for (int step = 0; step < n_steps; ++step) {
    double x = step == 0 ? 0.5 : step;
    surv_R[step] = 1. - R::pgamma(x, curr_shape, curr_scale, true, false);
  }
  
  for (int step = 0; step < n_steps; ++step) {
    double pICU = pICUs[step];
    for (int curr = step; curr < n_steps; ++curr) {
      inHosp_ave[curr] += pICU * iHosp[step] * surv_ICU[curr - step] +
        (1. - pICU) * iHosp[step] * surv_R[curr - step];
    }
    inHosp[step] = draw ? R::rpois(inHosp_ave[step]) : inHosp_ave[step];
  }
  
  return inHosp;
}


// [[Rcpp::export]]
IntegerVector sim_beds_ICU(NumericVector iICU, double mu_tauICU_R,
                           double cv_tauICU_R, bool draw = true) {
  int n_steps = iICU.size();
  NumericVector inICU_ave(n_steps);
  IntegerVector inICU(n_steps);
  
  // Survival function for ICU -> recovery
  NumericMatrix surv_R(n_steps);
  double curr_shape = 1.0 / pow(cv_tauICU_R, 2);
  double curr_scale = pow(cv_tauICU_R, 2) * mu_tauICU_R;
  for (int step = 0; step < n_steps; ++step) {
    double x = step == 0 ? 0.5 : step;
    surv_R[step] = 1. - R::pgamma(x, curr_shape, curr_scale, true, false);
  }
  
  for (int step = 0; step < n_steps; ++step) {
    for (int curr = step; curr < n_steps; ++curr) {
      inICU_ave[curr] += iICU[step] * surv_R[curr - step];
    }
    inICU[step] = draw ? R::rpois(inICU_ave[step]) : inICU_ave[step];
  }
  
  return inICU;
}


// [[Rcpp::export]]
IntegerVector sim_beds_Hosp_cps(NumericVector iHosp, NumericVector pICUs,
                                double tauHosp_ICU,
                                IntegerVector spoints_tauHosp_R,
                                NumericVector mu_tauHosp_R,
                                NumericVector cv_tauHosp_R,
                                bool draw = true) {
  if (mu_tauHosp_R.size() != (spoints_tauHosp_R.size() + 1)) {
    stop("ERROR: check change points for tauHosp_R (mu)!");
  }
  if (cv_tauHosp_R.size() != (spoints_tauHosp_R.size() + 1)) {
    stop("ERROR: check change points for tauHosp_R (cv)!");
  }
  if (iHosp.size() != pICUs.size()) {
    stop("ERROR: check length of iHosp and pICUs!");
  }
  int n_steps = iHosp.size();
  NumericVector inHosp_ave(n_steps);
  IntegerVector inHosp(n_steps);

  // Survival function for hosp -> ICU
  NumericVector surv_ICU(n_steps);
  for (int step = 0; step < n_steps; ++step) {
    double x = step == 0 ? 0.5 : step;
    surv_ICU[step] = 1. - R::pexp(x, tauHosp_ICU, true, false);
  }
  // Survival function for hosp -> recovery
  int n_tauHosp_R = mu_tauHosp_R.size();
  NumericMatrix surv_R(n_tauHosp_R, n_steps);
  for (int cp = 0; cp < n_tauHosp_R; ++cp) { // Time period
    double curr_shape = 1.0 / pow(cv_tauHosp_R[cp], 2);
    double curr_scale = pow(cv_tauHosp_R[cp], 2) * mu_tauHosp_R[cp];
    for (int step = 0; step < n_steps; ++step) {
      double x = step == 0 ? 0.5 : step;
      surv_R(cp, step) = 1. - R::pgamma(x, curr_shape, curr_scale, true, false);
    }
  }
  
  for (int step = 0; step < n_steps; ++step) {
    double pICU = pICUs[step];
    int cp = get_index(step, spoints_tauHosp_R); // Time period for tauHosp_R
    for (int curr = step; curr < n_steps; ++curr) {
      inHosp_ave[curr] += pICU * iHosp[step] * surv_ICU[curr - step] +
        (1. - pICU) * iHosp[step] * surv_R(cp, curr - step);
    }
    inHosp[step] = draw ? R::rpois(inHosp_ave[step]) : inHosp_ave[step];
  }
  
  return inHosp;
}


// [[Rcpp::export]]
IntegerVector sim_beds_ICU_cps(NumericVector iICU,
                               IntegerVector spoints_tauICU_R,
                               NumericVector mu_tauICU_R,
                               NumericVector cv_tauICU_R,
                               bool draw = true) {
  if (mu_tauICU_R.size() != (spoints_tauICU_R.size() + 1)) {
    stop("ERROR: check change points for tauICU_R (mu)!");
  }
  if (cv_tauICU_R.size() != (spoints_tauICU_R.size() + 1)) {
    stop("ERROR: check change points for tauICU_R (cv)!");
  }
  int n_steps = iICU.size();
  NumericVector inICU_ave(n_steps);
  IntegerVector inICU(n_steps);
  
  // Survival function for ICU -> recovery
  int n_tauICU_R = mu_tauICU_R.size();
  NumericMatrix surv_R(n_tauICU_R, n_steps);
  for (int cp = 0; cp < n_tauICU_R; ++cp) { // Time period
    double curr_shape = 1.0 / pow(cv_tauICU_R[cp], 2);
    double curr_scale = pow(cv_tauICU_R[cp], 2) * mu_tauICU_R[cp];
    for (int step = 0; step < n_steps; ++step) {
      double x = step == 0 ? 0.5 : step;
      surv_R(cp, step) = 1. - R::pgamma(x, curr_shape, curr_scale, true, false);
    }
  }
  
  for (int step = 0; step < n_steps; ++step) {
    int cp = get_index(step, spoints_tauICU_R); // Time period for tauHosp_R
    for (int curr = step; curr < n_steps; ++curr) {
      inICU_ave[curr] += iICU[step] * surv_R(cp, curr - step);
    }
    inICU[step] = draw ? R::rpois(inICU_ave[step]) : inICU_ave[step];
  }
  
  return inICU;
}

