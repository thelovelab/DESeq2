/*
 * DESeq2 C++ functions
 * 
 * Author: Michael I. Love
 * Last modified: September 30, 2015
 * License: LGPL (>= 3)
 *
 * Note: The canonical, up-to-date DESeq2.cpp lives in 
 * the DESeq2 library, the development branch of which 
 * can be viewed here: 
 *
 * https://github.com/Bioconductor-mirror/DESeq2/blob/master/src/DESeq2.cpp
 */

// include RcppArmadillo and Rcpp
#include "RcppArmadillo.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// user includes
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

// this function returns the log posterior of dispersion parameter alpha, for negative binomial variables
// given the counts y, the expected means mu, the design matrix x (used for calculating the Cox-Reid adjustment),
// and the parameters for the normal prior on log alpha
double log_posterior(double log_alpha, Rcpp::NumericMatrix::Row y, Rcpp::NumericMatrix::Row mu, arma::mat x, double log_alpha_prior_mean, double log_alpha_prior_sigmasq, bool use_prior) {
  double prior_part;
  double alpha = exp(log_alpha);
  Rcpp::NumericVector w_diag = pow(pow(mu, -1) + alpha, -1);
  arma::mat w = arma::diagmat(Rcpp::as<arma::vec>(w_diag));
  arma::mat b = x.t() * w * x;
  double cr_term = -0.5 * log(det(b));
  double alpha_neg1 = R_pow_di(alpha, -1);
  double ll_part = sum(lgamma(y + alpha_neg1) - Rf_lgammafn(alpha_neg1) - y * log(mu + alpha_neg1) - alpha_neg1 * log(1.0 + mu * alpha));
  if (use_prior) {
    prior_part = -0.5 * R_pow_di(log_alpha - log_alpha_prior_mean,2)/log_alpha_prior_sigmasq;
  } else {
    prior_part = 0.0;
  }
  double res =  ll_part + prior_part + cr_term;
  return(res);
}

// this function returns the derivative of the log posterior with respect to the log of the 
// dispersion parameter alpha, given the same inputs as the previous function
double dlog_posterior(double log_alpha, Rcpp::NumericMatrix::Row y, Rcpp::NumericMatrix::Row mu, arma::mat x, double log_alpha_prior_mean, double log_alpha_prior_sigmasq, bool use_prior) {
  double prior_part;
  double alpha = exp(log_alpha);
  Rcpp::NumericVector w_diag = pow(pow(mu, -1) + alpha, -1);
  arma::mat w = arma::diagmat(Rcpp::as<arma::vec>(w_diag));
  Rcpp::NumericVector dw_diag = -1.0 * pow(pow(mu, -1) + alpha, -2);
  arma::mat dw = arma::diagmat(Rcpp::as<arma::vec>(dw_diag));
  arma::mat b = x.t() * w * x;
  arma::mat db = x.t() * dw * x;
  double ddetb = ( det(b) * trace(b.i() * db) );
  double cr_term = -0.5 * ddetb / det(b);
  double alpha_neg1 = R_pow_di(alpha, -1);
  double alpha_neg2 = R_pow_di(alpha, -2);
  double ll_part = alpha_neg2 * sum(Rf_digamma(alpha_neg1) + log(1 + mu*alpha) - mu*alpha*pow(1.0 + mu*alpha, -1) - digamma(y + alpha_neg1) + y * pow(mu + alpha_neg1, -1));
  // only the prior part is w.r.t log alpha
  if (use_prior) {
    prior_part = -1.0 * (log_alpha - log_alpha_prior_mean)/log_alpha_prior_sigmasq;
  } else {
    prior_part = 0.0;
  }
  // Note: return dlog_post/dalpha * alpha because we take derivatives w.r.t log alpha
  double res = (ll_part + cr_term) * alpha + prior_part;
  return(res);
}

// this function returns the second derivative of the log posterior with respect to the log of the 
// dispersion parameter alpha, given the same inputs as the previous function
double d2log_posterior(double log_alpha, Rcpp::NumericMatrix::Row y, Rcpp::NumericMatrix::Row mu, arma::mat x, double log_alpha_prior_mean, double log_alpha_prior_sigmasq, bool use_prior) {
  double prior_part;
  double alpha = exp(log_alpha);
  Rcpp::NumericVector w_diag = pow(pow(mu, -1) + alpha, -1);
  arma::mat w = arma::diagmat(as<arma::vec>(w_diag));
  Rcpp::NumericVector dw_diag = -1 * pow(pow(mu, -1) + alpha, -2);
  arma::mat dw = arma::diagmat(as<arma::vec>(dw_diag));
  Rcpp::NumericVector d2w_diag = 2 * pow(pow(mu, -1) + alpha, -3);
  arma::mat d2w = arma::diagmat(as<arma::vec>(d2w_diag));
  arma::mat b = x.t() * w * x;
  arma::mat b_i = b.i();
  arma::mat db = x.t() * dw * x;
  arma::mat d2b = x.t() * d2w * x;
  double ddetb = ( det(b) * trace(b.i() * db) );
  double d2detb = ( det(b) * (R_pow_di(trace(b_i * db), 2) - trace(b_i * db * b_i * db) + trace(b_i * d2b)) );
  double cr_term = 0.5 * R_pow_di(ddetb/det(b), 2) - 0.5 * d2detb / det(b); 
  double alpha_neg1 = R_pow_di(alpha, -1);
  double alpha_neg2 = R_pow_di(alpha, -2);
  double ll_part = -2 * R_pow_di(alpha, -3) * sum(Rf_digamma(alpha_neg1) + log(1 + mu*alpha) - mu*alpha*pow(1 + mu*alpha, -1) - digamma(y + alpha_neg1) + y * pow(mu + alpha_neg1, -1)) + alpha_neg2 * sum(-1 * alpha_neg2 * Rf_trigamma(alpha_neg1) + pow(mu, 2) * alpha * pow(1 + mu*alpha, -2) + alpha_neg2 * trigamma(y + alpha_neg1) + alpha_neg2 * y * pow(mu + alpha_neg1, -2));
  // only the prior part is w.r.t log alpha
  if (use_prior) {
    prior_part = -1.0/log_alpha_prior_sigmasq; 
  } else {
    prior_part = 0.0;
  }
  // Note: return (d2log_post/dalpha2 * alpha^2 + dlog_post/dalpha * alpha) 
  //            = (d2log_post/dalpha2 * alpha^2 + dlog_post/dlogalpha)
  // because we take derivatives w.r.t log alpha
  double res = ((ll_part + cr_term) * R_pow_di(alpha, 2) + dlog_posterior(log_alpha, y, mu, x, log_alpha_prior_mean, log_alpha_prior_sigmasq, false)) + prior_part;
  return(res);
}

// Obtain the MLE or MAP dispersion estimate using line search.
// fitting occurs on the scale of log(alpha)
//
// [[Rcpp::export]]
Rcpp::List fitDisp(SEXP ySEXP, SEXP xSEXP, SEXP mu_hatSEXP, SEXP log_alphaSEXP, SEXP log_alpha_prior_meanSEXP, SEXP log_alpha_prior_sigmasqSEXP, SEXP min_log_alphaSEXP, SEXP kappa_0SEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP use_priorSEXP) {
  Rcpp::NumericMatrix y(ySEXP);
  arma::mat x = Rcpp::as<arma::mat>(xSEXP);
  int y_n = y.nrow();
  Rcpp::NumericVector log_alpha(clone(log_alphaSEXP));
  Rcpp::NumericMatrix mu_hat(mu_hatSEXP);
  Rcpp::NumericVector log_alpha_prior_mean(log_alpha_prior_meanSEXP);
  double log_alpha_prior_sigmasq = Rcpp::as<double>(log_alpha_prior_sigmasqSEXP);
  double min_log_alpha = Rcpp::as<double>(min_log_alphaSEXP);
  double kappa_0 = Rcpp::as<double>(kappa_0SEXP);
  int maxit = Rcpp::as<int>(maxitSEXP);
  double epsilon = 1.0e-4;
  double a, a_propose, kappa, lp, lpnew, dlp, theta_kappa, theta_hat_kappa, change;
  // record log posterior values
  Rcpp::NumericVector initial_lp(y_n);
  Rcpp::NumericVector initial_dlp(y_n);
  Rcpp::NumericVector last_lp(y_n);
  Rcpp::NumericVector last_dlp(y_n);
  Rcpp::NumericVector last_d2lp(y_n);
  Rcpp::NumericVector last_change(y_n);
  Rcpp::IntegerVector iter(y_n);
  Rcpp::IntegerVector iter_accept(y_n);
  double tol = Rcpp::as<double>(tolSEXP);
  bool use_prior = Rcpp::as<bool>(use_priorSEXP);

  for (int i = 0; i < y_n; i++) {
    Rcpp::checkUserInterrupt();
    Rcpp::NumericMatrix::Row yrow = y(i,_);
    Rcpp::NumericMatrix::Row mu_hat_row = mu_hat(i,_);
    // maximize the log likelihood over the variable a, the log of alpha, the dispersion parameter.
    // in order to express the optimization in a typical manner, 
    // for calculating theta(kappa) we multiple the log likelihood by -1 and seek a minimum
    a = log_alpha(i);
    // we use a line search based on the Armijo rule.
    // define a function theta(kappa) = f(a + kappa * d), where d is the search direction.
    // in this case the search direction is taken by the first derivative of the log likelihood
    lp = log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
    dlp = dlog_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
    kappa = kappa_0;
    initial_lp(i) = lp;
    initial_dlp(i) = dlp;
    change = -1.0;
    last_change(i) = -1.0;
    for (int t = 0; t < maxit; t++) {
      // iter counts the number of steps taken out of  maxit;
      iter(i)++;
      a_propose = a + kappa * dlp;
      // note: lgamma is unstable for values around 1e17, where there is a switch in lgamma.c
      // we limit log alpha from going lower than -30
      if (a_propose < -30.0) {
	kappa = (-30.0 - a)/dlp;
      }
      // note: we limit log alpha from going higher than 10
      if (a_propose > 10.0) {
	kappa = (10.0 - a)/dlp;
      }
      theta_kappa = -1.0 * log_posterior(a + kappa*dlp, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
      theta_hat_kappa = -1.0 * lp - kappa * epsilon * R_pow_di(dlp, 2);
      // if this inequality is true, we have satisfied the Armijo rule and 
      // accept the step size kappa, otherwise we halve kappa
      if (theta_kappa <= theta_hat_kappa) {
	// iter_accept counts the number of accepted proposals;
	iter_accept(i)++;
	a = a + kappa * dlp;
	lpnew = log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
	// look for change in log likelihood
	change = lpnew - lp;
	if (change < tol) {
	  lp = lpnew;
	  break;
	}
	// if log(alpha) is going to -infinity
	// break the loop
	if (a < min_log_alpha) {
	  break;
	}
	lp = lpnew;
	dlp = dlog_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
	// instead of resetting kappa to kappa_0 
	// multiple kappa by 1.1
	kappa = fmin(kappa * 1.1, kappa_0);
	// every 5 accepts, halve kappa
	// to prevent slow convergence
	// due to overshooting
	if (iter_accept(i) % 5 == 0) {
	  kappa = kappa / 2.0;
	}
      } else {
	kappa = kappa / 2.0;
      }
    }
    last_lp(i) = lp;
    last_dlp(i) = dlp;
    last_d2lp(i) = d2log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
    log_alpha(i) = a;
    // last change indicates the change for the final iteration
    last_change(i) = change;
  }

  return Rcpp::List::create(Rcpp::Named("log_alpha",log_alpha),
			    Rcpp::Named("iter",iter),
			    Rcpp::Named("iter_accept",iter_accept),
			    Rcpp::Named("last_change",last_change),
			    Rcpp::Named("initial_lp",initial_lp),
			    Rcpp::Named("initial_dlp",initial_dlp),
			    Rcpp::Named("last_lp",last_lp),
			    Rcpp::Named("last_dlp",last_dlp),
			    Rcpp::Named("last_d2lp",last_d2lp));
}

// fit the Negative Binomial GLM.
// note: the betas are on the natural log scale
//
// [[Rcpp::export]]
Rcpp::List fitBeta(SEXP ySEXP, SEXP xSEXP, SEXP nfSEXP, SEXP alpha_hatSEXP, SEXP contrastSEXP, SEXP beta_matSEXP, SEXP lambdaSEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP useQRSEXP) {
  arma::mat y = Rcpp::as<arma::mat>(ySEXP);
  arma::mat nf = Rcpp::as<arma::mat>(nfSEXP);
  arma::mat x = Rcpp::as<arma::mat>(xSEXP);
  int y_n = y.n_rows;
  int y_m = y.n_cols;
  int x_p = x.n_cols;
  arma::vec alpha_hat = Rcpp::as<arma::vec>(alpha_hatSEXP);
  arma::mat beta_mat = Rcpp::as<arma::mat>(beta_matSEXP);
  arma::mat beta_var_mat = arma::zeros(beta_mat.n_rows, beta_mat.n_cols);
  arma::mat contrast_num = arma::zeros(beta_mat.n_rows, 1);
  arma::mat contrast_denom = arma::zeros(beta_mat.n_rows, 1);
  arma::mat hat_matrix = arma::zeros(x.n_rows, x.n_rows);
  arma::mat hat_diagonals = arma::zeros(y.n_rows, y.n_cols);
  arma::colvec lambda = Rcpp::as<arma::colvec>(lambdaSEXP);
  arma::colvec contrast = Rcpp::as<arma::colvec>(contrastSEXP);
  int maxit = Rcpp::as<int>(maxitSEXP);
  arma::colvec yrow, nfrow, beta_hat, mu_hat, z;
  arma::mat w, ridge, sigma;
  // vars for QR
  bool useQR = Rcpp::as<bool>(useQRSEXP);
  arma::colvec gamma_hat, big_z;
  arma::vec big_w_diag;
  arma::mat weighted_x_ridge, q, r, big_w;
  // deviance, convergence and tolerance
  double dev, dev_old, conv_test;
  double tol = Rcpp::as<double>(tolSEXP);
  double large = 30.0;
  Rcpp::NumericVector iter(y_n);
  Rcpp::NumericVector deviance(y_n);
  // bound the estimated count, as weights include 1/mu
  double minmu = 0.5;
  for (int i = 0; i < y_n; i++) {
    Rcpp::checkUserInterrupt();
    nfrow = nf.row(i).t();
    yrow = y.row(i).t();
    beta_hat = beta_mat.row(i).t();
    mu_hat = nfrow % exp(x * beta_hat);
    for (int j = 0; j < y_m; j++) {
      mu_hat(j) = fmax(mu_hat(j), minmu);
    }
    ridge = diagmat(lambda);
    dev = 0.0;
    dev_old = 0.0;
    if (useQR) {
      // make an orthonormal design matrix including
      // the ridge penalty
      for (int t = 0; t < maxit; t++) {
	iter(i)++;
	w = diagmat(mu_hat/(1.0 + alpha_hat[i] * mu_hat));
	// prepare matrices
	weighted_x_ridge = join_cols(sqrt(w) * x, sqrt(ridge));
	qr(q, r, weighted_x_ridge);
	big_w_diag = arma::ones(y_m + x_p);
	big_w_diag(arma::span(0, y_m - 1)) = diagvec(w);
	big_w = diagmat(big_w_diag);
	big_z = arma::zeros(y_m + x_p);
	z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
	big_z(arma::span(0,y_m - 1)) = z;
	// IRLS with Q matrix for X    
	gamma_hat = q.t() * sqrt(big_w) * big_z;
	solve(beta_hat, r, gamma_hat);
	if (sum(abs(beta_hat) > large) > 0) {
	  iter(i) = maxit;
	  break;
	}
	mu_hat = nfrow % exp(x * beta_hat);
	for (int j = 0; j < y_m; j++) {
	  mu_hat(j) = fmax(mu_hat(j), minmu);
	}
	dev = 0.0;
	for (int j = 0; j < y_m; j++) {
	  // note the order for Rf_dnbinom_mu: x, sz, mu, lg
	  dev = dev + -2.0 * Rf_dnbinom_mu(yrow[j], 1.0/alpha_hat[i], mu_hat[j], 1);
	}
	conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
	if (std::isnan(conv_test)) {
	  iter(i) = maxit;
	  break;
	}
	if ((t > 0) & (conv_test < tol)) {
	  break;
	}
	dev_old = dev;
      }
    } else {
      // use the standard design matrix x
      // and matrix inversion
      for (int t = 0; t < maxit; t++) {
	iter(i)++;
	w = diagmat(mu_hat/(1.0 + alpha_hat[i] * mu_hat));
	z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
	solve(beta_hat, x.t() * w * x + ridge, x.t() * w * z);
	if (sum(abs(beta_hat) > large) > 0) {
	  iter(i) = maxit;
	  break;
	}
	mu_hat = nfrow % exp(x * beta_hat);
	for (int j = 0; j < y_m; j++) {
	  mu_hat(j) = fmax(mu_hat(j), minmu);
	}
	dev = 0.0;
	for (int j = 0; j < y_m; j++) {
	  // note the order for Rf_dnbinom_mu: x, sz, mu, lg
	  dev = dev + -2.0 * Rf_dnbinom_mu(yrow[j], 1.0/alpha_hat[i], mu_hat[j], 1);
	}
	conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
	if (std::isnan(conv_test)) {
	  iter(i) = maxit;
	  break;
	}
	if ((t > 0) & (conv_test < tol)) {
	  break;
	}
	dev_old = dev;
      }
    }
    deviance(i) = dev;
    beta_mat.row(i) = beta_hat.t();
    // recalculate w so that this is identical if we start with beta_hat
    w = diagmat(mu_hat/(1.0 + alpha_hat[i] * mu_hat));
    hat_matrix = sqrt(w) * x * (x.t() * w * x + ridge).i() * x.t() * sqrt(w);
    hat_diagonals.row(i) = diagvec(hat_matrix).t();
    // sigma is the covariance matrix for the betas
    sigma = (x.t() * w * x + ridge).i() * x.t() * w * x * (x.t() * w * x + ridge).i();
    contrast_num.row(i) = contrast.t() * beta_hat;
    contrast_denom.row(i) = sqrt(contrast.t() * sigma * contrast);
    beta_var_mat.row(i) = diagvec(sigma).t();
  }

  return Rcpp::List::create(Rcpp::Named("beta_mat",beta_mat),
			    Rcpp::Named("beta_var_mat",beta_var_mat),
			    Rcpp::Named("iter",iter),
			    Rcpp::Named("hat_diagonals",hat_diagonals),
			    Rcpp::Named("contrast_num",contrast_num),
			    Rcpp::Named("contrast_denom",contrast_denom),
			    Rcpp::Named("deviance",deviance));
}


// [[Rcpp::export]]
Rcpp::List fitDispGrid(SEXP ySEXP, SEXP xSEXP, SEXP mu_hatSEXP, SEXP disp_gridSEXP, SEXP log_alpha_prior_meanSEXP, SEXP log_alpha_prior_sigmasqSEXP, SEXP use_priorSEXP) {
  Rcpp::NumericMatrix y(ySEXP);
  arma::mat x = Rcpp::as<arma::mat>(xSEXP);
  int y_n = y.nrow();
  Rcpp::NumericMatrix mu_hat(mu_hatSEXP);
  arma::vec disp_grid = Rcpp::as<arma::vec>(disp_gridSEXP);
  int disp_grid_n = disp_grid.n_elem;
  Rcpp::NumericVector log_alpha_prior_mean(log_alpha_prior_meanSEXP);
  double log_alpha_prior_sigmasq = Rcpp::as<double>(log_alpha_prior_sigmasqSEXP);
  bool use_prior = Rcpp::as<bool>(use_priorSEXP);
  double a;
  double delta = disp_grid(1) - disp_grid(0);
  double a_hat;
  arma::vec disp_grid_fine;
  arma::vec logpostvec = arma::zeros(disp_grid.n_elem);
  arma::vec log_alpha = arma::zeros(y_n);
  arma::uword idxmax;
  for (int i = 0; i < y_n; i++) {
    Rcpp::checkUserInterrupt();
    Rcpp::NumericMatrix::Row yrow = y(i,_);
    Rcpp::NumericMatrix::Row mu_hat_row = mu_hat(i,_);
    for (int t = 0; t < disp_grid_n; t++) {
      // maximize the log likelihood over the variable a, the log of alpha, the dispersion parameter
      a = disp_grid(t);
      logpostvec(t) = log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
    }
    logpostvec.max(idxmax);
    a_hat = disp_grid(idxmax);
    disp_grid_fine = arma::linspace<arma::vec>(a_hat - delta, a_hat + delta, disp_grid_n);
    for (int t = 0; t < disp_grid_n; t++) {
      a = disp_grid_fine(t);
      logpostvec(t) = log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, use_prior);
    }
    logpostvec.max(idxmax);
    log_alpha(i) = disp_grid_fine(idxmax);
  }

  return Rcpp::List::create(Rcpp::Named("log_alpha",log_alpha));
}

