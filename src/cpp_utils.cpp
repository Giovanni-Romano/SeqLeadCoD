// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <iostream>
#include <gsl/gsl_sf_hyperg.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::uvec compatibility_check_cpp(int i, int H, const arma::uvec& gamma_tp1, 
                                   const arma::mat& c_t, const arma::mat& c_tp1) {
  // Step 1: If the i-th observation is free to move, compatibility is guaranteed
  if (gamma_tp1[i] == 0) {
    return arma::uvec(H + 1, fill::ones); // Return a vector of 1s of size H + 1
  }
  
  // Step 2: Find the constrained observations (R_tp1)
  arma::uvec R_tp1 = find(gamma_tp1 == 1); // Indices where gamma_tp1 == 1
  
  // Step 3: Remove i from R_tp1 to create R_tp1mi
  arma::uvec R_tp1mi = R_tp1.elem(find(R_tp1 != (unsigned)i)); // Exclude i
  
  // Rcpp::Rcout << "R_tp1mi size: " << R_tp1mi.n_elem << std::endl;
  
  // Step 4: If i is the only constrained observation, compatibility is guaranteed
  if (R_tp1mi.empty()) {
    return arma::uvec(H + 1, fill::ones); // Return a vector of 1s of size H + 1
  }
  
  // Step 5: Build the matrix for all possible allocations of the i-th observation
  arma::mat tmp1(H + 1, R_tp1mi.n_elem, fill::zeros);
  for (int row = 0; row < H; ++row) {
    tmp1(row, span(0, R_tp1mi.n_elem - 1)) = c_t.submat(R_tp1mi, arma::uvec{(unsigned)row}).t();
  }
  tmp1(H, span(0, R_tp1mi.n_elem - 1)) = arma::rowvec(R_tp1mi.n_elem, fill::zeros);
  
  // Step 6: Compute the i-th row of the co-clustering matrix at time t+1
  arma::rowvec tmp2 = c_tp1.row(i) * c_tp1.rows(R_tp1mi).t();
  
  // Step 7: Check compatibility for each possible allocation
  arma::uvec comp_checks(H + 1, fill::zeros);
  for (int row = 0; row < H + 1; ++row) {
    if (all(tmp1.row(row)== tmp2)) {
      comp_checks[row] = 1;
    }
  }
  
  return comp_checks;
}


// [[Rcpp::export]]
arma::vec dhamming_multicluster_cpp(const arma::vec& x, const arma::mat& mu, const arma::mat& sigma, 
                                    const arma::vec& m, bool logscale = false) {
  int H = mu.n_rows; // Number of rows
  
  // Step 1: Create matrix M by repeating vector m across rows
  arma::mat M = repmat(m.t(), H, 1); // Repeat m as a row vector across H rows
  
  // Step 2: Create matrix X by repeating vector x across rows
  arma::mat X = repmat(x.t(), H, 1); // Repeat x as a row vector across H rows
  
  // Step 3: Compute the inverse of sigma
  arma::mat inv_sigma = 1 / sigma;
  
  // Step 4: Compute the exponential of the inverse of sigma
  arma::mat exp_inv_sigma = exp(inv_sigma);
  
  // Step 5: Compute logdenom and lognum
  arma::mat logdenom = log1p((M - 1) / exp_inv_sigma); // log1p(x) = log(1 + x)
  arma::mat lognum = -((X != mu) % inv_sigma); // Element-wise comparison and multiplication
  
  // Step 6: Compute the row-wise sums of lognum and logdenom
  arma::vec sum_lognum = sum(lognum, 1); // Sum across columns (row-wise)
  arma::vec sum_logdenom = sum(logdenom, 1); // Sum across columns (row-wise)
  
  // Step 7: Compute the final output
  arma::vec logout = sum_lognum - sum_logdenom;
  
  // Step 8: If not in logscale, exponentiate the result
  if (!logscale) {
    logout = exp(logout);
  }
  
  return logout;
}



// [[Rcpp::export]]
double norm_const2(const double d ,const double c, const double m, const bool log_scale=false){
  
  double z= (m-1)/m;
  double alpha= d+c;
  double beta=1;
  double gamma= d+2;
  gsl_sf_result out;
  
  gsl_set_error_handler_off();
  
  int stat= gsl_sf_hyperg_2F1_e(alpha, beta, gamma, z, & out);
  
  //Rcpp::Rcout<<"Valore della hyper "<<out.val<<"\n";
  
  if (stat != GSL_SUCCESS)
  {
    //Rcpp::Rcout<<"Sono qui\n";
    return R_NaN;
  }
  
  
  if(log_scale){
    return -std::log(d+1)-(d+c)*std::log(m)+std::log(out.val);
    Rcpp::Rcout<<"Sono in log_scale\n";
  }
  else{
    return -std::exp(std::log(d+1)-(d+c)*std::log(m)+std::log(out.val));
    Rcpp::Rcout<<"Non sono in log_scale\n";
  }
  
} 


// [[Rcpp::export]]
arma::mat norm_const2_mat(const arma::ivec& clS,     // integer counts
                          const arma::imat& matches, // integer counts
                          const arma::vec& u,        // length p
                          const arma::vec& v,        // length p
                          const arma::ivec& m,        // length p
                          bool log_scale = true) {
  
  const int K = matches.n_rows, p = matches.n_cols;
  arma::mat out(K, p);
  
  
  for (int j = 0; j < p; ++j) {
    const double uj = u[j], vj = v[j];
    const int mj = m[j];
    for (int k = 0; k < K; ++k) {
      const int mat = matches(k,j);
      const int non = clS(k) - mat;
      const double c = uj + mat;
      const double d = vj + non;
      out(k,j) = norm_const2(d, c, mj, log_scale);
    }
  }
  return out;
}


// ---------- Build triangular tables for all s = 0..Smax ----------
// tabs is a List of length p; tabs[j] is a List of length Smax+1;
// tabs[j][s] is a NumericVector of length s+1 with values for x=0..s.
// [[Rcpp::export]]
Rcpp::List build_cached_norm_const(int Smax,
                                   const arma::vec& u,  // length p
                                   const arma::vec& v,  // length p
                                   const arma::vec& m)  // length p
{
  if (Smax < 0) stop("Smax must be >= 0.");
  const int p = u.n_elem;
  if ((int)v.n_elem != p || (int)m.n_elem != p) stop("u, v, m must have same length.");
  
  Rcpp::List tabs(p);
  for (int j = 0; j < p; ++j) {
    const double uj = u[j], vj = v[j], mj = m[j];
    Rcpp::List colTabs(Smax + 1);
    for (int s = 0; s <= Smax; ++s) {
      Rcpp::NumericVector tv(s + 1);
      for (int x = 0; x <= s; ++x) {
        const double d = vj + (s - x);
        const double c = uj + x;
        tv[x] = norm_const2(d, c, mj, /*log_scale=*/true);
      }
      colTabs[s] = tv;
    }
    tabs[j] = colTabs;
  }
  return tabs;
}


// ---------- O(1) evaluation using the tables ----------
// clS: K-vector (integers); matches: K x p (integers).
// [[Rcpp::export]]
arma::mat eval_norm_const(const arma::ivec& clS,
                          const arma::imat& matches,
                          const Rcpp::List& tabs)
{
  const int K = clS.n_elem;
  const int p = matches.n_cols;
  if (matches.n_rows != K) stop("matches must have K rows.");
  if (tabs.size() != p)    stop("tabs length must equal ncol(matches).");
  
  arma::mat out(K, p);
  for (int j = 0; j < p; ++j) {
    Rcpp::List colTabs = tabs[j];                 // length Smax+1
    for (int k = 0; k < K; ++k) {
      const int s = clS[k];
      if (s < 0 || s >= colTabs.size()) stop("clS[%d]=%d out of bounds.", k, s);
      Rcpp::NumericVector tv = colTabs[s];        // length s+1
      const int x = matches(k, j);
      if (x < 0 || x >= tv.size()) stop("matches(%d,%d)=%d out of bounds [0,%d].", k, j, x, tv.size()-1);
      out(k, j) = tv[x];
    }
  }
  return out;
}

// [[Rcpp::export]]
List count_matches_update_centers(const arma::imat& M, const arma::ivec& m) {
  int n = M.n_rows;
  int p = M.n_cols;
  
  if ((int)m.n_elem != p) stop("Length of m must equal number of columns of M.");
  
  List out(p);
  
  for (int j = 0; j < p; ++j) {
    int mj = m[j];
    IntegerVector counts(mj, 0);
    
    for (int i = 0; i < n; ++i) {
      int v = M(i, j);
      if (v == NA_INTEGER) continue;
      if (v >= 1 && v <= mj) counts[v - 1]++;
    }
    
    CharacterVector nm(mj);
    for (int k = 0; k < mj; ++k) nm[k] = std::to_string(k + 1);
    counts.attr("names") = nm;
    
    out[j] = counts;
  }
  
  return out;
}

// [[Rcpp::export]]
List logprob_update_centers(const arma::imat& M, const arma::ivec& m,
                            const int clS,
                            const arma::vec u, const arma::vec v) {
  int n = M.n_rows;
  int p = M.n_cols;
  
  if ((int)m.n_elem != p) stop("Length of m must equal number of columns of M.");
  
  List out(p);
  
  for (int j = 0; j < p; ++j) {
    int mj = m[j];
    IntegerVector counts(mj, 0);
    
    for (int i = 0; i < n; ++i) {
      int v = M(i, j);
      if (v == NA_INTEGER) continue;
      if (v >= 1 && v <= mj) counts[v - 1]++;
    }
    
    NumericVector logprobs(mj, 0.0);
    
    for (int l = 0; l < mj; ++l){
      logprobs[l] = norm_const2(v[j] + clS - counts[l], u[j] + counts[l],
                                mj, true);
    }
    
    CharacterVector nm(mj);
    for (int k = 0; k < mj; ++k) nm[k] = std::to_string(k + 1);
    logprobs.attr("names") = nm;
    
    out[j] = logprobs;
  }
  
  return out;
}


// [[Rcpp::export]]
List eval_logprob_update_centers(const arma::imat& M, const arma::ivec& m,
                            const int clS,
                            const Rcpp::List& tabs) {
  int n = M.n_rows;
  int p = M.n_cols;
  
  if ((int)m.n_elem != p) stop("Length of m must equal number of columns of M.");
  
  List out(p);
  
  for (int j = 0; j < p; ++j) {
    int mj = m[j];
    IntegerVector counts(mj, 0);
    
    for (int i = 0; i < n; ++i) {
      int v = M(i, j);
      if (v == NA_INTEGER) continue;
      if (v >= 1 && v <= mj) counts[v - 1]++;
    }
    
    NumericVector logprobs(mj, 0.0);
    Rcpp::List colTabs = tabs[j];
    Rcpp::NumericVector tv = colTabs[clS];
    
    for (int l = 0; l < mj; ++l){
      int x = counts[l];
      logprobs[l] = tv[x];
    }
    
    CharacterVector nm(mj);
    for (int k = 0; k < mj; ++k) nm[k] = std::to_string(k + 1);
    logprobs.attr("names") = nm;
    
    out[j] = logprobs;
  }
  
  return out;
}