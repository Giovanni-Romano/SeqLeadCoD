#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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