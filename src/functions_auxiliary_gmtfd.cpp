#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

NumericMatrix vec_fun_cpp(NumericMatrix x) {
  int rows = x.rows();
  int cols = x.cols();
  
  NumericMatrix x_temp(rows * cols, 1);
  
  for (int i = 0; i < cols; ++i) {
    for (int j = 0; j < rows; ++j) {
      x_temp(i * rows + j, 0) = x(j, i);
    }
  }
  
  return x_temp;
}

NumericMatrix generate_zeros(int rows, int cols) {
  NumericMatrix result(rows, cols);  // Tworzenie macierzy o rozmiarze rows x cols
  
  // Inicjalizacja macierzy zerami
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = 0.0;
    }
  }
  
  return result;
}

#include <iostream>
#include <vector>

// [[Rcpp::export]]
List gr_cov_cpp(std::vector<List> x, List gr_means, int kk, 
                int pp, int ntp, NumericVector n_vec) {
  std::vector<List> x_obs;
  for (int i = 0; i < kk; ++i) {
    List x_temp(n_vec[i]);
    List temp = x[i];
    for (int j = 0; j < n_vec[i]; ++j) {
      NumericMatrix x_temp_m = generate_zeros(pp, ntp);
      for (int l = 0; l < pp; ++l) {
        NumericMatrix temp_l = temp[l];
        NumericMatrix gr_means_l = gr_means[i];
        x_temp_m(l,_) = temp_l(j,_) - gr_means_l(l,_);
      }
      x_temp[j] = x_temp_m;
    }
    x_obs.push_back(x_temp);
  }
  
  List gr_cov(kk);
  for (int i = 0; i < kk; ++i) {
    List temp = x_obs[i];
    arma::mat cov_temp = as<arma::mat>(generate_zeros(pp * ntp, pp * ntp));
    for (int j = 0; j < n_vec[i]; ++j) {
      arma::mat x_ij_stack = as<arma::mat>(vec_fun_cpp(temp[j]));
      cov_temp = cov_temp + x_ij_stack * x_ij_stack.t();
    }
    gr_cov[i] = cov_temp / (n_vec[i] - 1);
  }
  return gr_cov;
}
