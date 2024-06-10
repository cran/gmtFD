
# multiple contrast tests

# For critical values calculation from multi_critical_values_Marc_ls.R
# half partition search, quicker for higher degree of dependence
# This is crit_values2, but without unecessary calculations. So it is faster.
crit_values_ls <- function(data_mat, alpha) {
  # the input is a matrix data_mat of dimenson dim_obs x n_obs containing the B (resampling)
  # observations X_b of dimension dim_obs
  # for our specific purpose X_i is the vector of the test statistics based on the different
  # contrast vectors, i.e. dim_obs = r in the pdf
  # each column represent one observation (i.e. one resampling iteration)
  # the output are the critical values for each dimension (large X_b lead to a rejection)
  # such that the overall level is alpha
  n <- length(data_mat[1, ])
  dimension <- length(data_mat[, 1])

  # First forget the connection of the coordinates, and sort the values per coordinate
  data_order <- t(apply(data_mat, 1, sort))

  # worst case1
  #for each point outside the box only one coordinate leads to a rejection.
  # Thus, for each coordinate (alpha/dim) * number of obs are outside the box  (Bonferoni correction)
  j_low <- ceiling((alpha / dimension) * n) - 1
  # A <- data_mat / data_order[, n - j_low]
  # alpha_low <- mean(apply(A, 2, max) > 1) # count the points being outside the box (where the box borders are given by the critical values)

  # worst case1
  # something like totally linear dependence (in the dimension two)
  j_high <- ceiling(alpha * n)
  # A <- data_mat / data_order[, n - j_low] # ???
  # alpha_high <- mean(apply(A, 2, max) > 1) # count the points being outside the box (where the box borders are given by the critical values)

  # now we search for values j_low and j_high = j_low + 1, such that alpha_low <= alpha and alpha_high > alpha
  while (j_high - j_low > 1) {
    # approx. middle between j_low and j_high
    j_mid <- ceiling(j_low + (j_high - j_low) / 2)
    A <- data_mat / data_order[, n - j_mid]
    alpha_sim <- mean(apply(A, 2, max) > 1)
    ifelse(alpha_sim <= alpha, j_low <- j_mid, j_high <- j_mid)
  }
  # critical values
  return(data_order[, n - j_low])
}

# group mean and covariance calculation
# x - list of k elements corresponding to groups. Each element representing a group 
#     is a list of p elements corresponding to functional variables, and each such element
#     representing a functional variable is a matrix n_i times ntp (number of design time points)
#     of descrete observations in design time points.
mean_cov_m_cpp <- function(x) {
  kk <- length(x)
  pp <- length(x[[1]])
  for (i in 2:kk) {
    if (pp != length(x[[i]])) stop("check the dimensions of functional data in different groups")
  }
  ntp <- ncol(x[[1]][[1]])
  for (i in seq_len(kk)) {
    for (j in seq_len(pp)) {
      if (ntp != ncol(x[[i]][[j]])) stop("check the number of design time points of functional data in different groups and dimensions")
    }
  }
  n_vec <- numeric(kk)
  for (i in seq_len(kk)) n_vec[i] <- nrow(x[[i]][[1]])
  for (i in seq_len(kk)) {
    for (j in seq_len(pp)) {
      if (n_vec[i] != nrow(x[[i]][[j]])) stop("check the number of observations of functional data in different groups")
    }
  }
  nn <- sum(n_vec)
  # sample vectors of mean functions
  gr_means <- vector("list", kk)
  for (i in seq_len(kk)) {
    mean_temp <- matrix(0, pp, ntp)
    for (j in 1:pp) {
      mean_temp[j, ] <- colMeans(x[[i]][[j]])
    }
    gr_means[[i]] <- mean_temp
  }
  gr_means_stack <- gr_means[[1]]
  for (i in 2:kk) {
    gr_means_stack <- rbind(gr_means_stack, gr_means[[i]])
  }
  
  temp_cpp_2 <- gr_cov_cpp(x, gr_means, kk, pp, ntp, n_vec)
  
  return(list(nn = nn, ntp = ntp, kk = kk, n_vec = n_vec, pp = pp,
              gr_means = gr_means_stack, gr_cov = temp_cpp_2))
}

# globalizing and supremum pointwise Hotelling's T^2 -test (GPH) statistic
# H - contrast matrix
# mc - the result of mean_cov_m_cpp() function
gsph_f_m <- function(H, mc) {
  # pointwise gph
  ntp <- mc$ntp
  n_vec <- mc$n_vec
  kk <- mc$kk
  nn <- mc$nn
  pp <- mc$pp
  gph <- numeric(ntp)
  H_gr_means <- H %*% mc$gr_means
  H_gr_means_t <- t(H_gr_means)
  H_t <- t(H)
  for (i in seq_len(ntp)) {
    gamma_hat_t <- matrix(0, kk * pp, kk * pp)
    for (i_k in seq_len(kk)) {
      gamma_hat_t[(i_k * pp + 1 - pp):(i_k * pp), (i_k * pp + 1 - pp):(i_k * pp)] <- 
        nn * mc$gr_cov[[i_k]][(i * pp + 1 - pp):(i * pp), (i * pp + 1 - pp):(i * pp)] / n_vec[i_k]
    }
    gph[i] <- nn * H_gr_means_t[i, ] %*% MASS::ginv(H %*% gamma_hat_t %*% H_t) %*%
      H_gr_means[, i]
  }
  return(c(sum(gph), max(gph)))
}

#' Pointwise Hotelling's \eqn{T^2}-test statistic
#'
#' The function \code{ph_test_statistic()} calculates the pointwise Hotelling's \eqn{T^2}-test statistic.
#'
#' @param x a list of \eqn{k} elements corresponding to groups. Each element representing a group 
#' is a list of \eqn{p} elements corresponding to functional variables, and each such element
#' (representing a functional variable) is a matrix of size \eqn{n_i\times ntp} of descrete observations 
#' in design time points. \eqn{ntp} denotes a number of design time points.
#' @param h contrast matrix. For contrast matrices based on Dunnett’s and Tukey’s contrasts, 
#' it can be created by the \code{contr_mat()} function from the package \code{GFDmcv} (see examples).
#'
#' @details For details, see the documentation of the \code{gmtFD()} function or
#' the papers Munko et al. (2023, 2024).
#'
#' @return A vector of values of the pointwise Hotelling's \eqn{T^2}-test statistic.
#'
#' @examples
#' # Some of the examples may run some time.
#' 
#' # Canadian weather data set
#' # There are three samples of mean temperature and precipitations for
#' # fifteen weather stations in Western Canada, another fifteen in Eastern Canada, 
#' # and the remaining five in Northern Canada.
#' 
#' # one functional variable - temperature
#' library(fda)
#' data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
#' # number of samples
#' k <- 3
#' # number of variables
#' p <- 1
#' # preparing data set
#' gr_label <- rep(c(1, 2, 3), c(15, 15, 5))
#' data_set <- list(list(data_set_t[gr_label == 1, ]),
#'                  list(data_set_t[gr_label == 2, ]),
#'                  list(data_set_t[gr_label == 3, ]))
#'
#' # Tukey's contrast matrix
#' h_tukey <- GFDmcv::contr_mat(k, type = "Tukey")
#' h_tukey_m <- kronecker(h_tukey, diag(p))
#' # plots for pointwise Hotelling's T^2-test statistics
#' oldpar <- par(mfrow = c(2, 2), mar = c(4, 2, 2, 0.1))
#' plot(ph_test_statistic(data_set, h_tukey_m), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Global hypothesis", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_tukey_m[1, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Contrast 1", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_tukey_m[2, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Contrast 2", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_tukey_m[3, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Contrast 3", xlab = "Day")
#' par(oldpar)
#' 
#' # Dunnett's contrast matrix
#' h_dunnett <- GFDmcv::contr_mat(k, type = "Dunnett")
#' h_dunnett_m <- kronecker(h_dunnett, diag(p))
#' # plots for pointwise Hotelling's T^2-test statistics
#' oldpar <- par(mfrow = c(3, 1), mar = c(4, 2, 2, 0.1))
#' plot(ph_test_statistic(data_set, h_dunnett_m), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_dunnett_m))),
#'      main = "Global hypothesis", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_dunnett_m[1, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_dunnett_m))),
#'      main = "Contrast 1", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_dunnett_m[2, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_dunnett_m))),
#'      main = "Contrast 2", xlab = "Day")
#' par(oldpar)
#' 
#' # two functional variables - temperature and precipitation
#' library(fda)
#' data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
#' data_set_p <- t(CanadianWeather$dailyAv[,, "Precipitation.mm"])
#' # number of samples
#' k <- 3
#' # number of variables
#' p <- 2
#' # preparing data set
#' gr_label <- rep(c(1, 2, 3), c(15, 15, 5))
#' data_set <- list(list(data_set_t[gr_label == 1, ], data_set_p[gr_label == 1, ]),
#'                  list(data_set_t[gr_label == 2, ], data_set_p[gr_label == 2, ]),
#'                  list(data_set_t[gr_label == 3, ], data_set_p[gr_label == 3, ]))
#'
#' # Tukey's contrast matrix
#' h_tukey <- GFDmcv::contr_mat(k, type = "Tukey")
#' h_tukey_m <- kronecker(h_tukey, diag(p))
#' # plots for pointwise Hotelling's T^2-test statistics
#' oldpar <- par(mfrow = c(2, 2), mar = c(4, 2, 2, 0.1))
#' plot(ph_test_statistic(data_set, h_tukey_m), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Global hypothesis", xlab = "Day")
#' plot(ph_test_statistic(data_set, h_tukey_m[1:2, ]), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Contrast 1", xlab = "Day")
#' plot(ph_test_statistic(data_set, h_tukey_m[3:4, ]), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Contrast 2", xlab = "Day")
#' plot(ph_test_statistic(data_set, h_tukey_m[5:6, ]), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_tukey_m))),
#'      main = "Contrast 3", xlab = "Day")
#' par(oldpar)
#' 
#' # Dunnett's contrast matrix
#' h_dunnett <- GFDmcv::contr_mat(k, type = "Dunnett")
#' h_dunnett_m <- kronecker(h_dunnett, diag(p))
#' # plots for pointwise Hotelling's T^2-test statistics
#' oldpar <- par(mfrow = c(3, 1), mar = c(4, 2, 2, 0.1))
#' plot(ph_test_statistic(data_set, h_dunnett_m), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_dunnett_m))),
#'      main = "Global hypothesis", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_dunnett_m[1, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_dunnett_m))),
#'      main = "Contrast 1", xlab = "Day")
#' plot(ph_test_statistic(data_set, matrix(h_dunnett_m[2, ], 1)), type = "l",
#'      ylim = c(0, max(ph_test_statistic(data_set, h_dunnett_m))),
#'      main = "Contrast 2", xlab = "Day")
#' par(oldpar)
#'
#' @references Dunnett C. (1955) A multiple comparison procedure for comparing several treatments
#' with a control. Journal of the American Statistical Association 50, 1096-1121.
#' 
#' Munko M., Ditzhaus M., Pauly M., Smaga L., Zhang J.T. (2023) General multiple tests for functional data. 
#' Preprint https://arxiv.org/abs/2306.15259
#'
#' Munko M., Ditzhaus M., Pauly M., Smaga L. (2024) Multiple comparison procedures for simultaneous inference 
#' in functional MANOVA. Preprint https://arxiv.org/abs/2406.01242
#'
#' Tukey J.W. (1953) The problem of multiple comparisons. Princeton University.
#'
#' @import MASS
#' @import fda
#' @import Rcpp
#'
#' @export
ph_test_statistic <- function(x, h) {
  mc <- mean_cov_m_cpp(x)
  ntp <- mc$ntp
  n_vec <- mc$n_vec
  kk <- mc$kk
  nn <- mc$nn
  pp <- mc$pp
  # pointwise gph
  gph <- numeric(ntp)
  H_gr_means <- h %*% mc$gr_means
  H_gr_means_t <- t(H_gr_means)
  H_t <- t(h)
  for (i in seq_len(ntp)) {
    gamma_hat_t <- matrix(0, kk * pp, kk * pp)
    for (i_k in seq_len(kk)) {
      gamma_hat_t[(i_k * pp + 1 - pp):(i_k * pp), (i_k * pp + 1 - pp):(i_k * pp)] <- 
        nn * mc$gr_cov[[i_k]][(i * pp + 1 - pp):(i * pp), (i * pp + 1 - pp):(i * pp)] / n_vec[i_k]
    }
    gph[i] <- nn * H_gr_means_t[i, ] %*% MASS::ginv(h %*% gamma_hat_t %*% H_t) %*%
      H_gr_means[, i]
  }
  return(gph)
}
