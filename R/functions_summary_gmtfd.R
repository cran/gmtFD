# na razie kopia z poprzedniego pakietu

na_przemian <- function(x, y, z, t) {
  temp <- c()
  for (i in seq_len(length(x))) {
    temp <- c(temp, x[i], y[i], z[i], t[i])
  }
  return(temp)
}

#' Print "multifmanova" object
#'
#' Prints the summary of the global and multiple contrasts testing for functional data.
#'
#' @param object a "multifmanova" object.
#' @param ... integer indicating the number of decimal places to be used to present the numerical results.
#' It can be named \code{digits} as in the \code{round()} function (see examples).
#'
#' @details The function prints out the information about the number of functional variables, 
#' number of samples, number of observations in each sample, number of design time points, 
#' contrasts used, test statistics, critical values, \eqn{p}-values of tests performed by the
#' \code{gmtFD()} function. It also gives the decisions.
#'
#' @return No return value, called for side effects.
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
#' # vector of blocks of contrasts labels
#' blocks_contrasts <- rep(1:(nrow(h_tukey_m) / p), each = p)
#' \donttest{
#' # testing without parallel computing
#' res <- gmtFD(data_set, h_tukey_m, blocks_contrasts)
#' summary(res, digits = 3)}
#' \dontshow{
#' data_set <- list(list(data_set_t[gr_label == 1, 1:5]),
#'                  list(data_set_t[gr_label == 2, 1:5]),
#'                  list(data_set_t[gr_label == 3, 1:5]))
#' # testing without parallel computing
#' res <- gmtFD(data_set, h_tukey_m, blocks_contrasts)
#' summary(res, digits = 3)}
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
#' # vector of blocks of contrasts labels
#' blocks_contrasts <- rep(1:(nrow(h_tukey_m) / p), each = p)
#' \donttest{
#' # testing without parallel computing
#' res <- gmtFD(data_set, h_tukey_m, blocks_contrasts)
#' summary(res, digits = 3)}
#'
#' @import fda
#' @import GFDmcv
#'
#' @method summary multifmanova
#' @export

# decisions based on p-values
summary.multifmanova <- function(object, ...) {
  res_global_round <- round(object$res_global, ...)
  temp_global <- data.frame(statistic = c(res_global_round[1, 1], res_global_round[1, 3],
                                          res_global_round[2, 1], res_global_round[2, 3]),
                            p.value = c(res_global_round[1, 2], res_global_round[1, 4], 
                                        res_global_round[2, 2], res_global_round[2, 4]),
                            decision = c(ifelse(object$res_global[1, 2] <= object$alpha, "H1", "H0"),
                                         ifelse(object$res_global[1, 4] <= object$alpha, "H1", "H0"),
                                         ifelse(object$res_global[2, 2] <= object$alpha, "H1", "H0"),
                                         ifelse(object$res_global[2, 4] <= object$alpha, "H1", "H0")))
  rownames(temp_global) <- c("GPH", "mGPH", "SPH", "mSPH")
  res_multi_round <- round(object$res_multi, ...)
  blocks_contrasts_num <- length(unique(object$blocks_contrasts))
  kontrasty <- rep("", 4 * blocks_contrasts_num)
  statystyki <- rep("", 4 * blocks_contrasts_num)
  for (i in seq_len(4 * blocks_contrasts_num)) {
    if (i %% 4 == 1) {
      # kontrasty[i] <- paste("(", paste(object$h[0.5 * i + 0.5, ], collapse = ", "), ")", sep = "")
      kontrasty[i] <- paste("Contrast", 0.25 * i + 0.75, sep = "_")
      statystyki[i] <- res_multi_round[0.25 * i + 0.75, 1]
    }
    if (i %% 4 == 3) {
      statystyki[i] <- res_multi_round[0.25 * i + 0.75, 2]
    }
  }
  temp_multi <- data.frame(contrast = kontrasty,
                           statistic = statystyki,
                           test = rep(c("GPH", "mGPH", "SPH", "mSPH"), times = blocks_contrasts_num),
                           critical_value = na_przemian(res_multi_round[, 3], res_multi_round[, 7], 
                                                        res_multi_round[, 4], res_multi_round[, 8]),
                           p.value = na_przemian(res_multi_round[, 5], res_multi_round[, 9],
                                                 res_multi_round[, 6], res_multi_round[, 10]),
                           decision = ifelse(na_przemian(res_multi_round[, 5], res_multi_round[, 9],
                                                         res_multi_round[, 6], res_multi_round[, 10]) <= object$alpha, "H1", "H0"))
  
  cat("#--- General multiple tests for multivariate functional data ------------#", "\n", "\n")
  cat("- Number of variables:", object$p, "\n")
  cat("- Number of samples:", object$k, "\n")
  cat("- Number of observations in samples:", object$n, "\n")
  cat("- Number of design time points:", object$ntp, "\n")
  cat("- Significance level:", object$alpha)
  cat("\n", "\n")
  cat("#--- Contrasts ----------------------------------------------------------#", "\n")
  if (object$p > 1) {
    gr_names <- rep("", ncol(object$h))
    i_temp <- 0
    for (i_s in seq_len(object$k)) {
      for (i_p in seq_len(object$p)) {
        i_temp <- i_temp + 1
        gr_names[i_temp] <- paste("S", "_", i_s, "_", "V", "_", i_p, sep = "")
      }
    }
    for (i_k in sort(unique(object$blocks_contrasts))) {
      cat(paste("Contrast", i_k, sep = "_"), "\n")
      temp_contrasts <- as.data.frame(object$h[object$blocks_contrasts == i_k, ])
      colnames(temp_contrasts) <- gr_names
      print(temp_contrasts)
    }
    cat("Legend: S_i_V_j - Sample_i_Variable_j")
  } else {
    gr_names <- rep("", ncol(object$h))
    i_temp <- 0
    for (i_s in seq_len(object$k)) {
      i_temp <- i_temp + 1
      gr_names[i_temp] <- paste("S", "_", i_s, sep = "")
    }
    for (i_k in sort(unique(object$blocks_contrasts))) {
      cat(paste("Contrast", i_k, sep = "_"), "\n")
      temp_contrasts <- as.data.frame(matrix(object$h[object$blocks_contrasts == i_k, ], nrow = 1))
      colnames(temp_contrasts) <- gr_names
      print(temp_contrasts)
    }
    cat("Legend: S_i - Sample_i")
  }
  cat("\n")
  cat("\n")
  cat("#--- Overall results ----------------------------------------------------#", "\n")
  print(temp_global)
  cat("\n")
  cat("#--- Multiple contrast testing results ----------------------------------#", "\n")
  print(temp_multi)
  cat("#------------------------------------------------------------------------#", "\n")
}
