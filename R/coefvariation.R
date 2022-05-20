
#' Coefficient of Variation (CV) based on log-transformed data
#'
#' Computes the coefficient of Variation (CV) based on log-transformed data
#'
#' @param x a numeric vector with values that are NOT log-transformed
#' @param na.rm logical. if TRUE then NA values are stripped from x before computation takes place
#'
#' @return numeric vector of length one. If x contains a zero, then NaN is returned. If x is not numeric or integer, NA_real_ is returned
#' @export
#'
#' @examples
#'
CVlog <- function(x, na.rm){
  #sqrt(exp(1)^(sd(log(x, ...))^2)-1) *100
  sqrt(10^(log(10) * stats::sd(log(x, 10), na.rm )^2) - 1) * 100
}
