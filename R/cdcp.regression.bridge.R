# Title: Monitoring Bridge
#
# Summary:
#
# TODO:
#
cdcp.regression.bridge.function <- function(data, index_val) {

  n         <- data$n
  y         <- data$y
  X         <- data$X
  Z         <- data$Z
  index     <- data$index
  Dummy     <- data$D # Make Global Dummy instead?

  # From the old code, not optimised for this problem.
  bb_estimation <- foreach::foreach(i = seq_along(index_val), .combine = rbind) %dopar% {

    index_i    <- index <= index_val[i]
    n_i        <- sum(index_i)

    y_left     <- y[index_i]
    X_left     <- X[index_i, ]
    Z_left     <- NULL; if (!is.null(Z)) Z_left <- Z[index_i, ]
    Dummy_left <- NULL; if (!is.null(Dummy)) Dummy_left <- Dummy[index_i, ]

    tmp_lm_i      <- lm(y_left ~ -1 + cbind(X_left, Z_left, Dummy_left))
    tmp_res_i     <- tmp_lm_i$residuals

    tmp_sigma2_i  <- mean(tmp_res_i^2)
    tmp_profile_i <- -(n_i/2)*(log(2*pi) + log(tmp_sigma2_i) + 1)

    c(index_val[i], n_i, tmp_profile_i)
  }

  bb_estimation <- bb_estimation[rowSums(bb_estimation) < Inf, ]

  year_val <- bb_estimation[, 1]
  year_n   <- bb_estimation[, 2]
  year_ll  <- bb_estimation[, 3]

  x0        <- year_val[1]
  x1        <- year_val[length(year_val)]
  y0        <- year_ll[1]
  y1        <- year_ll[length(year_val)]
  b         <- (y1 - y0)/(x1 - x0)
  a         <- -b*x1 + y1

  K         <- year_n[length(year_val)]
  k         <- year_n
  kappa_hat <- sqrt(mean((diff(year_ll)/diff(year_n) - data.table::last(year_ll/year_n))^2))
  # kappa_hat <- 1/sqrt(2)

  BBL       <- K^(-1/2)*(  year_ll[1:length(year_val)] - (a + b*year_val))/kappa_hat
  # BBL       <- K^(-1/2)*(year_ll[1:length(year_val)] - (k/K)*year_ll[length(year_val)])/kappa_hat


  # plot(year_val, year_ll, type = 'l')
  # abline(a = a, b = b, col = 2, lty = 2)

  # ylim <- c(min(c(-1.358, BBL)), max(c(1.358, BBL)))
  # plot(year_val, BBL, type = 'l', lwd = 1.3, main = paste("Region ", i), ylim = ylim)
  # abline(h = 1.358*c(-1, 0, 1), col = 'grey80', lty = c(2, 3, 2), lwd = 3)
  return(BBL)
}

#' Monitoring Bridge
#'
#' This function compiles data and can also be used to add common dummy variables.
#'
#' @param data Data from ..
#' @param index_val A sequence...
#' @return A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales.
#' A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales...
#' An object from `cdcp.regression.estimation(...)` is a list containing the following components:
#' \describe{
#'   \item{ln_max}{First item}
#'   \item{aic}{Second item}
#'   \item{Dn}{ddd}
#' }
#' \itemize{
#'  \item{"parameter 1"}{Stuff}
#'  \item{}{Stuff}
#' }
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' y <- rnorm(100)
#' mu <- 0.5*x + 0.1*5
#' @export
cdcp.regression.bridge <- function(data, index_val) {
  # TODO: It should be possible to have index_val = NULL and then get the maximum range

  if (!all(c("y", "X", "Z", "index", "group") %in% attributes(data)$names)) {
    stop("ERROR: data is not a proper cdcp regression data object. \n")
  }

  if (!all(index_val %in% unique(data$index))) {
    stop("ERROR: Values in index_val are not found in data index.\n")
  }

  return(cdcp.regression.bridge.function(data = data, index_val = index_val))
}

# cdcp.regression.bridge.plot <- function(model) {
#
#   index_val <- model$index_val
#   bbridge_val <- model$bbridge
#
#   ylim <- c(min(c(-1.358, bbridge_val)), max(c(1.358, bbridge_val)))
#   plot(index_val, bbridge_val, type = 'l', lwd = 1.3, ylim = ylim, xlab = "", ylab = "")
#   title(ylab = "Monitoring Bridge", line = 2.5)
#   abline(h = 1.358*c(-1, 0, 1), col = 'grey70', lty = c(2, 3, 2), lwd = 1.6)
#
# }
