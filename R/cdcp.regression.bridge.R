# Title: CDCPRegression - Monitoring Bridge
#
# Summary: Compute and plot the monitoring (brownian) bridge plot.
#
# TODO:

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

  return(list(index_val = index_val, bbridge = as.numeric(BBL)))
}

cdcp.regression.bridge.check <- function(bridge) {
  
  if (!all(c("index_val", "bbridge") %in% attributes(bridge)$names)) {
    stop("ERROR: model is not a proper bridge model object. \n")
  }
  
  if (!is.numeric(bridge$index_val)) {
    stop("ERROR: The index_val has to be numerical.\n")
  }
  
  if (!is.numeric(bridge$bbridge)) {
    stop("ERROR: The bbridge has to be numerical.\n")
  }
  
  if (length(bridge$bbridge) != length(bridge$index_val)) {
    stop("ERROR: The bbridge and index_val are of different length. \n")
  }
  
}

#' Monitoring Bridge
#'
#' This function compute the monitoring bridge for a set of indexes (index_val); see Hermansen (2021) for details. 
#'
#' @param data the output from the `cdcp.regression.data(...)` function.
#' @param index_val vector of indexes representing potential locations for a change point. 
#' @return a list consisting of two vectors: 
#' \describe{
#'   \item{index_val}{same as `index_val` argument}
#'   \item{bbridge_val}{the computed monitoring bridge values.}
#' }
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' data <- cdcp.regression.data.sim()
#' model <- cdcp.regression.bridge(data$data, data$index_val)
#' @export
#' @import foreach
#' @import doMC
cdcp.regression.bridge <- function(data, index_val, cores = 4) {
  
  # TODO: It should be possible to have index_val = NULL and then get the maximum range
  
  # Check data in data object
  cdcp.regression.data.check.data(data)
  
  # Check of parameters: 
  if (!is.numeric(index_val)) {
    stop("ERROR: The index_val has to be numerical.\n")
  }  
  
  if (!all(index_val %in% data$index)) {
    stop("ERROR: Some values of index_val are not among index in the data object.\n")
  }  
  
  doMC::registerDoMC(cores)

  return(cdcp.regression.bridge.function(data = data, index_val = index_val))
}

#' Plot of Monitoring Bridge  
#'
#' Plot the so-called monitoring bridge using the output of the `cdcp.regression.bridge(...)` function. 
#'
#' @param bridge The output from running the `cdcp.regression.bridge(...)` function.
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' data <- cdcp.regression.data.sim()
#' bridge <- cdcp.regression.bridge(data$data, data$index_val)
#' cdcp.regression.bridge.plot(bridge)
#' @export
cdcp.regression.bridge.plot <- function(bridge) {

  cdcp.regression.bridge.check(bridge)
  
  index_val <- bridge$index_val
  bbridge <- bridge$bbridge

  ylim <- c(min(c(-1.358, bbridge)), max(c(1.358, bbridge)))
  plot(index_val, bbridge, type = 'l', lwd = 1.3, ylim = ylim, xlab = "", ylab = "")
  title(ylab = "Monitoring Bridge", line = 2.5)
  abline(h = 1.358*c(-1, 0, 1), col = 'grey70', lty = c(2, 3, 2), lwd = 1.6)
}




