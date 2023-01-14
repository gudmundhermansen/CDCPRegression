# Title: Confidence Sets for the Change Point
#
# Summary:
#
# TODO:
#
cdcp.regression.simple <- function(data, model, index_val, boot = 100, index_val_sub = NULL, u = seq(0, 1, length.out = 100), boot_type = c("gaussian", "independent", "group")) {
  return(cdcp.regression.simple.function(data = data, model = model, index_val = index_val, boot = boot, index_val_sub = index_val_sub, u = u, boot_type = boot_type))
}

cdcp.regression.simple.function <- function(data, model, index_val, boot = 100, index_val_sub = NULL, u = NULL, boot_type = c("gaussian", "independent", "group")) {

  # The resolution of confidence (0 -> 1) used for the computation of the confidence sets
  if (is.null(u)) {
    u <- seq(0, 1, length.out = 100)
  }

  # If no subset of indices is given the whole set is used
  if (is.null(index_val_sub)) {
    index_val_sub <- index_val
  }


  y               <- data$y
  X               <- data$X
  Z               <- data$Z
  Dummy           <- data$D
  index           <- data$index
  n               <- length(y)

  p               <- data$p
  q               <- data$q; if (is.null(q)) {q <- 0}
  r               <- data$r; if (is.null(r)) {r <- 0}
  l               <- 2*p + q + r


  group           <- data$group
  group_val       <- unique(group)


  tau_betaL       <- model$betaL
  tau_betaR       <- model$betaR
  tau_gamma       <- model$gamma
  tau_dummy       <- model$dummy
  tau_beta        <- c(tau_betaL, tau_betaR)
  tau_sigma       <- model$sigma
  tau_Dn          <- model$Dn
  tau_res         <- model$res


  # Alternatives for simulating the residuals for the bootstrap
  if (boot_type[1] == "gaussian") {
    epsilon <- tau_sigma*rnorm(n = boot*n)
    epsilon <- matrix(epsilon, nrow = n, ncol = boot)
  } else if (boot_type[1] == "independent") {
    epsilon <- matrix(sample(x = tau_res, size = boot*n, replace = TRUE), nrow = n, ncol = boot)
  } else if (boot_type[1] == "group") {
    epsilon <- matrix(NA, nrow = n, ncol = boot, byrow = FALSE)
    for (group_i in group_val) {
      index_i            <- group == group_i
      n_i                <- sum(index_i)
      res_group_i        <- tau_res[index_i]
      epsilon[index_i, ] <- sample(res_group_i, size = boot*n_i, replace = TRUE)
    }
  }
  # ---------------------------------------------------------------------------------------




  # Pre-compute values  -------------------------------------------------------------------

  # The non-chang point part of the expectation
  mu_gamma <- cbind(Z, Dummy); if (is.null(mu_gamma)) mu_gamma <- rep(0, n) else mu_gamma <- mu_gamma%*%c(tau_gamma, tau_dummy)

  # Pre-computing the inverse of XtX
  tmp   <- cbind(Z, Dummy); if (is.null(tmp)) {tmp <- matrix(0, nrow = n, ncol = 1)}
  A     <- NULL
  B     <- tmp
  C     <- t(tmp)
  D     <- crossprod(tmp)
  if (length(D) == 1 & D[1, 1] == 0) {D_inv <- 0} else {D_inv <- solve(D)}

  XtX_k_inv <- foreach::foreach(k = seq_along(index_val_sub), .combine = cbind) %do% {
    if (k == 1) {
      k              <- 1
      index_k        <- index_val_sub[k]

      X_L            <- X*0
      index_L        <- index <= index_k
      X_L[index_L, ] <- X[index_L, ]

      X_R            <- X*0
      index_R        <- index > index_k
      X_R[index_R, ] <- X[index_R, ]

      X_k            <- cbind(X_L, X_R)
    } else {
      index_k        <- index_val_sub[k]
      tmp_k          <- which(index == index_k)

      x_k            <- X[tmp_k, ]
      X_k[tmp_k, ]   <- c(x_k, x_k*0)
    }

    tX_k    <- t(X_k)

    A_k     <- tX_k%*%X_k
    B_k     <- tX_k%*%B
    C_k     <- C%*%X_k

    A_k_inv <-  solve(A_k - B_k%*%(D_inv%*%C_k))
    B_k_inv <- -A_k_inv%*%B_k%*%D_inv
    C_k_inv <- -(D_inv%*%C_k)%*%A_k_inv
    D_k_inv <-  D_inv - C_k_inv%*%(B_k%*%D_inv)

    return(rbind(cbind(A_k_inv, B_k_inv), cbind(C_k_inv, D_k_inv))[1:l, 1:l])
  }
  # ---------------------------------------------------------------------------------------



  # pb <- txtProgressBar(min = 0, max = length(index_val_sub), style = 3)
  #
  # sim <- foreach::foreach(i = which(index_val %in% index_val_sub), .combine = cbind) %do% {
  #
  #   index_i       <- index_val[i]
  #
  #   mu_i          <- y + NA
  #
  #   index_L       <- index <= index_i
  #   mu_i[index_L] <- as.matrix(X[index_L, ])%*%tau_betaL
  #
  #   index_R       <- index > index_i
  #   mu_i[index_R] <- as.matrix(X[index_R, ])%*%tau_betaR
  #
  #   mu_i          <- mu_i + mu_gamma
  #
  #
  #   sim_j <- foreach::foreach(j = 1:boot, .combine = cbind) %dopar% {
  #
  #     y_j     <- mu_i + epsilon[, j]
  #     tmp_y_j <- NULL; if (!is.null(tmp)) tmp_y_j <- crossprod(tmp, y_j)
  #
  #     profile_k <- foreach::foreach(k = seq_along(index_val), .combine = c) %do% {
  #
  #       index_k <- index_val[k]
  #
  #       if (k == 1) {
  #         X_L            <- X*0
  #         index_L        <- index <= index_k
  #         X_L[index_L, ] <- X[index_L, ]
  #
  #         X_R            <- X*0
  #         index_R        <- index > index_k
  #         X_R[index_R, ] <- X[index_R, ]
  #
  #         X_k            <- cbind(X_L, X_R)
  #       } else {
  #         tmp_k          <- which(index == index_k)
  #         x_k            <- X[tmp_k, ]
  #         X_k[tmp_k, ]   <- c(x_k, x_k*0)
  #       }
  #
  #       tX_k      <- t(X_k)
  #
  #       A_k       <- tX_k%*%X_k
  #       B_k       <- tX_k%*%B
  #       C_k       <- C%*%X_k
  #
  #       A_k_inv   <- solve(A_k - B_k%*%D_inv%*%C_k)
  #       B_k_inv   <- -A_k_inv%*%B_k%*%D_inv
  #       C_k_inv   <- -D_inv%*%C_k%*%A_k_inv
  #       D_k_inv   <- D_inv - C_k_inv%*%(B_k%*%D_inv)
  #
  #       XtX_k_inv <- rbind(cbind(A_k_inv, B_k_inv), cbind(C_k_inv, D_k_inv))
  #
  #       beta_k    <- XtX_k_inv%*%rbind(tX_k%*%y_j, tmp_y_j)
  #       mu_k      <- cbind(X_k, tmp)%*%beta_k
  #
  #       res_k     <- y_j - mu_k
  #       sigma2_k  <- mean(res_k^2)
  #
  #       return(-(n/2)*(log(2*pi) + log(sigma2_k) + 1))
  #     }
  #     return(profile_k)
  #   }
  #
  #   setTxtProgressBar(pb, i)
  #
  #   return(apply(sim_j, 2, function(x) {2*(max(x, na.rm = TRUE) - x)})[index_val == index_i])
  # }
  # close(pb)

  tmp   <- cbind(Z, Dummy)


  sim <- foreach::foreach(i = seq_along(index_val_sub), .combine = cbind) %dopar% {

    index_i       <- index_val_sub[i]

    mu_i          <- y + NA

    index_L       <- index <= index_i
    mu_i[index_L] <- as.matrix(X[index_L, ])%*%tau_betaL

    index_R       <- index > index_i
    mu_i[index_R] <- as.matrix(X[index_R, ])%*%tau_betaR

    mu_i          <- mu_i + mu_gamma

    Y_i           <- matrix(mu_i, nrow = n, ncol = boot) + epsilon
    tmp_Y_i <- NULL; if (!is.null(tmp)) tmp_Y_i <- crossprod(tmp, Y_i)


    sim_i <- foreach::foreach(k = seq_along(index_val_sub), .combine = rbind) %do% {

      index_k <- index_val_sub[k]

      if (k == 1) {
        X_L            <- X*0
        index_L        <- index <= index_k
        X_L[index_L, ] <- X[index_L, ]

        X_R            <- X*0
        index_R        <- index > index_k
        X_R[index_R, ] <- X[index_R, ]

        X_k            <- cbind(X_L, X_R)
      } else {
        tmp_k          <- which(index == index_k)
        x_k            <- X[tmp_k, ]
        X_k[tmp_k, ]   <- c(x_k, x_k*0)
      }

      res_k     <- Y_i - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(t(X_k)%*%Y_i, tmp_Y_i))
      sigma2_k  <- colMeans(res_k^2)

      return(-(n/2)*(log(2*pi) + log(sigma2_k) + 1))
    }
    return(apply(sim_i, 2, function(x) {2*(max(x, na.rm = TRUE) - x)})[index_val_sub == index_i])
  }


  tau_cc          <- index_val_sub + NA
  tau_cc_half     <- index_val_sub + NA
  tau_cc_set      <- matrix(u, nrow = length(index_val_sub), ncol = length(u), byrow = TRUE)
  tau_cc_set_half <- matrix(u, nrow = length(index_val_sub), ncol = length(u), byrow = TRUE)


  for (i in seq_along(index_val_sub)) {
    tau_cc_half[i]            <- sum(sim[, i] < tau_Dn[index_val_sub[i] == index_val])/boot + 0.5*mean(sim[, i] == tau_Dn[i])
    index                     <- u < tau_cc_half[i]
    tau_cc_set_half[i, index] <- NA
  }


  # No half correction
  for (i in seq_along(index_val_sub)) {
    tau_cc[i]            <- sum(sim[, i] < tau_Dn[index_val_sub[i] == index_val])/boot
    index                <- u < tau_cc[i]
    tau_cc_set[i, index] <- NA
  }

  return(list(cc_set = tau_cc_set, cc_set_half = tau_cc_set_half, index_val = index_val_sub, u = u))
}


cdcp.regression.check <- function(cc) {
  
  if (!all(c("cc_set", "cc_set_half", "index_val", "u") %in% attributes(cc)$names)) {
    stop("ERROR: Provide a proper list of cc values, it must contain cc_set, cc_set_half, index_val and u. \n")
  }
  
  if (!(is.numeric(cc$cc_set) & min(dim(cc$cc_set)) > 0 & length(dim(cc$cc_set)) == 2)) {
    stop("ERROR: cc_set must be an numeric matrix \n")
  }
  
  if (!(is.numeric(cc$cc_set_half) & min(dim(cc$cc_set_half)) > 0 & length(dim(cc$cc_set_half)) == 2)) {
    stop("ERROR: cc_set_half must be an numeric matrix \n")
  }
  
  if (!(is.numeric(cc$index_val) & length(cc$index_val) == dim(cc$cc_set)[1])) {
    stop("ERROR: index_val must be an numeric vector with length equal the number of rows as cc_set \n")
  }

  if (!(is.numeric(cc$u) & length(cc$u) == dim(cc$cc_set)[2])) {
    stop("ERROR: u must be an numeric vector with length equal the number of columns as cc_set \n")
  }
  
}




#' Confidence Distribution for the Location of the Change Point
#'
#' `cdcp.regression(...)` compute a bootstrap for the change point
#'
#' @param data A data object/list that specify the structure of the model created by the `cdcp.regression.data(...)`.
#' @param model The estimated model object from `cdcp.regression.estimation(...)`.
#' @param index_val A consecutive sequence of indexes representing the locations for a potential change point.
#' @param boot The number of bootstrap samples.
#' @param index_val_sub A subset of `index_val` used to speed up the calculations.
#' @param u The resolution for the probability...
#' @param boot_type The method used to sample residuales for the bootstrap `= "group"` ...
#' @return A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales.
#' A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales...
#' An object from `cdcp.regression.estimation(...)` is a list containing the following components:
#' \describe{
#'   \item{cc_set}{First item}
#'   \item{cc_set_half}{Second item}
#'   \item{index_val}{ddd}
#'   \item{u}{ddd}
#' }
#' \itemize{
#'  \item{"parameter 1"}{Stuff}
#'  \item{}{Stuff}
#' }
#'
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @export
cdcp.regression <-  function(data, model, index_val, boot = 100, index_val_sub = NULL, u = seq(0, 1, length.out = 100), boot_type = c("gaussian", "independent", "group")) {

  cdcp.regression.data.check.data(data)
  
  cdcp.regression.estimate.check.model(model)
  
  if (!all(index_val %in% unique(data$index))) {
    stop("ERROR: Values in index_val are not found in data index.\n")
  }

  if (!is.null(index_val_sub)) {
    if (!all(index_val_sub %in% index_val)) {
      stop("ERROR: index_val_sub is not a subset of index_val.\n")
    }
  }
  
  # if (is.null(index_val_sub)) {
  #   cat("WARNING: It is often demanding to run bootstrap for all indexes. If possible, use index_val_sub to specify a subset of indexes.\n")
  # }
  
  if (!(is.numeric(boot) & (round(boot) - boot) == 0 & boot > 0)) {
    stop("ERROR: boot must be a positive integer. \n")
  }

  if (!is.numeric(u) | min(u) < 0 | min(u) > 1) {
    stop("ERROR: The values in u must be between 0 and 1. \n")
  }

  if (!(boot_type[1] %in% c("gaussian", "independent", "group"))) {
    stop("ERROR: Incorrect boot type, should be either gaussian, independent or group.\n")
  }

  return(cdcp.regression.simple(data, model, index_val, boot, index_val_sub, u, boot_type))
}





#' Plot ... Confidence Distribution for the Location of the Change Point
#'
#' `cdcp.regression(...)` compute a bootstrap for the change point
#'
#' @param data A data object/list that specify the structure of the model created by the `cdcp.regression.data(...)`.
#' @param model The estimated model object from `cdcp.regression.estimation(...)`.
#' @param index_val A consecutive sequence of indexes representing the locations for a potential change point.
#' @param boot The number of bootstrap samples.
#' @param index_val_sub A subset of `index_val` used to speed up the calculations.
#' @param u The resolution for the probability...
#' @param boot_type The method used to sample residuales for the bootstrap `= "group"` ...
#' @return A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales.
#' A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales...
#' An object from `cdcp.regression.estimation(...)` is a list containing the following components:
#' \describe{
#'   \item{cc_set}{First item}
#'   \item{cc_set_half}{Second item}
#'   \item{index_val}{ddd}
#'   \item{u}{ddd}
#' }
#' \itemize{
#'  \item{"parameter 1"}{Stuff}
#'  \item{}{Stuff}
#' }
#'
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @export
cdcp.regression.plot <- function(cc, half_correction = TRUE) {
  cdcp.regression.check(cc)
  
  if (!is.logical(half_correction)) {
    stop("ERROR: half_correction must be logical.\n")
  }
  
  plot_cc_set <- cc$cc_set

  if (half_correction) {
    plot_cc_set <- cc$cc_set_half
  }
  
  matplot(x = cc$index_val, y = plot_cc_set, pch = 21, col = "black", bg = "grey70", xlab = "", ylab = "")
  title(ylab = "Confidence Sets", xlab = "index", line = 2.5)
}







if (FALSE) {

  set.seed(1133462)

  par(bty = 'l')

  country_names <- c("Mirasius", "Oras")
  country_names <- sort(country_names)


  year_min <- 1900
  year_max <- 2000
  year_seq <- year_min:year_max


  n        <- length(year_seq)
  m        <- length(country_names)

  beta_0    <- 0.35
  beta_1    <- 0.05

  tau_0     <- 1950
  gamma_0   <- 0

  u           <- gamma_0*rep(year_seq, times = m)/year_max
  mu_year_seq <- beta_0 + beta_1*(year_seq > tau_0) + gamma_0*year_seq/year_max
  mu          <- rep(beta_0 + beta_1*(year_seq > tau_0), times = m) + u
  sigma       <- rep(c(0.01, 0.2), each = n)

  y         <- mu + sigma*rnorm(n*m)
  X         <- as.matrix(1 + mu*0)
  Z         <- NULL

  index     <- rep(year_seq, times = m)
  group     <- rep(country_names, each = n)
  country   <- rep(country_names, each = n)
  id        <- rep(1:m, each = n)

  index_val <- 1910:1990

  data      <- cdcp.regression.clean.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = FALSE, group_dummy = FALSE, lag = FALSE)
  model     <- cdcp.regression.estimate(data = data, index_val = index_val)


  par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
  sim <- cdcp.regression.simple.function(data, model, boot = 500, index_val = index_val, index_val_sub = 1930:1975, u = seq(0, 1, length.out = 65), boot_type = "group")

  #par(mar = c(2, 2, 1.6, 0.4), mfrow = c(1, 1))
  matplot(x = sim$tau_val, sim$cc_set, pch = 21, bg = 'grey70', col = 'black', axes = FALSE, xlab = '', ylab = '')

  # title(main = region_names[r])
  title(ylab = "Confidence Sets", line = 2.5)
  title(xlab = "Year", line = 2.8)
  axis(1)
  axis(2)
  box(bty = 'l')
  abline(h = 0.95, lty = 2, lwd = 1.7, col = 'black')
  abline(v = tau_0, lty = 1, lwd = 1.7, col = 'black')



}



if (FALSE) {

  # More stuff here for the bootstrap:
  #set.seed(1)
  # sd_group <- TRUE
  #
  # epsilon <- matrix(NA, nrow = n, ncol = boot, byrow = FALSE)
  #
  # if (sd_group) {
  #   for (group_i in group_val) {
  #     index_i            <- group_seq == group_i
  #     n_i                <- sum(index_i)
  #     res_group_i        <- tau_res[index_i]
  #     epsilon[index_i, ] <- sample(res_group_i, size = boot*n_i, replace = TRUE)
  #   }
  # }
  #

  # if (bboot) {
  #   k       <- 10
  #   epsilon_table <- matrix(NA, nrow = n - k + 1, ncol = k)
  #   for (i in 1:(n - k + 1)) {
  #     epsilon_table[i, ] <- res[i + 0:(k - 1)]
  #   }
  #
  #   index <- sample.int(n - k + 1, ceiling((n*B)/k), replace = TRUE)
  #
  #   epsilon <- NULL
  #   for (i in index) {
  #     epsilon <- c(epsilon, epsilon_table[i, ])
  #   }
  #   epsilon <- epsilon[1:(n*B)]
  # }

  # epsilon <- tau_sigma*rnorm(n = boot*n)
  # epsilon <- matrix(epsilon, nrow = n, ncol = boot, byrow = FALSE)


  # set.seed(1)
  # epsilon <- tau_sigma*rnorm(n = boot*n)
  # epsilon <- matrix(epsilon, nrow = n, ncol = boot, byrow = FALSE)



  data_raw         <- fread(input = '~/Dropbox/shared/Change points/Data/hvdem.csv')
  data_raw_country <- unique(data_raw$country_name)
  data_raw         <- data_raw[country_name %in% data_raw_country[33:35]]

  y         <- log(data_raw$v2x_polyarchy)
  X         <- cbind(data_raw$Maddison_gdppc_1990_estimate_lag)
  Z         <- cbind(data_raw$Maddison_pop_estimate)
  index     <- data_raw$year
  group     <- data_raw$country_name

  data      <- cdcp.regression.clean.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = TRUE, group_dummy = TRUE, lag = TRUE)
  index_val <- cdcp.regression.clean.data.find.index.val(data, threshold = 5)$index_val


  model     <- cdcp.regression.estimate(data = data, index_val = index_val)
  cc        <- cdcp.regression.simple.function(data = data, model = model, B = 100, index_val = index_val, tmp_val = index_val)


  # delta_val <- seq(0.05, 0.3, length.out = 50)
  # doc       <- cdcp.regression.beta.doc.function(data = data, model = model, k = 1, delta_val = delta_val, index_val = index_val, B = 100)

  # with(doc, plot(delta_val, cc_approx, type = 'l'))
  # with(doc, lines(delta_val, cc, type = 'l', col = "green"))

  # lines(doc$delta_val, cc, type = 'l', col = "green")


  #
  # print(system.time(bbridge <- cdcp.regression.bridge(data = data, index_val = index_val)))
  # model$bbridge <- bbridge


  # tmp <- cdcp.regression.robust.standard.errors(data = data, model = model)
  # tmp$sigma_gamma_robust
  # model$sigma_gamma



  # cdcp.regression.summary.plot.index(data, model)
  # cdcp.regression.summary.plot.group(data, model)
  # cdcp.regression.summary.plot.para(data, model)

  # cdcp.regression.summary.plot.residuales.index(data, model)
  # cdcp.regression.summary.plot.residuales.fit(data, model)
  # cdcp.regression.summary.plot.model(data, model)
  # cdcp.regression.summary.plot.data(data, model)

  # cdcp.regression.summary.plot.deviance(data, model)

  # print(system.time(sim <- cdcp.regression.function(data, model, B = 100, index_val = index_val, tmp_val = 1915:1950)))
  # image(x = sim$tau_val, y = sim$u, sim$cc_set)

  # print(system.time(foreach::foreach (i=1:5, .combine='c') %dopar% {Sys.sleep(2);i}))
  # print(system.time(foreach::foreach (i=1:5, .combine='c') %do% {Sys.sleep(2);i}))

}





# Old, but more optimised version, however, this has to be fixed.
# cdcp.regression.function <- function(data, model, B, index_val, tmp_val, u = NULL, bboot = FALSE) {
#
#   if (is.null(u)) {
#     u <- seq(0, 1, length.out = 100)
#   }
#
#   boot            <- B
#   y               <- data$y
#   X               <- data$X
#   Z               <- data$Z
#   Dummy           <- data$D
#
#   p               <- data$p
#   q               <- data$q; if (is.null(q)) {q <- 0}
#   r               <- data$r; if (is.null(r)) {r <- 0}
#   l               <- 2*p + q + r
#
#   group_seq       <- data$group
#   group_val       <- unique(group_seq)
#
#   tau_val         <- index_val
#   tau_seq         <- data$index
#   n               <- length(y)
#
#   tau_betaL       <- model$betaL
#   tau_betaR       <- model$betaR
#   tau_gamma       <- model$gamma
#   tau_dummy       <- model$dummy
#   tau_beta        <- c(tau_betaL, tau_betaR)
#   tau_sigma       <- model$sigma
#   tau_Dn          <- model$Dn
#   tau_res         <- model$res
#
#
#
#   # Sampeling Residuales ------------------------------------------------------------------
#   set.seed(1)
#   epsilon <- tau_sigma*rnorm(n = boot*n)
#   epsilon <- matrix(epsilon, nrow = n, ncol = boot, byrow = FALSE)
#   # ---------------------------------------------------------------------------------------
#
#
#   # That do not depend on the change point ------------------------------------------------
#   mu_gamma <- cbind(Z, Dummy); if (is.null(mu_gamma)) mu_gamma <- rep(0, n) else mu_gamma <- mu_gamma%*%c(tau_gamma, tau_dummy)
#
#   tmp   <- cbind(Z, Dummy); if (is.null(tmp)) {tmp <- matrix(0, nrow = n, ncol = 1)}
#   A     <- NULL
#   B     <- tmp
#   C     <- t(tmp)
#   D     <- crossprod(tmp)
#   if (length(D) == 1 & D[1, 1] == 0) {D_inv <- 0} else {D_inv <- solve(D)}
#   # ---------------------------------------------------------------------------------------
#
#
#
#
#
#   # Pre-computing the inverse of XtX ------------------------------------------------------
#   XtX_k_inv <- foreach::foreach(k = seq_along(tau_val), .combine = cbind) %do% {
#     if (k == 1) {
#       k              <- 1
#       tau_k          <- tau_val[k]
#
#       X_L            <- X*0
#       index_L        <- tau_seq <= tau_k
#       X_L[index_L, ] <- X[index_L, ]
#
#       X_R            <- X*0
#       index_R        <- tau_seq > tau_k
#       X_R[index_R, ] <- X[index_R, ]
#
#       X_k            <- cbind(X_L, X_R)
#     } else {
#       tau_k          <- tau_val[k]
#       index_k        <- which(tau_seq == tau_k)
#
#       x_k            <- X[index_k, ]
#       X_k[index_k, ] <- c(x_k, x_k*0)
#     }
#
#     tX_k      <- t(X_k)
#
#     A_k       <- tX_k%*%X_k
#     B_k       <- tX_k%*%B
#     C_k       <- C%*%X_k
#
#     A_k_inv   <-  solve(A_k - B_k%*%(D_inv%*%C_k))
#     B_k_inv   <- -A_k_inv%*%B_k%*%D_inv
#     C_k_inv   <- -(D_inv%*%C_k)%*%A_k_inv
#     D_k_inv   <-  D_inv - C_k_inv%*%(B_k%*%D_inv)
#
#     return(rbind(cbind(A_k_inv, B_k_inv), cbind(C_k_inv, D_k_inv))[1:l, 1:l])
#   }
#   # ---------------------------------------------------------------------------------------
#
#
#
#
#
#   ### TEST - computing the epsilons residuals... ###
#
#   tmp         <- cbind(Z, Dummy)
#   tmp_epsilon <- NULL; if (!is.null(tmp)) tmp_epsilon <- crossprod(tmp, epsilon)
#
#   k              <- 1
#   tau_k          <- tau_val[k]
#
#   X_L            <- X*0
#   index_L        <- tau_seq <= tau_k
#   X_L[index_L, ] <- X[index_L, ]
#
#   X_R            <- X*0
#   index_R        <- tau_seq > tau_k
#   X_R[index_R, ] <- X[index_R, ]
#
#   X_k            <- cbind(X_L, X_R)
#
#   res_epsilon <- foreach::foreach(k = seq_along(tau_val), .combine = rbind) %do% {
#     if (k > 1) {
#       tau_k          <- tau_val[k]
#       index_k        <- which(tau_seq == tau_k)
#       x_k            <- X[index_k, ]
#       X_k[index_k, ] <- c(x_k, x_k*0)
#     }
#     return(epsilon - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(crossprod(X_k, epsilon), tmp_epsilon)))
#   }
#   ####
#
#   # tmp_mu <- NULL; if (!is.null(tmp)) tmp_mu <- crossprod(tmp, mu)
#
#   # XtX_k_inv_tmp <- matrix(c(0, 0.1, 0.1, 1), 2, 2)
#   # XtX_k_inv_tmp <- XtX_k_inv[1:l, 1:l + (k - 1)*l]
#
#   res_mu <- foreach::foreach(k = seq_along(tau_val), .combine = rbind) %do% {
#     tau_i  <- tau_val[k]
#     if (k == 1) {
#       mu_i          <- y + NA
#
#       index_L       <- tau_seq <= tau_i
#       mu_i[index_L] <- as.matrix(X[index_L, ])%*%tau_betaL
#
#       index_R       <- tau_seq > tau_i
#       mu_i[index_R] <- as.matrix(X[index_R, ])%*%tau_betaR
#     } else {
#       tau_i         <- tau_val[i]
#       index_i       <- tau_seq == tau_i
#       mu_i[index_i] <- as.matrix(X[index_i, ])%*%tau_betaL
#     }
#
#     if (k > 1) {
#       tau_k          <- tau_val[k]
#       index_k        <- which(tau_seq == tau_k)
#       x_k            <- X[index_k, ]
#       X_k[index_k, ] <- c(x_k, x_k*0)
#     }
#
#     tmp_mu_i <- NULL; if (!is.null(tmp)) tmp_mu_i <- crossprod(tmp, mu_i)
#
#
#     return(mu_i - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(crossprod(X_k, mu_i), tmp_mu_i)))
#   }
#
#
#
#   mu_tau <- foreach::foreach(i = seq_along(tau_val), .combine = cbind) %do% {
#     tau_i  <- tau_val[i]
#     if (i == 1) {
#       mu_i          <- y + NA
#
#       index_L       <- tau_seq <= tau_i
#       mu_i[index_L] <- as.matrix(X[index_L, ])%*%tau_betaL
#
#       index_R       <- tau_seq > tau_i
#       mu_i[index_R] <- as.matrix(X[index_R, ])%*%tau_betaR
#     } else {
#       tau_i         <- tau_val[i]
#       index_i       <- tau_seq == tau_i
#       mu_i[index_i] <- as.matrix(X[index_i, ])%*%tau_betaL
#     }
#     return(mu_i)
#   }
#
#   mu_tau <- mu_tau + c(mu_gamma)
#
#   # New, we do not need the below anymore?
#   tmp    <- cbind(Z, Dummy)
#
#
#
#   sim <- foreach::foreach(i = which(tau_val %in% tmp_val), .combine = cbind) %dopar% {
#
#     tau_i          <- tau_val[i]
#
#     mu_i           <- mu_tau[, i]
#     y_i            <- mu_tau[, i] + epsilon[, i]
#
#     # New, we do not need the below anymore?
#     tmp_mu_i <- NULL; if (!is.null(tmp)) tmp_mu_i <- crossprod(tmp, mu_i)
#     tmp_y_i <- NULL; if (!is.null(tmp)) tmp_y_i <- crossprod(tmp, y_i)
#
#
#     k              <- 1
#     tau_k          <- tau_val[k]
#
#     X_L            <- X*0
#     index_L        <- tau_seq <= tau_k
#     X_L[index_L, ] <- X[index_L, ]
#
#     X_R            <- X*0
#     index_R        <- tau_seq > tau_k
#     X_R[index_R, ] <- X[index_R, ]
#
#     X_k            <- cbind(X_L, X_R)
#
#
#     sim_i <- foreach::foreach(k = seq_along(tau_val), .combine = rbind) %do% {
#
#       if (k > 1) {
#         tau_k          <- tau_val[k]
#         index_k        <- which(tau_seq == tau_k)
#         x_k            <- X[index_k, ]
#         X_k[index_k, ] <- c(x_k, x_k*0)
#       }
#
#       # From the estimation:
#       # res_i      <- y - cbind(X_i, tmp)%*%XtX_i_inv%*%rbind(tX_i%*%y, tmp_y)
#
#       # Original:
#       res_k_org  <- y_i - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(crossprod(X_k, y_i), tmp_y_i))
#       sigma2_k   <- colMeans(res_k_org^2)
#
#
#       # res_k <- mu_i - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(crossprod(X_k, mu_i), tmp_mu_i))
#
#       # res_k <- mu_i - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(crossprod(X_k, mu_i), tmp_mu_i))
#       # res_k <- res_epsilon[1:n + (k - 1)*n, 1:boot] + c(res_k)
#       # res_k <- res_epsilon[1:n + (k - 1)*n, 1:boot]
#       # sigma2_k   <- colMeans(res_k^2)
#
#       # res_k_tmp  <- epsilon + c(mu_i) - cbind(X_k, tmp)%*%(XtX_k_inv[1:l, 1:l + (k - 1)*l]%*%rbind(t(X_k), t(tmp)))%*%(epsilon + c(mu_i))
#       # sigma2_k   <- colMeans(res_k_tmp^2)
#
#       return(-(n/2)*(log(2*pi) + log(sigma2_k) + 1))
#
#     }
#
#
#     return(apply(sim_i, 2, function(x) {2*(max(x, na.rm = TRUE) - x)})[tau_val == tau_i])
#   }
#
#
#
#   tau_cc          <- tmp_val + NA
#   tau_cc_half     <- tmp_val + NA
#   tau_cc_set      <- matrix(u, nrow = length(tmp_val), ncol = length(u), byrow = TRUE)
#   tau_cc_set_half <- matrix(u, nrow = length(tmp_val), ncol = length(u), byrow = TRUE)
#
#
#   for (i in seq_along(tmp_val)) {
#     tau_cc_half[i]            <- sum(sim[, i] < tau_Dn[tmp_val[i] == tau_val])/boot + 0.5*mean(sim[, i] == tau_Dn[i])
#     index                     <- u < tau_cc_half[i]
#     tau_cc_set_half[i, index] <- NA
#   }
#
#   for (i in seq_along(tmp_val)) {
#     tau_cc[i]            <- sum(sim[, i] < tau_Dn[tmp_val[i] == tau_val])/boot + 0.0*mean(sim[, i] == tau_Dn[i]) # removed half correction
#     index                <- u < tau_cc[i]
#     tau_cc_set[i, index] <- NA
#   }
#
#
#   return(list(cc_set = tau_cc_set, cc_set_half = tau_cc_set_half, tau_val = tmp_val, u = u))
# }








