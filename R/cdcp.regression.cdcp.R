# Title: CDCPRegression - Confidence Sets for Location of the Change Point
#
# Summary: Functions to compute and plot the confidence sets for the location of a change point. 
#
# TODO:

cdcp.regression.simple <- function(data, model, index_val, boot = 100, index_val_sub = NULL, u = seq(0, 1, length.out = 100), boot_type = c("gaussian", "independent", "group"), cores = 4) {
  return(cdcp.regression.simple.function(data = data, model = model, index_val = index_val, boot = boot, index_val_sub = index_val_sub, u = u, boot_type = boot_type, cores = cores))
}

cdcp.regression.simple.function <- function(data, model, index_val, boot = 100, index_val_sub = NULL, u = NULL, boot_type = c("gaussian", "independent", "group"), cores = 4) {

  doMC::registerDoMC(cores)
  
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

  return(list(cd_set = tau_cc_set, cd_set_half = tau_cc_set_half, index_val = index_val_sub, u = u))
}

cdcp.regression.check <- function(cd) {
  
  if (!all(c("cd_set", "cd_set_half", "index_val", "u") %in% attributes(cd)$names)) {
    stop("ERROR: Provide a proper list of cd values, it must contain cd_set, cd_set_half, index_val and u. \n")
  }
  
  if (!(is.numeric(cd$cd_set) & min(dim(cd$cd_set)) > 0 & length(dim(cd$cd_set)) == 2)) {
    stop("ERROR: cd_set must be an numeric matrix \n")
  }
  
  if (!(is.numeric(cd$cd_set_half) & min(dim(cd$cd_set_half)) > 0 & length(dim(cd$cd_set_half)) == 2)) {
    stop("ERROR: cd_set_half must be an numeric matrix \n")
  }
  
  if (!(is.numeric(cd$index_val) & length(cd$index_val) == dim(cd$cd_set)[1])) {
    stop("ERROR: index_val must be an numeric vector with length equal the number of rows as cd_set \n")
  }

  if (!(is.numeric(cd$u) & length(cd$u) == dim(cd$cd_set)[2])) {
    stop("ERROR: u must be an numeric vector with length equal the number of columns as cd_set \n")
  }
  
}

#' Confidence Sets for the Location of a Change Point
#'
#' This function computes the confidence sets, using bootstrapping, for the location of a change point; see Hermansen (2021) for details. 
#'
#' @param data the output from the `cdcp.regression.data(...)` function.
#' @param model the output of the `cdcp.regression.estimation(...)` function.
#' @param index_val vector of indexes representing potential locations for a change point. 
#' @param boot the number of bootstrap samples.
#' @param index_val_sub a subset of `index_val` that can be used to speed up calculations if the approximate location of the change point is known.
#' @param u the resolution used to compute the confidence sets (a sequence of numbers between 0 and 1). 
#' @param boot_type the method used to sample residuals for the bootstrap: `gaussian` use use parametric 
#' bootstrap under the assumption of a Gaussian model, `independent` makes independent draws from the 
#' estimated residuals and `group` use the assumed group structure to draw independent residuals. 
#' @return A list containing the following:
#' \describe{
#'   \item{cc_set}{a matrix representation of the confidense set for all levels u and index_val.}
#'   \item{cc_set_half}{same as above but with half correction.}
#'   \item{index_val}{same as `index_val` argument}
#'   \item{u}{same as `u` argument}
#' }
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' data <- cdcp.regression.data.sim()
#' model <- cdcp.regression.estimate(data$data, data$index_val)
#' cd <- cdcp.regression(data$data, model, data$index_val)
#' @export
#' @import foreach
#' @import doMC
cdcp.regression <- function(data, model, index_val, boot = 100, index_val_sub = NULL, u = seq(0, 1, length.out = 100), boot_type = c("gaussian", "independent", "group"), cores = 4) {

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

  return(cdcp.regression.simple(data, model, index_val, boot, index_val_sub, u, boot_type, cores = 4))
}


#' Plot of Confidence Sets for the Location of a Change Point
#'
#' Plot the confidence sets for the location of a change point for panel regression. The 
#' plotting function is quite minimal, however, it should be easy to modify this to individual needs.
#'
#' @param cdcp the output of the `cdcp.regression(...)` function.
#' @param half_correction if `TRUE` the confidence sets are plotted with half correction.
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' data <- cdcp.regression.data.sim()
#' model <- cdcp.regression.estimate(data$data, data$index_val)
#' cd <- cdcp.regression(data$data, model, data$index_val)
#' cdcp.regression.plot(cd)
#' @export
cdcp.regression.plot <- function(cd, half_correction = TRUE) {
  cdcp.regression.check(cd)
  
  if (!is.logical(half_correction)) {
    stop("ERROR: half_correction must be logical.\n")
  }
  
  plot_cd_set <- cd$cd_set

  if (half_correction) {
    plot_cd_set <- cd$cd_set_half
  }
  
  matplot(x = cd$index_val, y = plot_cd_set, pch = 21, col = "black", bg = "grey70", xlab = "", ylab = "")
  title(ylab = "Confidence Sets", xlab = "index", line = 2.5)
}

