# Title: Confidence Distribution for Change Point for Panel Regression - Estimation 
#
# Summary: The main functions used for estimation of parameters for the CDCPRegression package. Some 
# clever tricks are needed to make the estimation efficient for large-ish data and high number of 
# parameters. 
#
# TODO:
# 1) Tests.
# 2) Test of indexes
# 2) Use similar method as for the CC and CS simulation to speed up calculations

cdcp.regression.estimate.check.model <- function(model) {
  if (!all(c("ln_max", "aic", "Dn", "beta", "betaL", "betaR", 
             "tau", "gamma", "dummy", "mu", 
             "sigma", "sigma_beta", "sigma_betaL", "sigma_betaLR", "sigma_gamma", "sigma_dummy", 
             "res", "index_val") %in% attributes(model)$names)) {
    stop("ERROR: data is not a proper cdcp regression data object. \n")
  }
}


cdcp.regression.estimate.function <- function(data, index_val = NULL) {

  m     <- length(index_val)

  y     <- data$y
  n     <- length(y)
  X     <- data$X

  Z     <- data$Z
  Dummy <- data$D

  index <- data$index

  tmp   <- cbind(Z, Dummy); if (is.null(tmp)) {tmp <- matrix(0, nrow = n, ncol = 1)}
  A     <- NULL
  B     <- tmp
  C     <- t(tmp)
  D     <- t(tmp)%*%tmp
  if (length(D) == 1 & D[1, 1] == 0) {D_inv <- 0} else {D_inv <- solve(D)}

  tmp_y  <- t(tmp)%*%y

  p      <- data$p
  q      <- data$q; if (is.null(q)) q <- 0
  r      <- data$r; if (is.null(r)) r <- 0

  index_para    <- matrix(NA, nrow = m, ncol = 2*p + q + r)
  index_para_sd <- matrix(NA, nrow = m, ncol = 2*p + q + r)

  index_profile <- index_val + NA
  index_mu      <- matrix(NA, nrow = m, ncol = n)



  for (i in seq_along(index_val)) {

    if (i == 1) {
      X_L            <- X*0
      index_L        <- index <= index_val[i]
      X_L[index_L, ] <- X[index_L, ]

      X_R            <- X*0
      index_R        <- index > index_val[i]
      X_R[index_R, ] <- X[index_R, ]

      X_i            <- cbind(X_L, X_R)
    } else {
      index_i        <- which(index == index_val[i])
      x_i            <- X[index_i, ]
      X_i[index_i, ] <- c(x_i, x_i*0)
    }

    tX_i      <- t(X_i)

    A_i       <- tX_i%*%X_i
    B_i       <- tX_i%*%B
    C_i       <- C%*%X_i

    A_i_inv   <- solve(A_i - B_i%*%D_inv%*%C_i)
    B_i_inv   <- -A_i_inv%*%B_i%*%D_inv
    C_i_inv   <- -D_inv%*%C_i%*%A_i_inv
    D_i_inv   <- D_inv - C_i_inv%*%(B_i%*%D_inv)

    XtX_i_inv <- rbind(cbind(A_i_inv, B_i_inv), cbind(C_i_inv, D_i_inv))

    beta_i    <- XtX_i_inv%*%rbind(tX_i%*%y, tmp_y)
    mu_i      <- cbind(X_i, tmp)%*%beta_i

    res_i     <- y - mu_i
    sigma2_i  <- mean(res_i^2)

    profile_i <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)

    index_para[i, ]    <- beta_i[1:(2*p + q + r)]
    index_para_sd[i, ] <- (sigma2_i*diag(XtX_i_inv))[1:(2*p + q + r)]^(1/2)
    index_mu[i, ]      <- mu_i
    index_profile[i]   <- profile_i

  }

  index_ln_max   <- max(index_profile)
  index_Dn       <- 2*(max(index_profile) - index_profile)
  index_aic      <- 2*(index_profile - (2*p + q + r + 1))

  index          <- as.numeric(which.min(index_Dn))
  index_index    <- index_val[index]

  index_beta     <- index_para[index, 1:(p + p)]
  index_betaL    <- index_para[index, 1:p]
  index_betaR    <- index_para[index, 1:p + p]
  index_gamma    <- NULL; if (q != 0) index_gamma <- index_para[index, 1:q + 2*p]
  index_dummy    <- NULL; if (r != 0) index_dummy <- index_para[index, 1:r + 2*p + q]

  index_beta_sd  <- index_para_sd[index, 1:(p + p)]
  index_betaL_sd <- index_para_sd[index, 1:p]
  index_betaR_sd <- index_para_sd[index, 1:p + p]
  index_gamma_sd <- NULL; if (q != 0) index_gamma_sd <- index_para_sd[index, 1:q + 2*p]
  index_dummy_sd <- NULL; if (r != 0) index_dummy_sd <- index_para_sd[index, 1:r + 2*p + q]


  # Compute mu, sigma and res for the minimum deviance

  index_mu       <- index_mu[index, ]
  index_res      <- y - index_mu
  index_sigma    <- mean(index_res^2)^(1/2)

  invisible(list(ln_max       = index_ln_max,
                 aic          = index_aic,
                 Dn           = index_Dn,
                 beta         = index_beta,
                 betaL        = index_betaL,
                 betaR        = index_betaR,
                 tau          = index_index,
                 gamma        = index_gamma,
                 dummy        = index_dummy,
                 mu           = index_mu,
                 sigma        = index_sigma,
                 sigma_beta   = index_beta_sd,
                 sigma_betaL  = index_betaL_sd,
                 sigma_betaR  = index_betaR_sd,
                 sigma_gamma  = index_gamma_sd,
                 sigma_dummy  = index_dummy_sd,
                 res          = index_res,
                 index_val    = index_val))

}

#' Estimating parameters for the Panale Regression Model with a Change Point
#'
#' This function estimate the parameters, including the location of the change point,
#' for the model specified by `cdcp.regression.data(...)`.
#'
#' @param data output list of data from `cdcp.regression.data(...)`.
#' @param index_val vector of indexes for the locations for a potential change points. 
#' @return A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales.
#' A list of estimated parameters, the maximised log-likelihood, aic, deviance and residuales...
#' An object from cdcp.regression.estimation(...) is a list containing the following components:
#' \describe{
#'   \item{ln_max}{maximised log-likelihood for each potential split in index_val}
#'   \item{aic}{Akaike Information Criterion (AIC) for each potential split in index_val}
#'   \item{Dn}{deviance for each potential split in index_val}
#'   \item{beta}{estimated beta (left and right) for the optimal change point}
#'   \item{betaL}{estimated beta to the left of the optimal change point}
#'   \item{betaR}{estimated beta to the right of the optimal change point}
#'   \item{tau}{optimal chage point}
#'   \item{gamma}{estimated gamma parameter}
#'   \item{dummy}{estimated dummy paramters}
#'   \item{mu}{estimated expecation}
#'   \item{sigma}{estimated sigma}
#'   \item{sigma_beta}{estimated sigma for estimated beta}
#'   \item{sigma_betaL}{estimated sigma for estimated beta to the left for optimal change point}
#'   \item{sigma_betaR}{estimated sigma for estimated beta to the right for optimal change point}
#'   \item{sigma_gamma}{estimated sigma for estimated gamma}
#'   \item{sigma_dummy}{estimated sigma for estimated dummy paramters}
#'   \item{res}{residuales}
#'   \item{index_val}{same as index_val argument}
#' }
#'
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' 
#' # Example 1: One individual 
#' n <- 100
#' x <- 1:n
#' X <- cbind(1, x)
#' y <- 3*(x <= 50) + rnorm(n)
#' data <- cdcp.regression.data(y = y, X = X, index = x, group = rep(1, n))
#' index_val <- cdcp.regression.data.find.index.val(data)$index_val
#' fit <- cdcp.regression.estimate(data, index_val)
#' 
#' # Example 2: 
#' data <- cdcp.regression.data.sim()$data
#' index_val <- cdcp.regression.data.find.index.val(data)$index_val
#' fit <- cdcp.regression.estimate(data, index_val)
#' @export
cdcp.regression.estimate <- function(data, index_val) {

  if (!all(c("y", "X", "Z", "index", "group") %in% attributes(data)$names)) {
    stop("ERROR: data is not a proper cdcp regression data object. \n")
  }

  if (!all(index_val %in% unique(data$index))) {
    stop("ERROR: Values in index_val are not found in data index.\n")
  }

  return(cdcp.regression.estimate.function(data = data, index_val = index_val))
}

cdcd.pannel.estimation.test <- function() {

  set.seed(441)

  country_names <- c("Mirasius", "Oras", "Anglinthius", "Olvion")
  country_names <- sort(country_names)


  year_min <- 1900
  year_max <- 2000
  year_seq <- year_min:year_max


  n           <- length(year_seq)
  m           <- length(country_names)

  beta_L      <- 0.35
  beta_R      <- 0.45

  tau_0       <- 1950
  gamma_0     <- 1

  u           <- rep(year_seq, times = m)/year_max
  mu          <- rep(beta_L + (beta_R - beta_L)*(year_seq > tau_0), times = m) + gamma_0*u
  sigma       <- 0.1

  y           <- mu + sigma*rnorm(n*m)
  X           <- as.matrix(1 + mu*0)
  Z           <- as.matrix(u)

  index       <- rep(year_seq, times = m)
  group       <- rep(country_names, each = n)
  country     <- rep(country_names, each = n)
  id          <- rep(1:m, each = n)

  index_val   <- 1910:1990

  data        <- cdcp.regression.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = FALSE, group_dummy = FALSE, lag = FALSE)
  model       <- cdcp.regression.estimate(data, index_val)

  # Standard Dummy
  # for (index_i in index_val[1:5]) {
  #
  #   cat("\nindex:", index_i, "\n")
  #
  #   #
  #   model_i <- cdcp.panel.estimate(data = data, index_val = index_i)
  #
  #
  #   # Using lm in R
  #   X_L <- X
  #   X_L[index <= index_i, ] <- 0
  #
  #   lm_i     <- lm(y ~ X_L + Z)
  #
  #   para_i   <- lm_i$coefficients
  #   res_i    <- lm_i$residuals
  #   sigma2_i <- mean(res_i^2)
  #
  #   loglik_i <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)
  #   loglik_i <-  logLik(lm_i)
  #
  #   betaL_i  <- para_i[1]
  #   betaR_i  <- para_i[1] + para_i[2]
  #
  #   gamma_i  <- para_i[3]
  #
  #   print(as.numeric(c(betaL_i, betaR_i, gamma_i)))
  #   print(c(model_i$beta, model_i$gamma))
  #
  #   cat("\n")
  #   print(sigma2_i^(1/2))
  #   print(model_i$sigma)
  #
  #   cat("\n")
  #   print(loglik_i[[1]])
  #   print(model_i$ln_max)
  #
  # }



  # # With group Dummy
  # data <- cdcp.panel.clean.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = FALSE, group_dummy = TRUE, lag = FALSE)
  # for (index_i in index_val[1:5]) {
  #
  #   cat("\nindex:", index_i, "\n")
  #
  #   #
  #   model_i <- cdcp.panel.estimate(data = data, index_val = index_i)
  #
  #
  #   # Using lm in R
  #   X_L <- X
  #   X_L[index <= index_i, ] <- 0
  #
  #   lm_i     <- lm(y ~ X_L + Z + as.factor(group))
  #
  #   para_i   <- lm_i$coefficients
  #   res_i    <- lm_i$residuals
  #   sigma2_i <- mean(res_i^2)
  #
  #   loglik_i <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)
  #   loglik_i <-  logLik(lm_i)
  #
  #   betaL_i  <- para_i[1]
  #   betaR_i  <- para_i[1] + para_i[2]
  #
  #   gamma_i  <- para_i[3]
  #
  #   print(as.numeric(c(betaL_i, betaR_i, gamma_i)))
  #   print(c(model_i$beta, model_i$gamma))
  #
  #   cat("\n")
  #   print(as.numeric(c(para_i[-(1:3)])))
  #   print(c(model_i$dummy))
  #
  #   cat("\n")
  #   print(sigma2_i^(1/2))
  #   print(model_i$sigma)
  #
  #   cat("\n")
  #   print(loglik_i[[1]])
  #   print(model_i$ln_max)
  #
  # }


  # With group Dummy and lag
  # data <- cdcp.panel.clean.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = FALSE, group_dummy = TRUE, lag = TRUE)
  # for (index_i in index_val[1:5]) {
  #
  #   cat("\nindex:", index_i, "\n")
  #
  #   #
  #   model_i <- cdcp.panel.estimate(data = data, index_val = index_i)
  #
  #
  #   # Using lm in R
  #   X_L <- X
  #   X_L[index <= index_i, ] <- 0
  #
  #   y_lag    <- matrix(y, ncol = m, nrow = n, byrow = FALSE)
  #   y_lag    <- as.vector(rbind(NA, y_lag[1:(n - 1), ]))
  #
  #   lm_i     <- lm(y ~ X_L + Z + y_lag + as.factor(group))
  #
  #   para_i   <- lm_i$coefficients
  #   res_i    <- lm_i$residuals
  #   sigma2_i <- mean(res_i^2)
  #
  #   loglik_i <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)
  #   loglik_i <-  logLik(lm_i)
  #
  #   betaL_i  <- para_i[1]
  #   betaR_i  <- para_i[1] + para_i[2]
  #
  #   rho_i    <- para_i[4]
  #
  #   gamma_i  <- para_i[3]
  #
  #   print(as.numeric(c(betaL_i, betaR_i, gamma_i, rho_i)))
  #   print(c(model_i$beta, model_i$gamma))
  #
  #   cat("\n")
  #   print(as.numeric(c(para_i[-(1:4)])))
  #   print(c(model_i$dummy))
  #
  #   cat("\n")
  #   print(sigma2_i^(1/2))
  #   print(model_i$sigma)
  #
  #   cat("\n")
  #   print(loglik_i[[1]])
  #   print(model_i$ln_max)
  #
  # }


  # With year and group Dummy and lag
  data <- cdcp.panel.clean.data(y = y, X = X, Z = NULL, index = index, group = group, index_dummy = TRUE, group_dummy = TRUE, lag = TRUE, intercept = TRUE)
  for (index_i in index_val[1:5]) {

    cat("\nindex:", index_i, "\n")

    #
    model_i <- cdcp.panel.estimate(data = data, index_val = index_i)


    # Using lm in R
    X_L <- X
    X_L[index <= index_i, ] <- 0

    y_lag    <- matrix(y, ncol = m, nrow = n, byrow = FALSE)
    y_lag    <- as.vector(rbind(NA, y_lag[1:(n - 1), ]))

    lm_i     <- lm(y ~ X_L + y_lag + as.factor(group) + as.factor(index))

    para_i   <- lm_i$coefficients
    res_i    <- lm_i$residuals
    sigma2_i <- mean(res_i^2)

    loglik_i <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)
    loglik_i <-  logLik(lm_i)

    betaL_i  <- para_i[1]
    betaR_i  <- para_i[1] + para_i[2]

    rho_i    <- para_i[3]

    print(as.numeric(c(betaL_i, betaR_i, rho_i)))
    print(c(model_i$beta, model_i$gamma))

    cat("\n")
    print(as.numeric(c(para_i[-(1:3)])))
    print(c(model_i$dummy))

    cat("\n")
    print(sigma2_i^(1/2))
    print(model_i$sigma)

    cat("\n")
    print(loglik_i[[1]])
    print(model_i$ln_max)

  }



}

# cdcd.pannel.estimation.test()


















