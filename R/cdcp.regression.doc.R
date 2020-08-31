# Title: Confidence Curve for the Degree of Change (DoC)
#
# Summary:
#
# TODO:
#
cdcp.regression.beta.doc.function <- function(data, model, k, delta_val, index_val, boot = 0, boot_type = c("gaussian", "independent", "group")) {

  y     <- data$y
  n     <- length(y)
  X     <- data$X

  p     <- data$p
  Z     <- data$Z
  p     <- data$p
  q     <- data$q; if (is.null(q)) {q <- 0}
  r     <- data$r; if (is.null(r)) {r <- 0}
  l     <- 2*p + q + r - 1


  Dummy <- data$D

  index <- data$index

  # From the model
  res    <- model$res
  sigma  <- model$sigma
  beta_L <- model$betaL
  beta_R <- model$betaR
  tau    <- model$tau

  group           <- data$group
  group_val       <- unique(group)


  # Old:
  # 1)
  # tmp   <- cbind(Z, Dummy)
  # A     <- NULL
  # B     <- tmp
  # C     <- t(tmp)
  # D     <- t(tmp)%*%tmp
  # D_inv <- solve(D)

  # New:
  # 2)
  tmp   <- cbind(Z, Dummy); if (is.null(tmp)) {tmp <- matrix(0, nrow = n, ncol = 1)}
  A     <- NULL
  B     <- tmp
  C     <- t(tmp)
  D     <- crossprod(tmp)
  if (length(D) == 1 & D[1, 1] == 0) {D_inv <- 0} else {D_inv <- solve(D)}

  X_0   <- X*0



  # Alternatives for simulating the residuals for the bootstrap
  if (boot_type[1] == "gaussian") {
    epsilon <- sigma*rnorm(n = boot*n)
    epsilon <- matrix(epsilon, nrow = n, ncol = boot)
  } else if (boot_type[1] == "independent") {
    epsilon <- matrix(sample(x = res, size = boot*n, replace = TRUE), nrow = n, ncol = boot)
  } else if (boot_type[1] == "group") {
    epsilon <- matrix(NA, nrow = n, ncol = boot, byrow = FALSE)
    for (group_i in group_val) {
      index_i            <- group == group_i
      n_i                <- sum(index_i)
      res_group_i        <- res[index_i]
      epsilon[index_i, ] <- sample(res_group_i, size = boot*n_i, replace = TRUE)
    }
  }
  # ---------------------------------------------------------------------------------------




  if (boot > 0) {


    pb <- txtProgressBar(min = 0, max = length(delta_val), style = 3)


    Kn_sim_delta_val <- foreach::foreach(i = seq_along(delta_val)) %do% {

      delta_i <- delta_val[i]

      X_L            <- X_0
      index_L        <- index <= tau
      X_L[index_L, ] <- X[index_L, ]

      X_R            <- X_0
      index_R        <- index > tau
      X_R[index_R, ] <- X[index > tau, ]

      mu_i <- X_L*beta_L + X_R*(beta_L - delta_i)

      Kn_sim_delta <- foreach::foreach(b = 1:boot, .combine = rbind) %dopar% {

        y_b <- as.vector(mu_i + epsilon[, b])

        profile_delta <- foreach::foreach(j = seq_along(index_val), .combine = cbind, index_val) %do% {

          X_L            <- X_0
          index_L        <- index <= index_val[j]
          X_L[index_L, ] <- X[index_L, ]

          X_R            <- X_0
          index_R        <- index > index_val[j]
          X_R[index_R, ] <- X[index_R, ]

          X_j            <- cbind(X_L, X_R)


          # for the next loop...
          X_i          <- X_j

          X_j[, p + k] <- X[, k]
          X_j          <- X_j[, -k]
          tX_j         <- t(X_j)

          A_i          <- tX_j%*%X_j
          B_i          <- tX_j%*%B
          C_i          <- C%*%X_j

          B_i_D_inv <- B_i%*%D_inv

          A_i_inv   <-  solve(A_i - B_i_D_inv%*%C_i)
          B_i_inv   <- -A_i_inv%*%(B_i_D_inv)
          C_i_inv   <- -(D_inv%*%C_i)%*%A_i_inv
          D_i_inv   <-  D_inv - C_i_inv%*%(B_i_D_inv)

          XtX_i_inv <- rbind(cbind(A_i_inv, B_i_inv), cbind(C_i_inv, D_i_inv))

          tmp_y_a   <- t(tmp)%*%y_b
          tmp_y_b   <- t(tmp)%*%X_i[, k]

          # New code
          y_i       <- y_b - X_i[, k]%*%t(delta_val)
          tmp_y_i   <- c(tmp_y_a) - c(tmp_y_b)%*%t(delta_val)

          res_i     <- y_i - cbind(X_j, tmp)%*%(XtX_i_inv%*%rbind(tX_j%*%y_i, tmp_y_i))
          sigma2_i  <- colMeans(res_i^2)

          profile_j <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)

          return(profile_j)
        }

        profile_delta <- as.numeric(apply(profile_delta, 1, max))

        Kn_sim_i <- 2*(max(profile_delta) - profile_delta)[delta_i == delta_val]

        return(c(Kn_sim_i, delta_i))
      }

      setTxtProgressBar(pb, i)

      return(data.table(Kn = Kn_sim_delta[, 1], delta = Kn_sim_delta[, 2]))
    }
    close(pb)

    Kn_sim <- rbindlist(Kn_sim_delta_val)
  } else {
    Kn_sim <- NULL
  }


  profile_delta <- foreach::foreach(j = seq_along(index_val), .combine = cbind, index_val) %do% {

    X_L            <- X_0
    index_L        <- index <= index_val[j]
    X_L[index_L, ] <- X[index_L, ]

    X_R            <- X_0
    index_R        <- index > index_val[j]
    X_R[index_R, ] <- X[index_R, ]

    X_j            <- cbind(X_L, X_R)


    # for the next loop...
    X_i          <- X_j

    X_j[, p + k] <- X[, k]
    X_j          <- X_j[, -k]
    tX_j         <- t(X_j)


    A_i       <- tX_j%*%X_j
    B_i       <- tX_j%*%B
    C_i       <- C%*%X_j

    B_i_D_inv <- B_i%*%D_inv

    A_i_inv   <- solve(A_i - B_i_D_inv%*%C_i)
    B_i_inv   <- -A_i_inv%*%(B_i_D_inv)
    C_i_inv   <- -(D_inv%*%C_i)%*%A_i_inv
    D_i_inv   <-  D_inv - C_i_inv%*%(B_i_D_inv)

    XtX_i_inv <- rbind(cbind(A_i_inv, B_i_inv), cbind(C_i_inv, D_i_inv))

    tmp_y_a   <- t(tmp)%*%y
    tmp_y_b   <- t(tmp)%*%X_i[, k]

    y_i       <- y - X_i[, k]%*%t(delta_val)

    tmp_y_i   <- c(tmp_y_a) - c(tmp_y_b)%*%t(delta_val)

    res_i     <- y_i - cbind(X_j, tmp)%*%(XtX_i_inv%*%rbind(tX_j%*%y_i, tmp_y_i))
    sigma2_i  <- colMeans(res_i^2)

    profile_j <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)

    return(profile_j)
  }

  profile_delta <- apply(profile_delta, 1, max)



  # Approximation:
  # diff_profile <- profile_delta
  # diff_Dn      <- 2*(max(diff_profile) - diff_profile)
  # diff_cc      <- pchisq(diff_Dn, df = 1)

  # Bootstrap:
  diff_profile <- profile_delta
  diff_Dn      <- 2*(max(diff_profile) - diff_profile)

  if (boot > 0) {
    diff_cc <- delta_val + NA
    for (j in seq_along(delta_val)) {
      Kn_j       <- ecdf(Kn_sim[delta == delta_val[j]]$Kn)
      diff_cc[j] <- Kn_j(diff_Dn[j])
    }
  } else {
    diff_cc <- NULL
  }

  invisible(list(cc = diff_cc, cc_approx = pchisq(diff_Dn, df = 1), Dn = diff_Dn, ln_profile = diff_profile, delta_val = delta_val, Kn = Kn_sim))
}

#' Confidence Distribution for the Degree of Change
#'
#' This function can compute both the bootstrap and the approximation...
#'
#' @param data A data object/list that specify the structure of the model created by the `cdcp.regression.data(...)`.
#' @param model The estimated model object from `cdcp.regression.estimation(...)`.
#' @param k
#' @param index_val A consecutive sequence of indexes representing the locations for a potential change point.
#' @param boot The number of bootstrap samples.
#' @param boot_type A subset of `index_val` used to speed up the calculations.
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
#'
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @export
cdcp.regression.beta.doc <- function(data, model, k, delta_val, index_val, boot = 0, boot_type = c("gaussian", "independent", "group")) {
  # Check size of k against dimension of X

  if (!all(c("y", "X", "Z", "index", "group") %in% attributes(data)$names)) {
    stop("ERROR: data is not a proper cdcp regression data object. \n")
  }

  if (!all(index_val %in% unique(data$index))) {
    stop("ERROR: Values in index_val are not found in data index.\n")
  }

  if (ncol(data$X) < k) {
    stop("ERROR: k is not a column in X. \n")
  }

  if (!(boot_type[1] %in% c("gaussian", "independent", "group"))) {
    stop("ERROR: Incorrect boot type, should be either gaussian, independent or group.\n")
  }


  return(cdcp.regression.beta.doc.function(data, model, k, delta_val, index_val, boot = 0, boot_type))

}






if (FALSE) {

  data_raw         <- fread(input = '~/Dropbox/shared/Change points/Data/hvdem.csv')
  data_raw_country <- unique(data_raw$country_name)
  data_raw         <- data_raw[country_name %in% data_raw_country[34:35]]

  y         <- log(data_raw$v2x_polyarchy)
  X         <- cbind(data_raw$Maddison_gdppc_1990_estimate_lag)
  Z         <- cbind(data_raw$Maddison_pop_estimate)
  index     <- data_raw$year
  group     <- data_raw$country_name

  data      <- cdcp.regression.clean.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = TRUE, group_dummy = TRUE, lag = TRUE)
  index_val <- cdcp.regression.clean.data.find.index.val(data, threshold = 5)$index_val


  model     <- cdcp.regression.estimate(data = data, index_val = index_val)

  delta_val <- seq(0.05, 0.3, length.out = 50)
  doc       <- cdcp.regression.beta.doc.function(data = data, model = model, k = 1, delta_val = delta_val, index_val = index_val, B = 100)

  with(doc, plot(delta_val, cc_approx, type = 'l'))
  with(doc, lines(delta_val, cc, type = 'l', col = "green"))

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

  # print(system.time(foreach (i=1:5, .combine='c') %dopar% {Sys.sleep(2);i}))
  # print(system.time(foreach (i=1:5, .combine='c') %do% {Sys.sleep(2);i}))

}






# Old:
if (FALSE) {
  # Degree of change:
  #
  #
  cdcp.regression.beta.doc.function <- function(data, k, delta_val, index_val) {

    y     <- data$y
    n     <- length(y)
    X     <- data$X

    p     <- data$p
    Z     <- data$Z

    Dummy <- data$D

    index <- data$index

    # Old:
    # 1)
    # tmp   <- cbind(Z, Dummy)
    # A     <- NULL
    # B     <- tmp
    # C     <- t(tmp)
    # D     <- t(tmp)%*%tmp
    # D_inv <- solve(D)
    # New:
    # 2)
    tmp   <- cbind(Z, Dummy); if (is.null(tmp)) {tmp <- matrix(0, nrow = n, ncol = 1)}
    A     <- NULL
    B     <- tmp
    C     <- t(tmp)
    D     <- crossprod(tmp)
    D_inv <- 0; if (D != 0) {D_inv <- solve(D)}

    X_0   <- X*0

    if (TRUE) {
      profile_delta <- foreach::foreach(j = seq_along(index_val), .combine = cbind, index_val) %dopar% {

        X_L            <- X_0
        index_L        <- index <= index_val[j]
        X_L[index_L, ] <- X[index_L, ]

        X_R            <- X_0
        index_R        <- index > index_val[j]
        X_R[index_R, ] <- X[index_R, ]

        X_j            <- cbind(X_L, X_R)


        # for the next loop...
        X_i            <- X_j

        X_j[, p + k]   <- X[, k]
        X_j            <- X_j[, -k]
        tX_j           <- t(X_j)

        # Here: Add test
        #
        # if (is.null(tmp)) {
        #   A_i            <- tX_j%*%X_j
        #   B_i            <- NULL
        #   C_i            <- NULL
        #   B_i_D_inv      <- NULL
        #   A_i_inv        <- solve(A_i) #
        #   B_i_inv        <- NULL
        #   C_i_inv        <- NULL
        #   D_i_inv        <- NULL
        #
        #   tmp_y_a        <- NULL
        #   tmp_y_b        <- NULL
        # } else {
        #
        # }


        A_i       <- tX_j%*%X_j
        B_i       <- tX_j%*%B
        C_i       <- C%*%X_j

        B_i_D_inv <- B_i%*%D_inv

        A_i_inv   <- solve(A_i - B_i_D_inv%*%C_i)
        B_i_inv   <- -A_i_inv%*%(B_i_D_inv)
        C_i_inv   <- -(D_inv%*%C_i)%*%A_i_inv
        D_i_inv   <-  D_inv - C_i_inv%*%(B_i_D_inv)

        XtX_i_inv <- rbind(cbind(A_i_inv, B_i_inv), cbind(C_i_inv, D_i_inv))

        tmp_y_a   <- t(tmp)%*%y
        tmp_y_b   <- t(tmp)%*%X_i[, k]

        profile_j <- foreach::foreach(i = seq_along(delta_val), .combine = c) %do% {

          # This part split into two parts...
          y_i       <- y - X_i[, k]*delta_val[i]
          # tmp_y_i   <- t(tmp)%*%y_i

          tmp_y_i   <- tmp_y_a - tmp_y_b*delta_val[i]

          # tmp_y_i <- t(tmp)%*%(y - X_j[, k]*delta_val[i])
          # tmp_y_i <- t(tmp)%*%y - (t(tmp)%*%X_i[, k])*delta_val[i]
          # res_i   <- y_i - cbind(X_j, tmp)%*%(XtX_i_inv%*%rbind(tX_j%*%y_i, tmp_y_i))

          res_i     <- y_i - cbind(X_j, tmp)%*%(XtX_i_inv%*%rbind(tX_j%*%y_i, tmp_y_i))
          sigma2_i  <- mean(res_i^2)

          profile_i <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)

          return(profile_i)
        }
        return(profile_j)
      }
      profile_delta <- apply(profile_delta, 1, max)
    }

    # Can this be made more efficient?
    # tmp_y <- foreach::foreach(i = seq_along(delta_val), .combine = c) %dopar% {
    #   tmp_y_j <- foreach::foreach(j = seq_along(index_val), .combine = c) %do% {
    #     X_L            <- X*0
    #     index_L        <- index <= index_val[j]
    #     X_L[index_L, ] <- X[index_L, ]
    #
    #     X_R            <- X*0
    #     index_R        <- index > index_val[j]
    #     X_R[index_R, ] <- X[index_R, ]
    #     X_j            <- cbind(X_L, X_R)
    #
    #     y_i            <- y - X_j[, k]*delta_val[i]
    #
    #     tmp_y_i        <- t(tmp)%*%y_i
    #     return(tmp_y_i)
    #   }
    #   return(tmp_y_j)
    # }
    #


    if (FALSE) {
      profile_delta <- foreach::foreach(i = seq_along(delta_val), .combine = c) %dopar% {

        profile_i <- foreach::foreach(j = seq_along(index_val), .combine = c) %do% {

          X_L            <- X*0
          index_L        <- index <= index_val[j]
          X_L[index_L, ] <- X[index_L, ]

          X_R            <- X*0
          index_R        <- index > index_val[j]
          X_R[index_R, ] <- X[index_R, ]
          X_j            <- cbind(X_L, X_R)

          y_i            <- y - X_j[, k]*delta_val[i]

          # This takes time to compute ...
          tmp_y_i        <- t(tmp)%*%y_i
          # tmp_y_i        <- as.matrix(tmp_y[, i])

          X_j[, p + k]   <- X[, k]
          X_j            <- X_j[, -k]
          tX_j           <- t(X_j)

          A_i            <- tX_j%*%X_j
          B_i            <- tX_j%*%B
          C_i            <- C%*%X_j

          B_i_D_inv      <- B_i%*%D_inv

          A_i_inv    <- solve(A_i - B_i_D_inv%*%C_i)
          B_i_inv    <- -A_i_inv%*%(B_i_D_inv)
          C_i_inv    <- -(D_inv%*%C_i)%*%A_i_inv
          D_i_inv    <-  D_inv - C_i_inv%*%(B_i_D_inv)

          XtX_i_inv  <- rbind(cbind(A_i_inv, B_i_inv), cbind(C_i_inv, D_i_inv))

          res_i      <- y_i - cbind(X_j, tmp)%*%(XtX_i_inv%*%rbind(tX_j%*%y_i, tmp_y_i))
          sigma2_i   <- mean(res_i^2)

          profile_i  <- -(n/2)*(log(2*pi) + log(sigma2_i) + 1)

          # profile_i <- 0

          return(profile_i)
        }
        return(max(profile_i))
      }
    }

    diff_profile <- profile_delta
    diff_Dn      <- 2*(max(diff_profile) - diff_profile)
    diff_cc      <- pchisq(diff_Dn, df = 1)


    invisible(list(cc = diff_cc, Dn = diff_Dn, ln_profile = diff_profile, delta_val = delta_val))
  }

  cdcp.regression.beta.doc <- function(data, k, delta_val, index_val) {


  }



}






















