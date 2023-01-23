# Title: CDCPRegression - Confidence Curve for the Degree of Change (DoC)
#
# Summary: Functions needed to compute the confidence curve (cc) for the so-called degree of change 
# for a specific beta coefficient. 
#
# TODO:

cdcp.regression.beta.doc.function <- function(data, model, k, delta_val, index_val, boot = 0, boot_type = c("gaussian", "independent", "group"), cores = 4) {
  
  doMC::registerDoMC(cores)
  
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

  return(list(cc = diff_cc, cc_approx = pchisq(diff_Dn, df = 1), Dn = diff_Dn, ln_profile = diff_profile, delta_val = delta_val, Kn = Kn_sim))
}

cdcp.regression.beta.doc.check <- function(cc) {
 
  if (!all(c("cc", "cc_approx", "delta_val") %in% attributes(cc)$names)) {
    stop("ERROR: Provide a proper list of cc values, it must contain cc, cc_approx and delta_val. \n")
  }
  
  if (!(is.numeric(cc$delta_val) & length(cc$delta_val) > 0)) {
    stop("ERROR: delta_val must be an numeric vector. \n")
  }
  
  if (!is.null(cc$cc)) {
    if (!(is.numeric(cc$cc) & length(cc$cc) == length(cc$delta_val))) {
      stop("ERROR: cc must be an numeric vector and of the same length as delta_val. \n")
    }
  }
  
  if (!is.null(cc$cc_approx)) {
    if (!(is.numeric(cc$cc_approx) & length(cc$cc_approx) == length(cc$delta_val))) {
      stop("ERROR: cc_approx must be an numeric vector and of the same length as delta_val. \n")
    }
  }
    
}

#' Confidence Curve for the Degree of Change
#'
#' This function compute the confidence curve for the so-called degree of change for one important 
#' \eqn{\beta_k} coefficients in the model. The degree of change for \eqn{\beta_k} is defined as \eqn{\delta_k = \beta_{k, {\rm Left}} - \beta_{k, {\rm Right}}}; see Hermansen (2021) for details.
#' Note that the full bootstrap computation is quite computer intensive, it is therefore recommended 
#' to run `boot = 0` first and then use `boot = 1000` on a coarse set to either verify or compute the 
#' full confidence curve on a subset of values in `delta_val`.
#'
#' @param data the output from the `cdcp.regression.data(...)` function.
#' @param model the output of the `cdcp.regression.estimation(...)` function.
#' @param k the index of the \eqn{\beta} to compute the confidence curve for the degree of change, e.g. `k = 1` computes the confidence curve for the degree of change for the intercept (i.e. \eqn{\beta_0}). 
#' @param delta_val a consecutive sequence of values representing potential values of change.
#' @param boot the number of bootstrap samples. Note that if `boot = 0` use instead large sample approximations to compute the confidence curve. 
#' @param boot_type boot_type The method used to sample residuals for the bootstrap: `gaussian` use use parametric 
#' bootstrap under the assumption of a Gaussian model, `independent` makes independent draws from the 
#' estimated residuals and `group` use the assumed group structure to draw independent residuals. 
#' @return A list containing the following:
#' \describe{
#'   \item{cc}{the confidence curve computed for all values of delta_val}
#'   \item{cc_approx}{the approximative confidence curve computed for all values of delta_val}
#'   \item{Dn}{the deviance for all value of delta_val}
#'   \item{ln_profile}{the maximised profile log-likelihood for all delta_val}
#'   \item{delta_val}{same as `delta_val` argument}
#'   \item{Kn}{the empirical cdf of the deviance evaluated at the observed deviance}
#' }
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' # Example 1: 
#' data <- cdcp.regression.data.sim()
#' model <- cdcp.regression.estimate(data$data, data$index_val)
#' cc <- cdcp.regression.beta.doc(data$data, model, k = 1, seq(-0.20, 0.15, length.out = 100), data$index_val)
#' @export
#' @import foreach
#' @import doMC
cdcp.regression.beta.doc <- function(data, model, index_val, k, delta_val, boot = 0, boot_type = c("gaussian", "independent", "group"), cores = 4) {

  cdcp.regression.data.check.data(data)

  cdcp.regression.estimate.check.model(model)

  if (ncol(data$X) < k | (round(k) - k) > 0) {
    stop("ERROR: k must represent a column in X. \n")
  }
    
  if (!(is.numeric(delta_val) & length(delta_val) > 0)) {
    stop("ERROR: delta_val must be an numeric vector. \n")
  }
  
  if (!all(index_val %in% unique(data$index))) {
    stop("ERROR: Values in index_val are not found in data index.\n")
  }
  
  if (!(is.numeric(boot) & (round(boot) - boot) == 0)) {
    stop("ERROR: boot must be an integer. \n")
  }

  if (!(boot_type[1] %in% c("gaussian", "independent", "group"))) {
    stop("ERROR: Incorrect boot type, should be either gaussian, independent or group.\n")
  }

  return(cdcp.regression.beta.doc.function(data, model, k, delta_val, index_val, boot = boot, boot_type, cores = cores))

}

#' Plot Confidence Curve for the Degree of Change
#'
#' Plot the confidence curve (cc) for the degree of change for a specific coefficient, see 
#' `cdcp.regression.beta.doc(...)` for details. The plotting function is quite minimal, however, it should be easy to modify this to individual needs.
#' 
#' @param cc the output from running the `cdcp.regression.data(...)` function.
#' @param approx should the approximation of the confidence curve be plotted (`approx = TRUE` will plot this).
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' # Example 1: 
#' data <- cdcp.regression.data.sim()
#' model <- cdcp.regression.estimate(data$data, data$index_val)
#' cc <- cdcp.regression.beta.doc(data$data, model, k = 1, seq(-0.20, 0.15, length.out = 100), data$index_val)
#' cdcp.regression.beta.doc.plot(cc)
#' @export
cdcp.regression.beta.doc.plot <- function(cc, approx = TRUE) {
  
  cdcp.regression.beta.doc.check(cc)
  
  if (!is.logical(approx)) {
    stop("ERROR: approx must be logical.\n")
  }
  
  if (!is.null(cc$cc)) {
    plot(cc$delta_val, cc$cc, type = 'l', lwd = 1.3, ylim = c(0, 1), xlab = "", ylab = "")
    title(ylab = "Confidence Curve", xlab = "delta", line = 2.5)
  } else {
    cat("WARNING: The confidence curve is missing.\n")
  }
  
  if (approx) {
    if (!is.null(cc$cc_approx)) {
      plot(cc$delta_val, cc$cc_approx, type = 'l', lwd = 1.3, ylim = c(0, 1), xlab = "", ylab = "")
      title(ylab = "Confidence Curve (approximation)", xlab = "delta", line = 2.5)
    } else {
      cat("WARNING: The approximative confidence curve is missing.\n")
    }
  }
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























