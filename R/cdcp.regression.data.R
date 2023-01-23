# Title: CDCPRegression - Data
#
# Summary: Collection of functions to pre-process data for the CDCPRegression package. This is needed  
# to sure a structured format for efficiently and correctly computations. 
#
# TODO:

cdcp.regression.data.add.dummy <- function(D) {

  # TODO: 
  # 1) Verfiy the output matrix

  n <- dim(D)[1]
  p <- dim(D)[2]

  if (p == 1) D <- as.matrix(D)

  D_val_return <- NULL
  D_return <- NULL

  for (i in 1:p) {

    D_i_val <- unique(D[, i])
    m       <- length(D_i_val)

    if (m > 1) {
      D_i <- matrix(0, nrow = n, ncol = m - 1)

      for (j in 2:m) {
        D_i[D[, i] == D_i_val[j], j - 1] <- 1
      }
    }

    colnames(D_i) <- D_i_val[-1]

    D_val_return <- c(D_val_return, D_i_val[-1])
    D_return <- cbind(D_return, D_i)
  }

  return(list(D = D_return, D_val = D_val_return))

}

cdcp.regression.data.add.group.dummy <- function(group) {

  n     <- length(group)
  D_val <- unique(group)
  m     <- length(D_val)

  if (m > 1) {
    D <- matrix(0, nrow = n, ncol = m - 1)

    for (i in 2:m) {
      D[group == D_val[i], i - 1] <- 1
    }
  }

  colnames(D) <- D_val[-1]

  return(list(D = D, D_val = D_val[-1]))
}

cdcp.regression.data.add.index.dummy.with.cp.intercept <- function(index) {

  n     <- length(index)
  D_val <- sort(unique(index))
  m     <- length(D_val)
  D     <- matrix(0, nrow = n, ncol = m - 2)

  # Note: This is important. Since we have a year dummy in
  # the model, and a change point, and since intercept is
  # one of the parameters that may change, we can only have
  # m - 2 year coeffecients

  for (i in 2:(m - 1)) {
    D[index == D_val[i], i - 1] <- 1
  }

  colnames(D) <- D_val[-c(1, m)]

  return(list(D = D, D_val = D_val[-c(1, m)]))
}

cdcp.regression.data.add.index.dummy.with.intercept <- function(index) {

  n     <- length(index)
  D_val <- sort(unique(index))
  m     <- length(D_val)
  D     <- matrix(0, nrow = n, ncol = m - 1)

  # Note: Remove one dummy since there is a intercept in Z

  for (i in 2:m) {
    D[index == D_val[i], i - 1] <- 1
  }

  colnames(D) <- D_val[-1]

  return(list(D = D, D_val = D_val[-1]))



}

cdcp.regression.data.add.index.dummy <- function(index) {

  n     <- length(index)
  D_val <- sort(unique(index))
  m     <- length(D_val)
  D     <- matrix(0, nrow = n, ncol = m)

  for (i in 1:m) {
    D[index == D_val[i], i] <- 1
  }

  colnames(D) <- D_val

  return(list(D = D, D_val = D_val))
}

cdcp.regression.data.function <- function(y, X, Z = NULL, D = NULL, index, group, index_dummy = FALSE, group_dummy = FALSE, lag = FALSE) {

  # Add lagged response variable as protected covariate
  if (lag) {
    tmp     <- data.table::data.table(y = y, index = index, group = group)
    tmp_lag <- tmp[, y_lag := c(NA, y[-.N]), by = 'group']$y_lag
    Z       <- cbind(Z, y_lag = tmp_lag)
  }


  # Remove NA rows, this finds the longest continuous sequence without missing for each group
  index_na <- data.table::data.table(index = apply(cbind(y, X, Z, index, group), 1, anyNA), group)
  index_na[, tmp := (c(0, diff(index)) != 0) + 0, by = "group"]
  index_na[, tmp := cumsum(tmp), by = "group"]
  index_na[, tmp := sum(index == 0), by = c("group", "tmp")]
  index_na_unique <- index_na[, first(tmp), by = c("group", "index")][, sum(V1 == max(V1)), by = "group"]


  if (any(index_na_unique$V1 > 1)) {
    cat(paste("ERROR: Longest continuous sequence without missing is not unique for group: \n "))
    print(index_na_unique[V1 > 1]$group)
    stop()
  }

  index_na[, index := (tmp == max(tmp)) & (index == 0), by = "group"]
  index_na <- index_na$index


  if (FALSE) {
    # Test:

    x <- c(1, 1, 0, 0, 0, 0, 1)
    g <- c(1, 1, 1, 1, 1, 1, 1)

    x <- c(x, c(0, 1, 0, 0, 0, 0, 1, 0))
    g <- c(g, c(2, 2, 2, 2, 2, 2, 2, 2))

    x <- c(x, c(0, 0, 0, 0, 1))
    g <- c(g, c(3, 3, 3, 3, 3))

    x <- c(x, c(1, 0, 0, 0, 0, 0))
    g <- c(g, c(4, 4, 4, 4, 4, 4))

    x <- c(x, c(1, 1, 0, 0, 1, 0, 0))
    g <- c(g, c(5, 5, 5, 5, 5, 5, 5))

    x <- c(x, c(1, 1, 1, 1, 1, 1, 1))
    g <- c(g, c(6, 6, 6, 6, 6, 6, 6))


    test <- data.table::data.table(y = 1, x = x, g = g)
    test[, tmp0 := (c(0, diff(x)) != 0) + 0, by = "g"]
    test[, tmp1 := cumsum(tmp0), by = "g"]
    test[, tmp2 := sum(x == 0), by = c("g", "tmp1")]

    test_unique <- test[, first(tmp2), by = c("g", "tmp1")][, sum(V1 == max(V1)), by = "g"]

    if (any(test_unique$V1 > 1)) {
      cat(paste("ERROR: Longest continuous sequence without missing is not unique for group: \n "))
      print(test_unique[V1 > 1]$g)
    }
    test[, tmp3 := ((tmp2 == max(tmp2)) & x == 0) + 0, by = "g"]
  }

  if (any(!index_na)) {
    y     <- y[index_na]
    X     <- as.matrix(X[index_na, ])
    if (!is.null(Z)) {Z <- as.matrix(Z[index_na, ])}
    index <- index[index_na]
    group <- group[index_na]
  }

  # TODO:
  # 1) Stop if data is of length 0.

  # Check for data gaps in the data
  if (length(index) > 0) {
    index_gap <- data.table::data.table(index, group)
    index_gap[,  tmp1 := c(0, diff(index)), by = "group"]
    index_gap[1:2, tmp1 := max(tmp1), by = "group"]
    index_gap[, tmp2 := cumsum(tmp1 > 1), by = "group"]
    index_gap[, tmp3 := length(index), by = c("group", "tmp2")]

    index_gap_unique <- index_gap[, first(tmp3), by = c("group", "tmp2")][, sum(V1 == max(V1)), by = "group"]

    if (any(index_gap_unique$V1 > 1)) {
      cat(paste("ERROR: Longest continuous sequence without missing is not unique for group: \n "))
      print(index_gap_unique[V1 > 1]$group)
    }

    index_gap[, tmp4 := (tmp3 == max(tmp3)), by = "group"]
    index_gap <- index_gap$tmp4


    if (any(!index_gap)) {
      y     <- y[index_gap]
      X     <- as.matrix(X[index_gap, ])
      if (!is.null(Z)) {Z <- as.matrix(Z[index_gap, ])}
      index <- index[index_gap]
      group <- group[index_gap]
    }
  }

  n    <- length(y)
  p    <- ncol(X)
  q    <- ncol(Z)

  data <- list(y = y, X = X, Z = Z, index = index, group = group)

  # Add dummy variables
  if (!is.null(D)) {
    Dummy       <- cdcp.regression.data.add.dummy(D)
    D           <- Dummy$D
    D_val       <- list(D = Dummy$D_val)
    D_val_index <- list(D = 1:length(Dummy$D_val))
    nr_dummy    <- length(Dummy$D_val)
  } else {
    D           <- NULL
    D_val       <- list()
    D_val_index <- list()
    nr_dummy    <- 0
  }


  if (index_dummy) {
    # If there is a change point in the intercept, we need to correct the number
    # of index dummies, similar if there is a intercept among the protected variables

    intercept_in_X <- any(apply(X, 2, function(x) {length(unique(x))}) == 1)
  
    intercept_in_Z <- FALSE
    try(intercept_in_Z <- any(apply(Z, 2, function(x) {length(unique(x))}) == 1), silent = TRUE)

    if (intercept_in_X & intercept_in_Z) {
      stop("ERROR: Intercept found among X and Z. Remove that column from either X or Z to avoid problems with identifiability with index dummies in the model.\n")
    } else if (intercept_in_X) {
      Dummy <- cdcp.regression.data.add.index.dummy.with.cp.intercept(index)
      cat("WARNING: Intercept found in X. Adjust index dummies to avoid problems with identifiability.\n")
    } else if (intercept_in_Z) {
      Dummy <- cdcp.regression.data.add.index.dummy.with.intercept(index)
      cat("WARNING: Intercept found in Z. Adjust index dummies to avoid problems with identifiability.\n")
    } else {
      Dummy <- cdcp.regression.data.add.index.dummy(index)
    }

    D                 <- cbind(D, Dummy$D)
    D_val$index       <- Dummy$D_val
    D_val_index$index <- nr_dummy + 1:length(Dummy$D_val)
    nr_dummy          <- nr_dummy +   length(Dummy$D_val)
  }

  if (group_dummy) {
    Dummy             <- cdcp.regression.data.add.group.dummy(group)
    D                 <- cbind(D, Dummy$D)
    D_val$group       <- Dummy$D_val
    D_val_index$group <- nr_dummy + 1:length(Dummy$D_val)
    nr_dummy          <- nr_dummy +   length(Dummy$D_val)
  }

  r <- ncol(D)

  data$D <- D
  data$D_val <- D_val
  data$D_val_index <- D_val_index

  data$n <- n
  data$p <- p
  data$q <- q
  data$r <- r

  invisible(data)
}

cdcp.regression.data.check.data <- function(data) {
  
  if (!all(c("y", "X", "Z", "index", "group") %in% attributes(data)$names)) {
    stop("ERROR: data is not a proper cdcp regression data object. \n")
  }

  y <- data$y
  X <- data$X
  Z <- data$Z
  index <- data$index
  group <- data$group
  D <- data$D
  
  cdcp.regression.data.check(y, X, Z, index, group, D)
}

cdcp.regression.data.check <- function(y, X, Z, index, group, D) {
  
  # ERRORS:
  if (!is.numeric(y)) {
    stop("ERROR: The denpendent/response variables in y has to be numeric.\n")
  }
  
  if (!is.numeric(X)) {
    stop("ERROR: The independent variables in X has to be numeric.\n")
  }
  
  if (!is.matrix(X)) {
    stop("ERROR: X must be matrix. \n")
  }
  
  if (length(y) != nrow(X)) {
    stop("ERROR: X and y have different number of observations. \n")
  }
  
  if (length(y) != length(index)) {
    stop("ERROR: index and y have different number of observations. \n")
  }
  
  if (length(y) != length(group)) {
    stop("ERROR: group and y have different number of observations. \n")
  }
  
  if (!is.null(Z)) {
    if (!is.numeric(Z)) {
      stop("ERROR: The independent variables in Z has to be numeric.\n")
    }
    if (!is.matrix(Z)) {
      stop("ERROR: Z must be matrix. \n")
    }
    if (length(y) != nrow(Z)) {
      stop("ERROR: Z and y have different number of observations. \n")
    }
  }
  
  if (!is.null(D)) {
    if (!is.matrix(D)) {
      stop("ERROR: D must be a matrix. \n")
    }
    if (length(y) != nrow(D)) {
      stop("ERROR: D and y have different number of observations. \n")
    }
  }
}

#' Find Suggested Indexes for a Change Point
#'
#' For a model specified by the `cdcp.regression.data(...)` function, it finds a set of suggested 
#' locations for a change point that ensures a given number of data points per estimated parameter is 
#' at least `threshold` to the left and right of each potential change point. This function is most useful 
#' for large panel data with several dummy variables.
#'
#' @param data output from the `cdcp.regression.data(...)` function.
#' @param threshold minimum number of data points per estimated parameter.
#' @param index_dummy if dummy variables are used in the model, set `index_dummy = TRUE` since this adds additional 
#' constraints to the number of data points; see `cdcp.regression.data` for some additional details.
#' @return A list containing the suggested change points and the number of data points per parameter 
#' on each side of the change point: 
#' \describe{
#'   \item{index_val}{suggested locations (indexes) for a change point}
#'   \item{nr_para_left}{number of data points per paramter to the left for the corresponding index}
#'   \item{nr_para_right}{number of data points per paramter to the right for the corresponding index}
#'}
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' # Example 1: One individual 
#' n <- 100
#' x <- 1:n
#' X <- cbind(1, x)
#' y <- 3*(x <= 50) + rnorm(n)
#' data <- cdcp.regression.data(y = y, X = X, index = x, group = rep(1, n))
#' index_val <- cdcp.regression.data.find.index.val(data)$index_val
#' 
#' # Example 2: Simulated Data 
#' data <- cdcp.regression.data.sim()$data
#' index_val <- cdcp.regression.data.find.index.val(data)$index_val
#'
#' @export
cdcp.regression.data.find.index.val <- function(data, index_dummy = FALSE, threshold = 10) {

  # Check data in data object
  cdcp.regression.data.check.data(data)
  
  # Check of parameters: 
  if (!is.logical(index_dummy)) {
    stop("ERROR: The index_dummy variable has to be logical.\n")
  }  
  
  if (!is.numeric(threshold)) {
    stop("ERROR: The threshold variable has to be numerical.\n")
  }  
  
  index <- data$index
  group <- data$group
  p     <- data$p
  q     <- data$q; if (is.null(q)) q <- 0
  r     <- data$r; if (is.null(r)) r <- 0

  index_val <- range(index)
  index_val <- (index_val[1] + 1):index_val[2]

  data_left  <- index_val*0
  para_left  <- index_val*0

  data_right <- index_val*0
  para_right <- index_val*0

  for (i in seq_along(index_val)) {
    data_left[i] <- sum(index <= index_val[i])
    para_left[i] <- sum(index_val <= index_val[i])*(index_dummy)
    para_left[i] <- para_left[i] + length(unique(group[index_val <= index_val[i]]))

    data_right[i] <- sum(index > index_val[i])
    para_right[i] <- sum(index_val > index_val[i])*(index_dummy)
    para_right[i] <- para_right[i] + length(unique(group[index_val <= index_val[i]]))

  }

  data_para_left <- data_left/(para_left + (2*p + q)/2)
  data_para_left_max <- max(data_para_left, na.rm = TRUE)

  data_para_right <- data_right/(para_right + (2*p + q)/2)
  data_para_right_max <- max(data_para_right, na.rm = TRUE)

  
  if ((data_para_left_max < threshold) | (data_para_right_max < threshold)) {
    threshold_max <- min(data_para_left_max, data_para_right_max)
    threshold <- max(floor(threshold_max), 0.5)
    cat("WARNING: It is not possible to obtain the desiered number of effective parameters.")
    cat("Maxium threshold is", round(threshold_max, 2), "threshold is set to", threshold, "\n")
  }

  index_val <- min(index_val[data_para_left > threshold]):max(index_val[data_para_right > threshold])

  return(list(index_val = index_val, nr_para_left = data_para_left, nr_para_right = data_para_right))

}

#' Pre-process Data for the CDCPRegression Package
#'
#' This function pre-process and creates a structured set of data to be used with the various functions in the CDCPRegression
#' package. The model must be within the class \eqn{y_{i, j} = x_{i, j}^{\rm t} \beta_{\rm L} + x_{i, j}^{\rm t} (\beta_{\rm R} - \beta_{\rm L}) I(t_i < \delta) + \gamma Z_{i ,j} + \alpha D_{i, j} + y_{i - 1, j} + \epsilon_{i, j}}; see Hermansen (2021) for more details.
#'
#' @param y the response variable, a \eqn{(n*m)} numeric vector, where \eqn{m} is the number of individuals and \eqn{n} is the number of observations per individual.
#' @param X a numeric covariate matrix of size \eqn{(n*m) \times p}. The matrix X repent the part of the model specification we believe is affected by a change point. 
#' @param Z protected covariates, i.e. covariates we believe are constant and are not affected by the change point. 
#' @param D protected categorical covariates.
#' @param index a \eqn{(n*m)}-vector of indexes, e.g. time or year, such that index = 1 represent the first observation for each individual in the panel.
#' @param group a vector of length \eqn{(n*m)} representing a potential group structure in the data. If all observations belong to the same group, this should be a \eqn{(n*m)}-vector with the value 1.
#' @param index_dummy set `index_dummy = TRUE` to include a categorical variable for each index (e.g. year dummy) in the model.
#' @param group_dummy set `group_dummy = TRUE` to include a categorical variable for each group. 
#' @param lag set `lag = TRUE` to add the lagged response \eqn{y_{i - 1}} to the model.
#' @return A list of data used as input to the various functions in the CDCPRegression Package. 
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' # Example 1: Simple Illustration with one Individual 
#' 
#' index <- 1:n
#' group <- rep(1, n)
#'  
#' X <- matrix(1, nrow = n, ncol = 1)
#' y <- 1.5*(index <= 50) + rnorm(n)
#' 
#' data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
#' 
#'
#' # Example 2: Simple Illustration with m Individuals 
#' 
#' m <- 2
#' n <- 100 
#' 
#' index <- rep(1:n, times = m)
#' group <- rep(1:m, each = n)
#' 
#' X <- matrix(1, nrow = n*m)
#' y <- 1.5*(index <= 50) + rnorm(n = m*n)
#' 
#' data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
#' @export
#' @import data.table
cdcp.regression.data <- function(y, X, Z = NULL, D = NULL, index, group, index_dummy = FALSE, group_dummy = FALSE, lag = FALSE) {

  # Check of parameters: 
  if (!is.logical(index_dummy)) {
    stop("ERROR: The index_dummy variable has to be logical.\n")
  }
  
  if (!is.logical(group_dummy)) {
    stop("ERROR: The group_dummy variable has to be logical.\n")
  }
  
  if (!is.logical(lag)) {
    stop("ERROR: The lag variable has to be logical.\n")
  }
  
  
  # y: numeric of length n
  # X: numeric matrix of dimension n x p
  # Z: numeric matrix of dimension n x q
  # D: matrix of dimension n x r
  # index: vector of length n (e.g. year, weeks, days or time within each group)
  # group: vector of length n (e.g. country or region such that all rows where group is of type A includes all data for group A)

  # Check of data: 
  cdcp.regression.data.check(y, X, Z, index, group, D)
  
  # if (!is.numeric(y)) {
  #   stop("ERROR: The denpendent/response variables in y has to be numeric.\n")
  # }
  # 
  # if (!is.numeric(X)) {
  #   stop("ERROR: The independent variables in X has to be numeric.\n")
  # }
  # 
  # if (!is.matrix(X)) {
  #   stop("ERROR: X must be matrix. \n")
  # }
  # 
  # if (length(y) != nrow(X)) {
  #   stop("ERROR: X and y have different number of observations. \n")
  # }
  # 
  # if (length(y) != length(index)) {
  #   stop("ERROR: index and y have different number of observations. \n")
  # }
  # 
  # if (length(y) != length(group)) {
  #   stop("ERROR: group and y have different number of observations. \n")
  # }
  # 
  # if (!is.null(Z)) {
  #   if (!is.numeric(Z)) {
  #     stop("ERROR: The independent variables in Z has to be numeric.\n")
  #   }
  #   if (!is.matrix(Z)) {
  #     stop("ERROR: Z must be matrix. \n")
  #   }
  #   if (length(y) != nrow(Z)) {
  #     stop("ERROR: Z and y have different number of observations. \n")
  #   }
  # }
  # 
  # if (!is.null(D)) {
  #   if (!is.matrix(D)) {
  #     stop("ERROR: D must be a matrix. \n")
  #   }
  #   if (length(y) != nrow(D)) {
  #     stop("ERROR: D and y have different number of observations. \n")
  #   }
  # }

  # WARNINGS:
  if (!is.null(D)) {
    cat("WARNING: It is not recommended to have time related covariates in D and also in X or Z, at the same time.")
    cat(" This can lead to numerical instability and problems of identifiability in the model. \n")
  }
  
  
  # WARNINGS:
  if (index_dummy) {
    cat("WARNING: If X, Z or D contains time-dependent covariates, the use of index dummies are not recommended.")
    cat(" This can lead to numerical instability and problems of identifiability in the model. \n")
  }


  return(cdcp.regression.data.function(y, X, Z, D, index, group, index_dummy, group_dummy, lag))
}


#' Simulate Test Data for the CDCPRegression Package
#' 
#' Simulate data from the model \eqn{y_{i, j} = \beta_{\rm L} + (\beta_{\rm R} - \beta_{\rm L}) I(t_i > 1950) + \gamma z_i + \epsilon_{i, j}}, for  
#' \eqn{t_i = 1900, \ldots , 2000}, \eqn{z_i = i/n}, \eqn{j = 1, \ldots m}, \eqn{n = 100} and \eqn{m = 4}, and with independent \eqn{\epsilon_{i, j} \sim {\rm N}(0, \sigma^2)}. The simulated data set is then aggregate this using `cdcp.regression.data(...)` function to create a toy data 
#' set that can be used to test the functionality of the CDCPRegression package. 
#' 
#' @param seed a seed for the random number generator. 
#' @param beta the intercept to the left and right of the change point.    
#' @param sigma the standard deviation used in simulating data.
#' @return Output of `cdcp.regression.data(...)` for the simulated data. 
#' @references Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504
#' @references Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
#' @examples
#' data <- cdcp.regression.data.sim()
#' @export
cdcp.regression.data.sim <- function(seed = 1, beta = c(0.30, 0.35), sigma = 0.1) {
  
  if (!is.numeric(beta)) {
    stop("ERROR: beta must be numeric.\n")
  }
  
  if (!is.numeric(sigma)) {
    stop("ERROR: sigma must be numeric.\n")
  }
  
  if (length(beta) != 2) {
    stop("ERROR: beta must be of length 2.\n")
  }
  
  set.seed(seed)
  
  country_names <- c("Mirasius", "Oras", "Anglinthius", "Olvion")
  country_names <- sort(country_names)
  
  year_min <- 1900
  year_max <- 2000
  year_seq <- year_min:year_max
  
  n <- length(year_seq)
  m <- length(country_names)
  
  beta_L <- beta[1]
  beta_R <- beta[2]
  
  tau_0 <- 1950
  gamma_0 <- 1
  
  u <- rep(year_seq, times = m)/year_max
  mu <- rep(beta_L + (beta_R - beta_L)*(year_seq > tau_0), times = m) + gamma_0*u
  sigma <- sigma
  
  y <- mu + sigma*rnorm(n*m)
  X <- as.matrix(1 + mu*0)
  Z <- as.matrix(u)
  
  index <- rep(year_seq, times = m)
  group <- rep(country_names, each = n)
  country <- rep(country_names, each = n)
  id <- rep(1:m, each = n)
  
  index_val <- 1910:1990
  
  data <- cdcp.regression.data(y = y, X = X, Z = Z, index = index, group = group, index_dummy = FALSE, group_dummy = FALSE, lag = FALSE)
  
  return(list(data = data, index_val = index_val))
}



