rm(list = ls())

library(data.table)
library(readstata13)
library(foreign)
library('doMC')
registerDoMC(4)

# HVDEM: Change point for each county
#
#


source("R/cdcp.panel.functions.R")

data_raw <- read.dta13(file = "~/Dropbox/shared/Change points/Data/HVDEM_Polyarchy_with_region_covariates_updated_update.dta")
data_raw <- as.data.table(data_raw)
data_raw[, s_far_Maddison_gdp_1990_estimate_lag := c(NA, s_far_Maddison_gdp_1990_estimate[-.N]), by = 'country_name']

# data_raw <- as.data.table(data_raw)
data_raw[, s_mil_region := mean(s_mil_region, na.rm = TRUE), by = "country_name"]


# SEA= Sotheast Asia,
# SA=South Asia,
# EA=East Asia (including pacific here),
# ASS= Africa sub-Sahara,
# LA= Latin America,
# MENA= Middle East North Africa,
# EE= East Europe and Post Soviet,
# WE= West Europe + .


data_tmp <-  matrix(c(2, "LA", "Suriname",
                      4, "ASS", "South Sudan",
                      1, "EE", "Kosovo",
                      4, "ASS", "Cape Verde",
                      8, "SA", "Maldives",
                      3, "MENA", "Palestine/West Bank",
                      3, "MENA", "Palestine/Gaza",
                      4, "ASS", "Somaliland",
                      2, "LA", "Barbados",
                      6, "EA", "Hong Kong",
                      5, "WE", "Iceland",
                      5, "WE", "Luxembourg",
                      5, "WE", "Malta",
                      4, "ASS", "Sao Tome and Principe",
                      4, "ASS", "Seychelles",
                      6, "EA", "Vanuatu",
                      3, "MENA", "Palestine/British Mandate",
                      4, "ASS", "Zanzibar",
                      5, "WE", "Hanover",
                      5, "WE", "Hesse-Kassel",
                      5, "WE", "Hesse-Darmstadt",
                      5, "WE", "Mecklenburg Schwerin",
                      5, "WE", "Hamburg",
                      5, "WE", "Brunswick",
                      5, "WE", "Oldenburg",
                      5, "WE", "Saxe-Weimar-Eisenach",
                      5, "WE", "Nassau"), nrow = 27, ncol = 3, byrow = TRUE)


for (i in unique(data_raw[is.na(s_mil_region)]$country_name)) {
  tmp_i <- data_tmp[which(data_tmp[, 3] == i), 1]
  data_raw[country_name == i, s_mil_region := as.numeric(tmp_i)]
}


C        <- 1/200
y        <- data_raw$v2x_polyarchy
X        <- as.matrix(cbind(data_raw$s_far_Maddison_gdp_1990_estimate_lag))
Z        <- cbind(1, C*(data_raw$year - 1789), (C*(data_raw$year - 1789))^2, (C*(data_raw$year - 1789))^3, data_raw$s_far_Maddison_pop_estimate)
index    <- data_raw$year
group    <- data_raw$country_name
region   <- data_raw$s_mil_region
country  <- data_raw$country_name
id       <- data_raw$country_id


L               <- length(unique(id))
country_hat_seq <- rep("", L)
tau_hat_seq     <- rep(NA, L)
tau_min_hat_seq <- rep(NA, L)
tau_max_hat_seq <- rep(NA, L)
tau_cc_seq      <- list()
beta_95_seq     <- rep(NA, L)
beta_cc_seq     <- list()
beta_hat_seq    <- rep(NA, L)
beta_ln_seq     <- list()

region_seq      <- rep(NA, L)


country_tmp <- NULL

l <- 1

for (i in sort(unique(id))) {

  index_tmp <- id == i
  y_i       <- y[index_tmp]
  X_i       <- as.matrix(X[index_tmp, ])
  Z_i       <- as.matrix(Z[index_tmp, ])
  group_i   <- group[index_tmp]
  index_i   <- index[index_tmp]
  region_i  <- first(region[index_tmp])

  country_i <- unique(country[index_tmp])


  # Removing with less than 100 observations
  # tmp <- data.table(group = group_r, index = index_r)[, list(nr = length(index)), by = "group"]
  # tmp <- tmp[tmp$nr < 100]$group
  # tmp <- which(!(group_r %in% tmp))
  #
  # y_r       <- y_r[tmp]
  # X_r       <- as.matrix(X_r[tmp, ])
  # Z_r       <- as.matrix(Z_r[tmp, ])
  # group_r   <- group_r[tmp]
  # index_r   <- index_r[tmp]

  # y_r       <- data.table(y_r, group_r)[, c(NA, diff(y_r)), by = 'group_r']$V1

  data_i      <- cdcp.panel.clean.data(y = y_i, X = X_i, Z = Z_i, index = index_i, group = group_i, index_dummy = FALSE, group_dummy = FALSE, lag = TRUE)
  n_i         <- data_i$n


  # In the current model there are 1 parameter for population, 3 parameters for time, 1 intercept,
  # and then 1 for the GDP (on each side). In total 6 parameters on each side. If we are going to
  # have at least 4 data points for each parameter, we need more than 2*4*6 = 48 data points in order
  # too look for a change point

  if (n_i - 2*6*4 > 50) {

    cat(country_i, ": \n", sep = "")
    cat(" n =", n_i, "\n")


    # index_val_r <- range(data_r$index) + c(10, -10)
    # index_val_r <- index_val[1]:index_val[2]

    # index_val_r <- cdcp.panel.clean.data.find.index.val(data = data_r, threshold = 4)$index_val

    # for (j in 1:4) {
    #   try(index_val_i <- cdcp.panel.clean.data.find.index.val(data = data_i, threshold = j)$index_val)
    # }


    # Note:
    # I think this is correct, since we spil y_{t_i} for t_i <= \tau,
    # then there should be the same number to the left and right for the
    # edge cases

    index_val_i <- data_i$index[(6*4):(n_i - 6*4)]


    model_i        <- cdcp.panel.estimate(data = data_i, index_val = index_val_i)


    country_hat_seq[l] <- country_i
    tau_hat_seq[l]     <- model_i$tau
    tau_min_hat_seq[l] <- min(data_i$index)
    tau_max_hat_seq[l] <- max(data_i$index)

    beta_hat_seq[l]    <- model_i$betaL - model_i$betaR

    region_seq[l] <- region_i

    #
    # print(max(model_i$aic))

    if (FALSE) {
      bbridge_i <- cdcp.panel.bridge(data = data_i, index_val = index_val_i)

      model_i$bbridge <- bbridge_i

      cdcp.panel.bridge.plot(model_i)
      # title(main = region_names[r])
    }

    if (TRUE) {

      index_tmp <- index_val_i[model_i$Dn < 25]
      index_tmp <- c(min(index_tmp) + (-2):(-1), index_tmp, max(index_tmp) + 1:2)
      index_tmp <- range(index_tmp)

      index_tmp[1] <- max(min(index_val_i), min(index_tmp))
      index_tmp[2] <- min(max(index_val_i), max(index_tmp))

      index_tmp <- index_tmp[1]:index_tmp[2]

      # data, model, B, index_val, tmp_val, u = NULL, bboot = FALSE

      # print(system.time(sim <- cdcp.panel.function(data_i, model_i, B = 100, index_val = index_val_i, tmp_val = index_tmp, u = seq(0, 1, length.out = 65))))

      print(system.time(sim <- cdcp.panel.function(data_i, model_i, B = 200, index_val = index_val_i, tmp_val = index_tmp, u = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))))


      sim_95 <- ((!is.na(sim$cc_set) + 0)*sim$tau_val)[, 6]
      # sim_95 <- sim_95[sim_95 > 0]

      tau_cc_seq[[l]] <- sim_95

      par(mar = c(2, 2, 1.6, 0.4), mfrow = c(1, 2))
      matplot(x = sim$tau_val, sim$cc_set, pch = 21, bg = 'grey70', col = 'black', axes = FALSE, xlab = '', ylab = '')
      title(main = country_i)
      title(ylab = "Confidence Sets", line = 2.5)
      axis(1)
      axis(2)
      box(bty = 'l')
      abline(h = 0.95, lty = 2, lwd = 1.4)
      # abline(v = model_r$tau, lty = 1.4)

      # Not for r %in% c(8, 7, 6, 5, 4, 1)
      #
      # if (!(r %in% c(8, 7, 6, 5, 4, 1))) {
      #   axis(1, at = model_r$tau, model_r$tau)
      # }
      # text(model_r$tau - 2, 0, model_r$tau, srt = 0, pos = 2)
    }

    # Degree of Change
    if (TRUE) {


      delta_val_i   <- sort(c(0, seq(-0.05, 0.050, length.out = 500)))
      country_doc_i <- cdcp.panel.beta.doc.function(data = data_i, k = 1, delta_val = delta_val_i, index_val = index_val_i)

      beta_ln_seq[[l]] <- country_doc_i$ln_profile

      beta_95_seq[l] <- country_doc_i$cc[country_doc_i$delta_val == 0]

      beta_cc_seq_i <- country_doc_i$delta_val
      beta_cc_seq_i[country_doc_i$cc > 0.95 + 1e-6] <- NA
      beta_cc_seq[[l]] <- beta_cc_seq_i

      with(country_doc_i, plot(delta_val, cc, type = 'l', xlab = '', ylab = '', lwd = 1.6))
      # title(main = region_names[r])
      title(ylab = "Confidence Curve", line = 2.5)
      title(xlab = expression(beta[L] - beta[R]), line = 2.5)
      abline(h = 0.95, lty = 2, col = 'grey50', lwd = 1.6)
      abline(v = 0.00, lty = 2, col = 'grey50', lwd = 1.6)

    }
    l <- l + 1
  } else {
    country_tmp <- c(country_tmp, country_i)
  }
}



tmp_index       <- which(!is.na(tau_hat_seq))
tau_hat_seq     <- tau_hat_seq[tmp_index]
country_hat_seq <- country_hat_seq[tmp_index]

tau_min_hat_seq <- tau_min_hat_seq[tmp_index]
tau_max_hat_seq <- tau_max_hat_seq[tmp_index]

tmp_index_1700  <- which(tau_min_hat_seq < 1800)
tmp_index_1800  <- which(1800 <= tau_min_hat_seq & tau_min_hat_seq < 1900)
tmp_index_1900  <- which(1900 <= tau_min_hat_seq & tau_min_hat_seq < Inf)

region_col <- rainbow(8)
tmp_order <- order(tau_hat_seq)
tmp_n <- length(tmp_order)

par(mar = c(4, 9, 2.2, 4.5), mfrow = c(1, 1))
matplot(matrix(rep(c(1700, 2020), times = l - 1), ncol = l - 1), matrix(rep(1:(l - 1), each = 2), ncol = l - 1), type = 'l', col = 'grey86', lty = 1, xlim = c(1800, 2000), xlab = "", ylab = "",  ylim = c(3, tmp_n - 2), axes = FALSE)



# Plot of CC set for tau
i <- 1
for (j in tmp_order) {

  plot_x <- tau_cc_seq[[j]]
  plot_x <- plot_x[first(which(plot_x > 0)):last(which(plot_x > 0))]

  plot_y <- ifelse(plot_x == 0, NA, i)
  plot_x <- range(plot_x[plot_x > 0])
  plot_x <- plot_x[1]:plot_x[2]

  lines(plot_x, plot_y, type = 'l', lwd = 2.5)
  points(tau_hat_seq[j], i, type = 'p', pch = 21, cex = 0.9, bg = ifelse(beta_95_seq[j] > 0.95, "black", "white"))
  axis(side = 2, at = i, labels = country_hat_seq[j], las = 2, cex.axis = 0.80, col.axis = region_col[region_seq[j]])

  i <- i + 1
}

box(bty = 'l')
axis(side = 1)
# axis(side = 2, at = 1:tmp_n, labels = country_hat_seq[tmp_order], las = 2, cex.axis = 0.9)

title(xlab = "year", line = 2.5, cex.lab = 1.1)
title(main = "0.95 confidence sets", line = 1.0, cex.lab = 1.3)

par(xpd=TRUE)
legend("topright", inset = c(-0.20, 0), col = region_col, legend = paste("Region", 1:8), pch = 19, cex = 0.9, x.intersp = 0.3, y.intersp = 0.2, bty = 'n')
par(xpd=FALSE)





# Degree of Change at 0.95 level
par(mar = c(4, 9, 2.2, 4.5), mfrow = c(1, 1))
matplot(matrix(rep(c(-0.047, 0.047), times = l - 1), ncol = l - 1), matrix(rep(1:(l - 1), each = 2), ncol = l - 1), type = 'l', col = 'grey86', lty = 1, xlim = c(-0.047, 0.047), xlab = "", ylab = "", axes = FALSE, ylim = c(3, tmp_n - 2))

tmp_order <- order(beta_hat_seq[!is.na(beta_hat_seq)])

i <- 1
for (j in tmp_order) {

  plot_x <- delta_val_i
  plot_y <- beta_cc_seq[[j]]
  plot_y[!is.na(plot_y)] <- i
  # plot_x <- plot_x[first(which(plot_x > 0)):last(which(plot_x > 0))]

  # plot_y <- ifelse(plot_x == 0, NA, i)
  # plot_x <- range(plot_x[plot_x > 0])
  # plot_x <- plot_x[1]:plot_x[2]

  lines(plot_x, plot_y, type = 'l', lwd = 2.5)
  points(beta_hat_seq[j], i, type = 'p', pch = 21, cex = 0.9, bg = ifelse(beta_95_seq[j] > 0.95, "black", "white"))
  axis(side = 2, at = i, labels = country_hat_seq[j], las = 2, cex.axis = 0.80, col.axis = region_col[region_seq[j]])

  i <- i + 1
}

box(bty = 'l')
axis(side = 1)
# axis(side = 2, at = 1:tmp_n, labels = country_hat_seq[tmp_order], las = 2, cex.axis = 0.9, col.axis = region_col[region_seq[tmp_order]])
abline(v = 0, col = "black", lty = 2)
title(xlab = expression(beta[L] - beta[R]), line = 2.5, cex.lab = 1.1)
title(main = "0.95 confidence interval", line = 1.0, cex.lab = 1.3)

par(xpd=TRUE)
legend("topright", inset = c(-0.17, 0), col = region_col, legend = paste("Region", 1:8), pch = 19, cex = 0.9, x.intersp = 0.3, y.intersp = 0.2, bty = 'n')
par(xpd=FALSE)







# Combined Confidence Curve
par(mar = c(4, 4, 2, 8.9), xpd=FALSE, bty = 'l')
plot(delta_val_i, delta_val_i*0, type = 'l', ylim = c(0, 1), lwd = 1, col = "white", xlab = "", ylab = "", xlim = c(-0.045, 0.025))

beta_ln_all <- delta_val_i*0
for (j in tmp_order) {
  beta_ln_all <- beta_ln_all + beta_ln_seq[[j]]

  lines(delta_val_i, pchisq(2*(max(beta_ln_seq[[j]]) - beta_ln_seq[[j]]), df = 1), lwd = 1.5, lty = 1, type = 'l', ylim = c(0, 1), col = region_col[region_seq[j]])
}

lines(delta_val_i, pchisq(2*(max(beta_ln_all) - beta_ln_all), df = 1), lwd = 2.2, type = 'l', ylim = c(0, 1), col = "black")
abline(h = 0, lwd = 1.8)
abline(h = 0.95, lwd = 1.8, lty = 2)
abline(v = 0, lwd = 1.8, lty = 2)
title(xlab = expression(beta[L] - beta[R]), line = 2.5, cex.lab = 1.3)
title(ylab = "confidence curve", line = 2.5, cex.lab = 1.3)

par(xpd=TRUE)
legend("topright", inset = c(-0.318, 0), col = c(region_col, "black"), legend = c(paste("Region", 1:8), "Combined"), lwd = 1.7, cex = 1.00, bty = 'n', x.intersp = 0.4, y.intersp = 0.7, seg.len = 0.9)
par(xpd=FALSE)






if (FALSE) {
  plot_country <- country_hat_seq[tmp_index_1700]
  plot_tau_hat <- tau_hat_seq[tmp_index_1700]
  plot_order <- order(plot_tau_hat)
  plot_n <- length(plot_tau_hat)


  par(mar = c(2, 10, 1, 1))
  plot(plot_tau_hat[plot_order], 1:plot_n, axes = FALSE, xlab = "", ylab = "", pch = 21, bg = "grey70", cex = 2, main = "Countries with Data Dating Back to the 18th Century")
  box(bty = 'l')
  axis(side = 1)
  axis(side = 2, at = 1:plot_n, labels = plot_country[plot_order], las = 2, cex.axis = 0.9)

  plot_country <- country_hat_seq[tmp_index_1800]
  plot_tau_hat <- tau_hat_seq[tmp_index_1800]
  plot_order <- order(plot_tau_hat)
  plot_n <- length(plot_tau_hat)

  par(mar = c(2, 10, 1, 1))
  plot(plot_tau_hat[plot_order], 1:plot_n, axes = FALSE, xlab = "", ylab = "", pch = 21, bg = "grey70", cex = 2, main = "Countries with Data Dating Back to the 19th Century")
  box(bty = 'l')
  axis(side = 1)
  axis(side = 2, at = 1:plot_n, labels = plot_country[plot_order], las = 2, cex.axis = 0.9)


  plot_country <- country_hat_seq[tmp_index_1900]
  plot_tau_hat <- tau_hat_seq[tmp_index_1900]
  plot_order <- order(plot_tau_hat)
  plot_n <- length(plot_tau_hat)

  par(mar = c(2, 10, 1, 1))
  plot(plot_tau_hat[plot_order], 1:plot_n, axes = FALSE, xlab = "", ylab = "", pch = 21, bg = "grey70", cex = 2, main = "Countries with Data Dating Back to the 20th Century")
  box(bty = 'l')
  axis(side = 1)
  axis(side = 2, at = 1:plot_n, labels = plot_country[plot_order], las = 2, cex.axis = 0.9)
}






# par(mar = c(2, 10, 1, 1))
# plot(tau_hat_seq[tmp_order], 1:tmp_n, axes = FALSE, xlab = "", ylab = "")
# box(bty = 'l')
# axis(side = 1)
# axis(side = 2, at = 1:tmp_n, labels = country_hat_seq[tmp_order], las = 2, cex.axis = 0.6)

# par(mfrow = c(1, 1))
# plot(tau_hat_seq[tmp_index_1700], tau_hat_seq[tmp_index_1700]*0)
# text(tau_hat_seq[tmp_index_1700], y = 0.4, labels = country_hat_seq[tmp_index_1700], srt = 90)
#
# par(mfrow = c(1, 1))
# plot(tau_hat_seq[tmp_index_1800], tau_hat_seq[tmp_index_1800]*0)
# text(tau_hat_seq[tmp_index_1800], y = 0.4, labels = country_hat_seq[tmp_index_1800], srt = 90)
#
# par(mfrow = c(1, 1))
# plot(tau_hat_seq[tmp_index_1900], tau_hat_seq[tmp_index_1900]*0)
# text(tau_hat_seq[tmp_index_1900], y = 0.4, labels = country_hat_seq[tmp_index_1900], srt = 90)











































