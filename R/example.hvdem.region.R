

# rm(list = ls())

library(data.table)
library(readstata13)
library(foreign)
library('doMC')
registerDoMC(4)

# Read Raw Data
# data_raw <- read.dta13(file = "~/Dropbox/shared/Change points/Data/HVDEM_Polyarchy_with_region_covariates.dta")
# data_raw <- as.data.table(data_raw)
# data_raw[, s_far_Maddison_gdp_1990_estimate_lag := c(NA, s_far_Maddison_gdppc_1990_estima[-.N]), by = 'country_name']

# data_raw <- read.dta13(file = "~/Dropbox/shared/Change points/Data/HVDEM_Polyarchy_with_region_covariates_updated.dta")
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



# data_raw <- read.table(file   = '~/Dropbox/shared/Change points/Data/hvdem_region_sub.csv',
#                         header = TRUE,
#                         sep    = ',',
#                         stringsAsFactors = FALSE)
#
# y      <- data_raw$v2x_polyarchy
#
# X      <- cbind(1, data_raw$s_far_Maddison_gdppc_1990_estimate_lag)
# Z      <- cbind(data_raw$s_far_Maddison_pop_estimate)
# index  <- data_raw$year
# group  <- data_raw$country_name
# region <- data_raw$s_mil_region


region_names <- c('Eastern Europe and Soviet Space',
                  'Latin America',
                  'Middle East and North Africa',
                  'Sub-Saharan Africa',
                  'Western Europe and North America',
                  'East Asia',
                  'South-East Asia',
                  'South Asia')


# par(mfrow = c(2, 2), mar = c(4, 4, 1.6, 0.4), bty = 'l')
par(mfrow = c(1, 3), mar = c(4, 4, 1.6, 0.4), bty = 'l', cex.lab = 1.4)

for (r in sort(unique(region))) {
  # For region 6 only use 1 data point??
  #
  #

  index_tmp <- region == r
  y_r       <- y[index_tmp]
  X_r       <- as.matrix(X[index_tmp, ])
  Z_r       <- as.matrix(Z[index_tmp, ])
  group_r   <- group[index_tmp]
  index_r   <- index[index_tmp]

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

  data_r      <- cdcp.panel.clean.data(y = y_r, X = X_r, Z = Z_r, index = index_r, group = group_r, index_dummy = FALSE, group_dummy = TRUE, lag = TRUE)
  # index_val_r <- range(data_r$index) + c(10, -10)
  # index_val_r <- index_val[1]:index_val[2]

  # index_val_r <- cdcp.panel.clean.data.find.index.val(data = data_r, threshold = 4)$index_val

  for (j in 1:4) {
     try(index_val_r <- cdcp.panel.clean.data.find.index.val(data = data_r, threshold = j)$index_val)
  }

  model_r <- cdcp.panel.estimate(data = data_r, index_val = index_val_r)
  print(max(model_r$aic))

  if (TRUE) {
    bbridge_r <- cdcp.panel.bridge(data = data_r, index_val = index_val_r)

    model_r$bbridge <- bbridge_r

    cdcp.panel.bridge.plot(model_r)
    # title(main = region_names[r])
  }


  if (TRUE) {
    # cdcp.panel.summary.plot.deviance(data_r, model_r)

    index_tmp <- index_val_r[model_r$Dn < 20]
    index_tmp <- c(min(index_tmp) + (-2):(-1), index_tmp, max(index_tmp) + 1:2)
    index_tmp <- range(index_tmp)

    index_tmp[1] <- max(min(index_val_r), min(index_tmp))
    index_tmp[2] <- min(max(index_val_r), max(index_tmp))

    index_tmp <- index_tmp[1]:index_tmp[2]

    print(system.time(sim <- cdcp.panel.function(data_r, model_r, B = 100, index_val = index_val_r, tmp_val = index_tmp, u = seq(0, 1, length.out = 65))))

    #par(mar = c(2, 2, 1.6, 0.4), mfrow = c(1, 1))
    matplot(x = sim$tau_val, sim$cc_set, pch = 21, bg = 'grey70', col = 'black', axes = FALSE, xlab = '', ylab = '')
    title(main = region_names[r])
    title(ylab = "Confidence Sets", line = 2.5)
    axis(1)
    axis(2)
    box(bty = 'l')
    abline(h = 0.95, lty = 2, lwd = 1.4)
    # abline(v = model_r$tau, lty = 1.4)

    # Not for r %in% c(8, 7, 6, 5, 4, 1)
    #
    if (!(r %in% c(8, 7, 6, 5, 4, 1))) {
      axis(1, at = model_r$tau, model_r$tau)
    }
    # text(model_r$tau - 2, 0, model_r$tau, srt = 0, pos = 2)
  }

  # Degree of Change
  if (TRUE) {
    # r == 1 => 1922 or 1988 [-0.07, 0.11]
    # r == 2 => No clear change point [-0.05, 0.11]
    # r == 3 => No clear change point [-0.05, 0.05]
    # r == 4 => 2004 [-0.1, 0.1]
    # r == 5 => 1973 [-0.1, 0.1]
    # r == 6 => 1951 (not a lot of data, limit to 1) [-0.1, 0.1]
    # r == 7 => No clear change point (not a lot of data, limit set to 2) [-0.05, 0.2] with limit 1 => [-0.40, 0.2]
    # r == 8 => 2003 (not a lot of data, limit set to 2) [-0.07, 0.11]

    delta_range <- matrix(c(-0.015, 0.05,
                            -0.005, 0.06,
                            -0.060, 0.00,
                            -0.000, 0.04,
                             0.015, 0.06,
                            -0.010, 0.04,
                            -0.040, 0.00,
                            -0.005, 0.06), byrow = TRUE, nrow = 8, ncol = 2)


    # delta_val_r <- seq(min(delta_range[, 1]), max(delta_range[, 2]), length.out = 100)
    delta_val_r <- seq(-0.014, 0.003, length.out = 100)
    region_doc <- cdcp.panel.beta.doc.function(data = data_r, k = 1, delta_val = delta_val_r, index_val = index_val_r)

    with(region_doc, plot(delta_val, cc, type = 'l', xlab = '', ylab = '', lwd = 1.6))
    # title(main = region_names[r])
    title(ylab = "Confidence Curve", line = 2.5)
    title(xlab = expression(beta[L] - beta[R]), line = 2.5)
    abline(h = 0.95, lty = 2, col = 'grey50', lwd = 1.6)
    abline(v = 0.00, lty = 2, col = 'grey50', lwd = 1.6)

  }
}











