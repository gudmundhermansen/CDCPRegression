# CDCPRegression
Confidence Distribution for Change Points for Panel Regression

# Introduction and Summary 

This package contains a simplified version of the main methods used in Hermansen (2021) for analysis of change point in panel regression models. The package is designed for models of the type $$Y_{i, j} = x_{i, j}^{\rm t} \beta_{\rm L} + x_{i, j}^{\rm t} (\beta_{\rm R} - \beta_{\rm L}) + z_{i, j}^{\rm t} \gamma + d_{i, j}^{\rm t} \delta + \rho y_{i - 1, j} + \epsilon_{i, j},$$ where $x_{i, j}$ are the explanatory variables we believe are affected by a change point. The explanatory variables in $z_{i, j}$ and $\delta_{i, j}$ are assumed to not be affected by the change point. In the model above $\delta_{i, j}$ represent categorical (or dummy) explanatory variables, and the model is written in this way to mimic the implementation below, since categorical requires some special treatment to ensure efficient estimation. Finally, we assume that $\epsilon_{i, j}$ are independent and typically $\epsilon_{i, j} \sim {\rm N}(0, \sigma^2)$ (some robust options are available). Note that the package has some special functions to handle categorical for each index, i.e. for each index $j$ in the model above (e.g. year), and also group categorical variables, for a given partition of $i$ into groups (collection of countries). 


# Installation 
To install the `CDCPRegression` package, install `devtools`, `data.table`, `foreach` and `doMC` below and run the following set of commands in R: 
```
library(devtools)
library(data.table)
library(foreach)
library(doMC)

install_git("https://github.com/gudmundhermansen/CDCPRegression.git")

library(CDCPRegression)
```

# Note

The package is still under development, it must therefore be used with a bit of caution and it is not guaranteed to be completely error-free. The bootstrap routine used to compute the confidence sets for the location of a change point, and the non-approximation confidence curve for the so-called degree of change, are both quite computer intensive. To reduce the computational complexity, it is advised to do an initial search for reasonable values of `index_val_sub` and `delta_val` with a low number of bootstrap samples (`boot = 25`) before increasing this to a more reasonable number of bootstrap samples (100+). 
 

# Examples: 

In order to have an efficient implementation, and also to treat categorical, or dummy, variables correctly with a change point, a specific structure is required to compute the models in the `CDCPRegression` package. To properly structure the data for a given model, data must be pre-processed using the `cdcp.regression.data(...)` function. The examples below are intended to make this setup more or less self-explanatory.

**Example 1 -- 2:** Two simple illustrations with (simulated) data from the model $$y_{i, j} = x_{i, j}^{\rm t} \beta_{\rm L} + x_{i, j}^{\rm t} (\beta_{\rm R} - \beta_{\rm L}) + z_{i, j} + \epsilon_{i, j},$$ first for a single and then for a pair of individuals (e.g. countries), with a change point in the intercept. 
```
# Example 1: One individual with Change Point in the Intercept

# Data

set.seed(1) # set seed to make sure delta_val spans the relevant range of values

n <- 100
index <- 1:n
group <- rep(1, n)
mu <- 1.5*(index <= 50)

y <- mu + rnorm(n = n)
X <- matrix(1, nrow = n)

# Plot of data
plot(index, y)
lines(index, mu, type = 's')

data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
index_val <- cdcp.regression.data.find.index.val(data)$index_val

# Plot Monitoring Bridge 
bridge <- cdcp.regression.bridge(data, index_val)
cdcp.regression.bridge.plot(bridge)

# Estimation 
model <- cdcp.regression.estimate(data, index_val)

# Plot Confidence Sets 
cd <- cdcp.regression(data, model, index_val)
cdcp.regression.plot(cd)

# Plot Confidence Curve for Degree of Change in the Intercept 
delta_val <- seq(-0.5, 2.5, length.out = 100) # set of plausible values for the the change in 
                                              # intercept at the change point 
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2) # the 0.95 confidence set
abline(v = 0.00, lty = 2) # zero 
abline(v = 1.50, lty = 1) # true change in intercept



# Example 2: Two individuals with Change Point in the Intercept

set.seed(1)

m <- 2
n <- 100
index <- rep(1:n, times = m)
group <- rep(1:m, each = n)
mu <- 1.5*(index <= 50)

y <- mu + rnorm(n = m*n)
X <- matrix(1, nrow = n*m)

# Plot of data
plot(index, y, col = rep(1:m, each = n))
lines(index[1:n], mu[1:n], type = 's')

data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
index_val <- cdcp.regression.data.find.index.val(data)$index_val

# Plot Monitoring Bridge 
bridge <- cdcp.regression.bridge(data, index_val)
cdcp.regression.bridge.plot(bridge)

# Estimation 
model <- cdcp.regression.estimate(data, index_val)

# Plot Confidence Sets 
cd <- cdcp.regression(data, model, index_val)
cdcp.regression.plot(cd)

# Plot Confidence Curve for Degree of Change in the Intercept 
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2)
abline(v = 0.00, lty = 2) # zero 
abline(v = 1.50, lty = 1) # true change in intercept 
```

**Example 3:** Simulated data of one individual, however, this time the change point is in the slope and not the intercept. 

```
# Example 3: One individual with Change Point in the Slope

set.seed(1) # set seed to make sure delta_val spans the relevant range of values

# Data
n <- 100
index <- 1:n
group <- rep(1, n)
mu <- 1.0 + 3.0*(index <= 50)*(index/n)

y <- mu + rnorm(n = n)
X <- cbind(1, index/n)

# Plot of data
plot(index, y)
lines(index, mu, type = 's')

data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
index_val <- cdcp.regression.data.find.index.val(data)$index_val

# Plot Monitoring Bridge 
bridge <- cdcp.regression.bridge(data, index_val)
cdcp.regression.bridge.plot(bridge)

# Estimation 
model <- cdcp.regression.estimate(data, index_val)

# Plot Confidence Sets 
cd <- cdcp.regression(data, model, index_val)
cdcp.regression.plot(cd)

# Plot Confidence Curve for Degree of Change in the Intercept (where there is no change)
delta_val <- seq(-3, 3.5, length.out = 100) # finding a relevant range of values is usually an iterative process
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2) # the 0.95 confidence set
abline(v = 0.00, lty = 2) # zero 

# Plot Confidence Curve for Degree of Change in the Slope (where there is a change point) 
delta_val <- seq(-1, 7, length.out = 100)
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 2, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2) # the 0.95 confidence set
abline(v = 0.00, lty = 2) # zero 
abline(v = 3.00, lty = 1) # true change 
```

**Example 4:** Simulated data from a simple linear regression model, for one individual, with a lagged response and a change point in the intercept. However, here we first fit two misspecified models, both missing the slope and one also missing the lagged response in the model. In the final model, the linear increasing variable is included as a protected variable, since this is "believed" to not be affected by the change point.  

```
# Example 4: One individual, Slope, Lagged Response and Change Point in the Intercept

set.seed(1)

# Data
n <- 100
index <- 1:n
group <- rep(1, n)
mu <- 3.0*(index <= 50) + 1.0*(index/n)

y <- rep(NA, n)
y[1] <- mu[1] + rnorm(1)

for (i in 2:n) {
  y[i] <- mu[i] + 0.3*y[i - 1] + 1.0*rnorm(1)
}

# Plot of data
plot(index, y)
lines(index, mu, type = 's')

# Model 1:

X <- matrix(1, nrow = n, ncol = 1)

# Note that model is missing the slope and lagged response
data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
index_val <- cdcp.regression.data.find.index.val(data)$index_val

# Plot Monitoring Bridge 
bridge <- cdcp.regression.bridge(data, index_val)
cdcp.regression.bridge.plot(bridge)

# Estimation 
model <- cdcp.regression.estimate(data, index_val)

# Plot Confidence Sets 
cd <- cdcp.regression(data, model, index_val)
cdcp.regression.plot(cd)

# Plot Confidence Curve for Degree of Change in the Intercept 
delta_val <- seq(-0.5, 7.0, length.out = 100) # set of plausible values for the the change in intercept at the change point 

cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2)
abline(v = 0.00, lty = 2)
abline(v = 3.00, lty = 2)


# Model 2: 

# same data, but now with the lagged response included
data <- cdcp.regression.data(y = y, X = X, index = index, group = group, lag = TRUE)
index_val <- cdcp.regression.data.find.index.val(data)$index_val

# Plot Monitoring Bridge 
bridge <- cdcp.regression.bridge(data, index_val)
cdcp.regression.bridge.plot(bridge)

# Estimation 
model <- cdcp.regression.estimate(data, index_val)

# Plot Confidence Sets 
cd <- cdcp.regression(data, model, index_val)
cdcp.regression.plot(cd)

# Plot Confidence Curve for Degree of Change in the Intercept 
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2)
abline(v = 0.00, lty = 2)
abline(v = 3.00, lty = 2)


# Model 3: 

# Protected variable, i.e. a variable that we think is not affected by the change point 
Z <- matrix(index/n)

# the correctly specified model
data <- cdcp.regression.data(y = y, X = X, Z = Z, index = index, group = group, lag = TRUE)
index_val <- cdcp.regression.data.find.index.val(data)$index_val

# Plot Monitoring Bridge 
bridge <- cdcp.regression.bridge(data, index_val)
cdcp.regression.bridge.plot(bridge)

# Estimation 
model <- cdcp.regression.estimate(data, index_val)

# Plot Confidence Sets 
cd <- cdcp.regression(data, model, index_val)
cdcp.regression.plot(cd)

# Plot Confidence Curve for Degree of Change in the Intercept 
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2)
abline(v = 0.00, lty = 2)
abline(v = 3.00, lty = 2)
```


# Refrences  

Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504

Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.





