# CDCPRegression
Confidence Distribution for Change Points for Panel Regression

# Introduction and Summary 

This package contain the main methods used in ??? for change point analysis in panel regression models. 


# Note

This package is still under development, and must therefore be used with a bit of caution and is not guaranteed to be completely error-free. The bootstrap routine used to compute the confidence sets for the change point, and the non-approximation confidence curve for the so-called degree of change is quite computer intensive. To reduce the computational complexity, it is advised to do an initial search for reasonable values of `index_val_sub` and `delta_val` with a low number of bootstrap samples before increasing this to 100-1000 samples. 
 


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

# Examples: 
```
# Example 1: Simulated Data of One individual with Change Point in the Intercept

set.seed(1)

# Data
n <- 100
index <- 1:n
group <- rep(1, n)
mu <- 3.0*(index <= 50) + 1.0*(index/n)

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
delta_val <- seq(-0.5, 4.5, length.out = 100) # set of plausible values for the the change in intercept at the change point 
cc <- cdcp.regression.beta.doc(data, model, index_val, k = 1, delta_val = delta_val)
cdcp.regression.beta.doc.plot(cc)
abline(h = 0.95, lty = 2) # the 0.95 confidence set
abline(v = 0.00, lty = 2) # zero 
abline(v = 3.00, lty = 1) # true change in intercept

```
