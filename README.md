# plmR: Partial Linear Model for High-Dimensional Data Analysis

<!--[![R-CMD-check](https://github.com/SangkyuStat/plmR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SangkyuStat/plmR/actions/workflows/R-CMD-check.yaml)```-->

This R package, plmR, implements a high-dimensional partial linear model as well as a high-dimensional partial linear quantile regression model using LASSO and Quantile LASSO. 
 
Our package utilize **trend filtering** and **smoothing spline** for the nonparametric term, and is also adapted to the partial linear quantile regression with these two methods. Our model can be simply reduced to the models without the penalties, by letting penalty to be 0. For both methods, `lambda1` is a penalty related to LASSO and `lambda2` is a penalty related to the nonparametric term.

All of our functions are built on Block-wise Coordinate Descent algorithm (Backfitting algorithm). However, `plqrtf` and `plqrss` are designed with one more nested loop for approximation, therefore, it has to be used with care because of the time consumption. 

### Installation

The current version can be installed from source using the package `devtools`
```R
devtools::install_github("SangkyuStat/plmR")
```

### Usage Examples

#### (1) `pltf` and `plqrtf` functions
`pltf` function provides partial linear trend filtering. Its argument incorporates `y`, `x` and `z`, which stand for response variable, covariates for linear term and a covariate for nonparametric term. Note that `z` should be a variable, and `x` can be a matrix. `k` is the order for the bounded total variation of the derivative of function. `lambda1` is the tuning parameter for the linear part, and `lambda2` is the tuning parameter for the nonparametric part.

`plqrtf` function provides partial linear quantile trend filtering. Its argument incorporates `y`, `x` and `z`, which stand for response variable, covariates for linear term and a covariate for nonparametric term. Note that `z` should be a variable, and `x` can be a matrix. `k` is the order for the bounded total variation of the derivative of function. `lambda1` is the tuning parameter for the linear part, and `lambda2` is the tuning parameter for the nonparametric part. Additionally, there is one more argument `tau`, which stands for the specific conditional quantile level of interest, if nothing is given, then median is set as a default.

```R
plqr.fit(y = response, x = design matrix, z = covariate, k = order,
         lambda1 = some_value1, lambda2 = some_value2)

plqrtf.fit(y = response, x = design matrix, z = covariate, k = order, tau = 0.5)
           lambda1 = some_value1, lambda2 = some_value2)
```

<!--the user has to indicate whether the variable is time-varying or not. If there is no time-varying variable, then user can perform the function as below:
```R
vc.pb(formula = response ~ variable + 
                any modifier, 
                id,
                data = input_data, 
                group = disparity_group, 
                modifier = "any modifier")
```-->
#### (2) `plss` and `plqrss` functions
`plss` function provides partial linear smoothing spline. Its argument incorporates `y`, `x` and `z`, which stand for response variable, covariates for linear term and a covariate for nonparametric term. Note that `z` should be a variable, and `x` can be a matrix. By now, order is fixed as a cubic smoothing spline, but this will be updated soon. `lambda1` is the tuning parameter for the linear part, and `lambda2` is the tuning parameter for the nonparametric part.

`plqrss` function provides partial linear quantile smoothing spline. Its argument incorporates `y`, `x` and `z`, which stand for response variable, covariates for linear term and a covariate for nonparametric term. Note that `z` should be a variable, and `x` can be a matrix. By now, order is fixed as a cubic smoothing spline, but this will be updated soon. `lambda1` is the tuning parameter for the linear part, and `lambda2` is the tuning parameter for the nonparametric part. Additionally, there is one more argument `tau`, which stands for the specific conditional quantile level of interest, if nothing is given, then median is set as a default. 

```R
plss.fit(y = response, x = design matrix, z = covariate,
         lambda1 = some_value1, lambda2 = some_value2)

plqrss.fit(y = response, x = design matrix, z = covariate, tau = 0.5,
           lambda1 = some_value1, lambda2 = some_value2)
```

### Developing

- The Rcpp version is being tested now for stability.
- The cross-validation function for choosing the tuning parameters will be updated soon.
- We are trying to develop other methods as well.

### Version

- 0.0.1: Initial Version (Current Version)

<!--
### References

Peters, C. C. (1941) A method of matching groups for experiment with no loss of population. Journal of Educational Research, 34, 606-612.

Belson, W. A. (1956) A Technique for Studying the Effects of a Television 
Broadcast.  JRSSC, 5(3), 195-202.

Lee, S. K., Kim, S., Kim, M.-O., Grantz, K. L., and Hong, H. G. (2024) Decomposition of Longitudinal Disparities: An Application to the Fetal Growth-Singletons Study. *submitted*.
```-->