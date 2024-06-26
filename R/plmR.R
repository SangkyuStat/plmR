#' Partial Linear Quantile Trend Filtering for a High-dimensional Heterogeneous Data
#'
#' The \code{plqr.fit} performs partial linear quantile trend filtering model on the high-dimensional data by combining Quantile LASSO and trend filtering.
#' @param y a vector of a response variable for the model.
#' @param x a matrix of covariates related to the linear part.
#' @param z a vector of covariate related to the linear part.
#' @param tau a real number for a quantile level ranged between 0 and 1.
#' @param lambda1 a real number which is the tuning parameter for LASSO
#' @param lambda2 a real number which is the tuning parameter for trend filering
#' @param k an integer for order of the bounded total variation, default is 2.
#' @param iter.max a maximum iteration numer for the blockwise coordinate descent (BCD).
#' @param w.iter.max a maximum iteration numer for the MM algorithm for quantile regression.
#' @param tol1 a tolerance error when the BCD loop break.
#' @param tol2 a tolerance error when the MM loop break.
#' @param init.g (optional) a vector used for initialization of nonparametric term.
#' @param ... used for controlling the others.
#' @author Sang Kyu Lee
#' @return \code{plqr.fit} returns an object of class \code{"list"}, which is a list containing
#' following components:
#' @return
#' \item{beta}{fitted coefficients for the linear term.}
#' \item{g}{fitted values for nonparametric terms.}
#' \item{model.g}{fitted nonparmetric object.}
#' @importFrom rqPen rq.lasso.fit
#' @import glmgen
#' @export
plqr.fit = function(y, x, z, tau, lambda1, lambda2, k=2,
                    iter.max = 1000, w.iter.max = 1000, 
                    tol1 = 1e-3, tol2 = 1e-4, 
                    init.g = NULL, ...){
  n = length(y)
  p = dim(x)[2]
  for(iter in 1:iter.max){
    if(iter == 1){
      if(is.null(init.g)){
        g.fit = rep(0, n)
      } else {
        g.fit = init.g
      }
    } else {
      beta.old = beta.new
      g.old.total = g.fit
    }     
    y_temp1 = y - g.fit
    lasso.fit = rqPen::rq.lasso.fit(x = x, y = y_temp1, tau = tau, lambda = lambda1,
                                     penalty = "LASSO", alg = "br", intercept = FALSE, scalex = FALSE)

    beta.new = lasso.fit$coefficients 
    y_temp2 = y - x%*%beta.new
    for(w.iter in 1:w.iter.max){
      if(w.iter == 1 & iter == 1){
        g.old = g.fit
      } else if(w.iter == 1 & iter > 1){
        g.old = g.old.total
      } else {
        g.old = g.fit
      }
      weight.tf = weight.func(yw = y_temp2, g = g.old, tau = tau)
      tf.fit = trendfilter(x = z, y = y_temp2, weights = weight.tf,
                                   family = "gaussian", k=k,
                                   lambda = lambda2)
      g.fit = predict(tf.fit, x.new = (z), lambda = lambda2)

      if(sum((g.fit - g.old)^2) < tol2) break
    }
    if(iter!=1){
      check.error = abs(check.sum(y = y, x = x, beta = beta.new, g = g.fit, tau = tau) - check.sum(y = y, x = x, beta = beta.old, g = g.old.total, tau = tau))

      if(sum((x%*%(beta.new - beta.old) + g.fit - g.old.total)^2) < tol1) break
    } else {
      next
    }
  }
  res = list(beta = beta.new, g = g.fit, model.g = tf.fit)
  res
}

weight.func = function(yw, g, eta = 1e-4, tau){
  (tau * I(yw > as.numeric(g)) + (1-tau)*(yw <= as.numeric(g)))/
    sqrt((yw - as.numeric(g))^2 + eta^2)
}

lm.loss = function(y, x, beta, g){
  n = length(y)
  sum((y - x%*%beta - g)^2)/n
}

check.loss = function(x, tau){
  x*(tau - I(x<0))
}

check.sum = function(y, x, beta, g, tau){
  n = length(y)
  sum(check.loss(y - x%*%as.vector(beta) - g, tau))/n
}

#' Partial Linear Quantile Smoothing Spline for a High-dimensional Heterogeneous Data
#'
#' The \code{plqrss.fit} performs partial linear quantile trend filtering model on the high-dimensional data by combining Quantile LASSO and smoothing spline.
#' @param y a vector of a response variable for the model.
#' @param x a matrix of covariates related to the linear part.
#' @param z a vector of covariate related to the linear part.
#' @param tau a real number for a quantile level ranged between 0 and 1.
#' @param lambda1 a real number which is the tuning parameter for LASSO
#' @param lambda2 a real number which is the tuning parameter for trend filering
#' @param iter.max a maximum iteration numer for the blockwise coordinate descent (BCD).
#' @param w.iter.max a maximum iteration numer for the MM algorithm for quantile regression.
#' @param tol1 a tolerance error when the BCD loop break.
#' @param tol2 a tolerance error when the MM loop break.
#' @param init.g (optional) a vector used for initialization of nonparametric term.
#' @param ... used for controlling the others.
#' @author Sang Kyu Lee
#' @return \code{plqrss.fit} returns an object of class \code{"list"}, which is a list containing
#' following components:
#' @return
#' \item{beta}{fitted coefficients for the linear term.}
#' \item{g}{fitted values for nonparametric terms.}
#' \item{model.g}{fitted nonparmetric object.}
#' @importFrom rqPen rq.lasso.fit
#' @import stats
#' @export
plqrss.fit = function(y, x, z, tau, lambda1, lambda2, 
                      iter.max = 1000, w.iter.max = 100, 
                      tol1 = 1e-3, tol2 = 1e-4, 
                      init.g = NULL, ...){
  n = length(y)
  p = dim(x)[2]
  for(iter in 1:iter.max){
    if(iter == 1){
      if(is.null(init.g)){
        g.fit = rep(0, n)
      } else {
        g.fit = init.g
      }
    } else {
      beta.old = beta.new
      g.old.total = g.fit
    }     
    y_temp1 = y - g.fit

    lasso.fit = rqPen::rq.lasso.fit(x = x, y = y_temp1, tau = tau, lambda = lambda1,
                                     penalty = "LASSO", alg = "br", intercept = FALSE, scalex = FALSE)
    beta.new = lasso.fit$coefficients 
    y_temp2 = y - x%*%beta.new
    for(w.iter in 1:w.iter.max){
      if(w.iter == 1 & iter == 1){
        g.old = g.fit
      } else if(w.iter == 1 & iter > 1){
        g.old = g.old.total
      } else {
        g.old = g.fit
      }
      weight.tf = weight.func(yw = y_temp2, g = g.old, tau = tau)
      ss.fit = smooth.spline(x = z, y = y_temp2, w = weight.tf,
                             lambda = lambda2)
      
      g.fit = predict(ss.fit, x = (z), lambda = lambda2)$y
      
      if(sum((g.fit - g.old)^2) < tol2) break
    }
    if(iter!=1){
      check.error = abs(check.sum(y = y, x = x, beta = beta.new, g = g.fit, tau = tau) - check.sum(y = y, x = x, beta = beta.old, g = g.old.total, tau = tau))
      if(sum((x%*%(beta.new - beta.old) + g.fit - g.old.total)^2) < tol1) break
    } else {
      next
    }
  }
  res = list(beta = beta.new, g = g.fit, model.g = ss.fit)
  res
}



#' Partial Linear Trend Filtering for a High-dimensional Heterogeneous Data
#'
#' The \code{pltf.fit} performs partial linear quantile trend filtering model on the high-dimensional data by combining Quantile LASSO and trend filtering.
#' @param y a vector of a response variable for the model.
#' @param x a matrix of covariates related to the linear part.
#' @param z a vector of covariate related to the linear part.
#' @param lambda1 a real number which is the tuning parameter for LASSO
#' @param lambda2 a real number which is the tuning parameter for trend filering
#' @param k an integer for order of the bounded total variation, default is 2.
#' @param iter.max a maximum iteration numer for the blockwise coordinate descent (BCD).
#' @param tol1 a tolerance error when the BCD loop break.
#' @param ... used for controlling the others.
#' @author Sang Kyu Lee
#' @return \code{pltf.fit} returns an object of class \code{"list"}, which is a list containing
#' following components:
#' @return
#' \item{beta}{fitted coefficients for the linear term.}
#' \item{g}{fitted values for nonparametric terms.}
#' \item{model.g}{fitted nonparmetric object.}
#' @importFrom glmnet glmnet
#' @import glmgen
#' @export
pltf.fit = function(y, x, z, lambda1, lambda2, k=2,
                    iter.max = 1000,
                    tol1 = 1e-3, 
                    ...){
  n = length(y)
  p = dim(x)[2]
  for(iter in 1:iter.max){
    if(iter == 1){
      g.fit = rep(0, n)
    } else {
      beta.old = beta.new
      g.old.total = g.fit
    }     
    y_temp1 = y - g.fit

    lasso.fit = glmnet::glmnet(x = x, y = y_temp1, lambda = lambda1,
                               standardize = FALSE,
                               intercept = FALSE)
    beta.new = c(as.numeric(lasso.fit$beta)) # for tau = 0.5
    y_temp2 = y - x%*%beta.new
    
    tf.fit = trendfilter(x = z, y = y_temp2,
                                 family = "gaussian", k=k,
                                 lambda = lambda2)
    g.fit = predict(tf.fit, x.new = (z), lambda = lambda2)
    if(iter!=1){
      if(sum((x%*%(beta.new - beta.old) + g.fit - g.old.total)^2) < tol1) break
    } else {
      next
    }
  }
  res = list(beta = beta.new, g = g.fit, model.g = tf.fit)
  res
}



#' Partial Linear Smoothing Spline for a High-dimensional Heterogeneous Data
#'
#' The \code{plss.fit} performs partial linear quantile trend filtering model on the high-dimensional data by combining Quantile LASSO and smoothing spline.
#' @param y a vector of a response variable for the model.
#' @param x a matrix of covariates related to the linear part.
#' @param z a vector of covariate related to the linear part.
#' @param lambda1 a real number which is the tuning parameter for LASSO
#' @param lambda2 a real number which is the tuning parameter for trend filering
#' @param iter.max a maximum iteration numer for the blockwise coordinate descent (BCD).
#' @param tol1 a tolerance error when the BCD loop break.
#' @param ... used for controlling the others.
#' @author Sang Kyu Lee
#' @return \code{plss.fit} returns an object of class \code{"list"}, which is a list containing
#' following components:
#' @return
#' \item{beta}{fitted coefficients for the linear term.}
#' \item{g}{fitted values for nonparametric terms.}
#' \item{model.g}{fitted nonparmetric object.}
#' @importFrom rqPen rq.lasso.fit
#' @import stats
#' @export
plss.fit = function(y, x, z, lambda1, lambda2,
                    iter.max = 1000,
                    tol1 = 1e-3, 
                    ...){
  n = length(y)
  p = dim(x)[2]
  for(iter in 1:iter.max){
    if(iter == 1){
      g.fit = rep(0, n)
    } else {
      beta.old = beta.new
      g.old.total = g.fit
    }     
    y_temp1 = y - g.fit

    lasso.fit = glmnet::glmnet(x = x, y = y_temp1, lambda = lambda1,
                               standardize = FALSE,
                               intercept = FALSE)
    beta.new = c(as.numeric(lasso.fit$beta)) 
    y_temp2 = y - x%*%beta.new
    
    ss.fit = smooth.spline(x = z, y = y_temp2,
                           lambda = lambda2)
    g.fit = predict(ss.fit, x = (z), lambda = lambda2)$y
    if(iter!=1){
      if(sum((x%*%(beta.new - beta.old) + g.fit - g.old.total)^2) < tol1) break
    } else {
      next
    }
  }
  res = list(beta = beta.new, g = g.fit, model.g = ss.fit)
  res
}