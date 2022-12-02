#' Weighted Quantile Regression for Longitudinal Big Data Based on Multi-Block ADMM
#' This is a non-parallel implementation.
#'
#' @param x The design matrix (without intercept)
#' @param y The response vector
#' @param rep The repeat observation number for each subject
#' @param tau The quantile of interest
#' @param intercept Whether to include the intercept into the model
#' @param esttype The method for estimating the regression coefficients (conventional quantile regression ("CQR") or weighted quantile regression ("WQR"))
#' @param warmstart Whether to employ a warm start for the algorithm
#' @param corrtype The working correlation matrix to use when computing the weights, currently support ("exchangeable" and "ar1")
#' @param rhoCQR The augmentation parameter for ADMM when computing CQR estimator
#' @param rhoWQR The augmentation parameter for ADMM when computing WQR estimator
#' @param eps The tolerance parameter for convergence (the WQR-ADMM algorithm for computing CQR or WQR estimator)
#' @param epsw The tolerance parameter for convergence (the Newton-Raphson algorithm for computing the weights)
#' @param maxstep Maximum number of iterations allowed for the WQR-ADMM algorithm
#' @param maxstepw Maximum number of iterations allowed for the Newton-Raphson algorithm
#' @return The coefficient estimation, the number of iterations and the running time for the WQR-ADMM algorithm, and the total time cost for computing the WQR estimator (when the parameter esttype = "WQR")
#' @examples
#' gcov = function(p, rho, type){
#'   if(type == "exchangeable"){
#'     cov = matrix(rho, p, p)
#'     diag(cov) = rep(1, p)
#'   }
#'   else{
#'     cov = diag(p)
#'     for(i in 1:p){
#'       for(j in 1:p){
#'         if(i < j) cov[i,j] = rho^{j-i}
#'         else cov[i,j] = cov[j,i]
#'       }
#'     }
#'   }
#'   cov
#' }
#'
#' N = 10000
#' p = 100
#' n = 10
#' rep = rep(n, N)
#' nsum = sum(rep)
#' d = 0.75*p
#' rho_X = 0.5
#' rho_e = 0.5
#' tau = 0.75
#' set.seed(999)
#' X = matrix(rnorm(nsum*p), nsum, p)
#' cov_X = gcov(p, rho_X, "ar1")
#' X = X%*%chol(cov_X)
#' for(i in 1:d){
#'   X[,i] = pnorm(X[,i])
#' }
#' set.seed(1000)
#' e = matrix(rt(N*n, 3), N, n)
#' cov_e = gcov(n, rho_e, "ar1")
#' e = as.vector(t(e%*%chol(cov_e)))
#' sigma = 0.5
#' e = sigma*e
#' beta0 = rnorm(1)
#' beta = rnorm(p)
#' Y = beta0+X%*%beta+apply(X[,1:d]*e/d, 1, sum)
#' beta_true = c(beta0, quantile(e/d, tau)+beta[1:d], beta[(d+1):p])
#'
#' WQR = WQRADMM(X, Y, rep, tau, TRUE, "WQR", TRUE)
#' beta_WQR = WQR$Estimation_WQR
#' AE_WQR = sum(abs(beta_WQR-beta_true))
#' Time_WQR = WQR$Time_WQR
#' Time_WQRADMM = WQR$Time_total
#' @export
#'

WS <- function(x, y, rep, tau, betahat, intercept, corrtype = "ar1", epsw = 1e-06, maxstepw = 100) {
  .Call('_WQRADMM_WS', PACKAGE = 'WQRADMM', x, y, rep, tau, betahat, intercept, corrtype, epsw, maxstepw)
}

WQRADMM <- function(x, y, rep, tau, intercept, esttype, warmstart = TRUE, corrtype = "ar1", rhoCQR = 1, rhoWQR = 1, eps = 1e-04, epsw = 1e-06, maxstep = 5000, maxstepw = 100) {
  .Call('_WQRADMM_WQRADMM', PACKAGE = 'WQRADMM', x, y, rep, tau, intercept, esttype, warmstart, corrtype, rhoCQR, rhoWQR, eps, epsw, maxstep, maxstepw)
}





