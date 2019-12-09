#' Weighted Quantile Regression for Longitudinal Data Based on Multi-Block ADMM
#' This is a non-parallel implementation.
#'
#' @param x The design matrix (without intercept)
#' @param y The response vector
#' @param rep The repeat observation number for each subject
#' @param tau The quantile of interest
#' @param intercept Whether to include the intercept into the model
#' @param esttype The method for estimating the regression coefficients (conventional quantile regression ("CQR") or weighted quantile regression ("WQR"))
#' @param corrtype The working correlation matrix to use when computing the weights, currently support ("ar1" and "exchangeable")
#' @param rhoCQR The augmentation parameter for ADMM when computing CQR estimator
#' @param rhoWQR The augmentation parameter for ADMM when computing WQR estimator
#' @param eps The tolerance parameter for convergence (the WQR-ADMM algorithm for computing CQR or WQR estimator)
#' @param epsw The tolerance parameter for convergence (the Newton-Raphson algorithm for computing the weights)
#' @param maxstep Maximum number of iterations allowed for the WQR-ADMM algorithm
#' @param maxstepw Maximum number of iterations allowed for the Newton-Raphson algorithm
#' @return The coefficient estimation, the number of iterations and the running time for the WQR-ADMM algorithm, and the total time cost for computing the WQR estimator (when the parameter esttype = "WQR")
#' @examples
#' N = 10000
#' p = 100
#' n = 10
#' rep = rep(n, N)
#' nsum = sum(rep)
#' rho_X = 0.5
#' rho_e = 0.5
#' tau = 0.7
#'
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
#' set.seed(66)
#' X = matrix(rnorm(nsum*p), nsum, p)
#' cov_X = gcov(p, rho_X, "ar1")
#' X = X%*%chol(cov_X)
#' X[,1] = pnorm(X[,1])
#' beta_true = rnorm(p)
#' beta_true[1] = qnorm(tau)
#' cov_e = gcov(n, rho_e, "ar1")
#' e = matrix(rnorm(N*n), N, n)
#' e = as.vector(t(e%*%chol(cov_e)))
#' Y = X[,2:p]%*%beta_true[2:p]+X[,1]*e
#'
#' CQR = WQRADMM(X, Y, rep, tau, FALSE, "CQR", "ar1")
#' beta_CQR = CQR$Estimation_CQR
#' AE_CQR = sum(abs(beta_CQR-beta_true))
#' Iteration_CQR = CQR$Iteration_CQR
#' Time_CQR = CQR$Time_CQR
#'
#' WQR = WQRADMM(X, Y, rep, tau, FALSE, "WQR", "ar1")
#' beta_WQR = WQR$Estimation_WQR
#' AE_WQR = sum(abs(beta_WQR-beta_true))
#' Iteration_WQR = WQR$Iteration_WQR
#' Time_WQR = WQR$Time_WQR
#' Time_total = WQR$Time_total
#' @export
#'

WQRADMM <- function(x, y, rep, tau, intercept, esttype, corrtype, rhoCQR = 5, rhoWQR = 5, eps = 1e-03, epsw = 1e-06, maxstep = 1000, maxstepw = 100) {
  .Call('_WQRADMM_WQRADMMCPP', PACKAGE = 'WQRADMM', x, y, rep, tau, intercept, esttype, corrtype, rhoCQR, rhoWQR, eps, epsw, maxstep, maxstepw)
}




