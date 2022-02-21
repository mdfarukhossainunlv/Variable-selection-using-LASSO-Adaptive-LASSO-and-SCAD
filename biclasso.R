library(leaps)
library(MASS)
set.seed(1980)
n=200
p=10
sigmaE=1.5
rho=.75
betaEst=c(1,-1,2,0,0,0,-1.5,0,1,0)
Rep=500
mux=rep(0,p)
sigmax=matrix(NA,p,p)
for ( i in seq_len(p)) {
  for (j in seq_len(p)){
    sigmax[i,j]=rho^abs(i-j)
  }
}
sigmax
x=mvrnorm(n,mux,sigmax)
y=matrix(0,nrow=n,ncol=Rep)
dim(y)
for (r in seq_len(Rep)){
  y[,r]=x%*%betaEst+rnorm(n)*sigmaE
}
y[,1]
dim(y)
dim(x)
x=data.frame(x)
attach(x)
## lasso
XX = as.matrix(x)
X=as.matrix(x)
lambdas_to_try<- 10^seq(1, -3, by = -.1)

library(glmnet)
X_scaled <- scale(XX)
aic <- c()
bic <- c()
for (r in seq_len(Rep)) {
for (lambda in seq(lambdas_to_try)) {
  # Run model
  model <- glmnet(X, y[,r], alpha = 0, lambda = lambdas_to_try[lambda], standardize = TRUE)
  # Extract coefficients and residuals (remove first row for the intercept)
  betas <- as.vector((as.matrix(coef(model))[-1, ]))
  resid <- y[,r] - (X_scaled %*% betas)
  # Compute hat-matrix and degrees of freedom
  ld <- lambdas_to_try[lambda] * diag(ncol(X_scaled))
  H <- X_scaled %*% solve(t(X_scaled) %*% X_scaled + ld) %*% t(X_scaled)
  df <- tr(H)
  # Compute information criteria
  aic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df
  bic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df * log(nrow(X_scaled))
  lambda_aic <- lambdas_to_try[which.min(aic)]
  lambda_bic <- lambdas_to_try[which.min(bic)]
  lasso_model <- glmnet(X, y[,r], alpha = 1, lambda = lambda_bic, standardize = TRUE)
  coeffsel.lasso[r,] = as.vector(unlist(coefficients(lasso_model)))[-1]
   }
     }





