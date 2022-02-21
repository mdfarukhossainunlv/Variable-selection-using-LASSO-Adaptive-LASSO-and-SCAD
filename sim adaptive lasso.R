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
x=as.matrix(x)
X = as.matrix(x)

lambdas <- 10^seq(1, -3, by = -.1)

library(glmnet)

coeffsel.lasso=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.lasso) = paste0("X",seq_len(p))
weight=matrix(0,nrow=Rep,ncol=p)
Betaols = matrix(0, nrow = Rep, ncol = p)
for (r in 1:Rep) {
  Betaols[r,] = solve(t(x)%*%x) %*% t(x) %*% y[,r]
}

Betaols


coeffsel.adaplasso=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.adaplasso) = paste0("X",seq_len(p))
for (r in seq_len(Rep)){
weight[r,] = 1/(abs(Betaols[r,])+.00001)
adaplasso <- cv.glmnet(X, y = y[,r], alpha = 1, 
                   family = "gaussian", intercept=FALSE,lambda=lambdas,standardize=FALSE,
                   nfolds = 5, penalty.factor=weight)
lambda_best1 <- adaplasso$lambda.min
adaplasso_model <- glmnet(X, y[,r], intercept=FALSE, alpha = 1, lambda = lambda_best1,penalty.factor=weight, standardize = FALSE)

coeffsel.adaplasso[r,] = as.vector(unlist(coefficients(adaplasso_model)))[-1]
}
coeffsel.adaplasso
