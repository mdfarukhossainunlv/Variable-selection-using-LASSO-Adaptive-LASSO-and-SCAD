library(leaps)
library(MASS)
set.seed(1980)
n=200
p=10
sigmaE=1.5
sigmae=1.5
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

### backward for repetation
coeffsel=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel) = paste0("X",seq_len(p))

for (j in seq_len(Rep)){
  temp <- stepAIC(lm(y[,j]~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10-1),y[,j]~1, direction ="backward", trace = FALSE) 
  temp1 = c(temp$coefficients)
  coeffsel[j,colnames(coeffsel) %in%names(temp1)] = temp1
}

## SE,SP
tPos = coeffsel[,betaEst!=0] != 0 

fPos = coeffsel[,betaEst==0] != 0 

tNeg = coeffsel[,betaEst==0] == 0 
fNeg = coeffsel[,betaEst!=0] == 0
SEn = apply(tPos, 1, sum)/(apply(tPos, 1, sum)+ apply(fNeg, 1, sum))
mean(SEn)*100

SPe = apply(tNeg, 1, sum)/(apply(tNeg, 1, sum)+ apply(fPos, 1, sum))
mean(SPe)*100

## MSE of beta and sigma
beta.matrix=as.matrix(coeffsel)
dim(beta.matrix)
beta.mean=apply(beta.matrix,2,mean)
estimation.bias <- beta.mean - betaEst
beta.mse<- apply(beta.matrix,2,var)
beta.mse


si=sigma(temp)
si.bias <- si - sigmae
sigma.mse=(si.bias^2)
sigma.mse


##LASSO

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
    model <- glmnet(X, y[,r], alpha = 1, lambda = lambdas_to_try[lambda], standardize = TRUE)
    # Extract coefficients and residuals (remove first row for the intercept)
    betas <- as.vector((as.matrix(coef(model))[-1, ]))
    resid <- y[,r] - (X_scaled %*% betas)
    
    # Compute information criteria
    
    bic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * log(nrow(X_scaled))
    
    lambda_bic <- lambdas_to_try[which.min(bic)]
    lasso_model <- glmnet(X, y[,r], alpha = 1, lambda = lambda_bic, standardize = TRUE)
    coeffsel.lasso[r,] = as.vector(unlist(coefficients(lasso_model)))[-1]
  }
}
## SE,SP

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


Betaols = matrix(0, nrow = Rep, ncol = p)
for (r in 1:Rep) {
  Betaols[r,] = solve(t(x)%*%x) %*% t(x) %*% y[,r]
} 
lambda1 =lambda_bic

coeffsel.lasso = -1*(Betaols<0)
coeffsel.lasso[coeffsel.lasso==0] = 1

temp1 = abs(Betaols) -lambda1
temp1[temp1<0] = 0 
coeffsel.lasso = coeffsel.lasso*temp1
 
tPos = coeffsel.lasso[,betaEst!=0] != 0 

fPos = coeffsel.lasso[,betaEst==0] != 0 

tNeg = coeffsel.lasso[,betaEst==0] == 0 
fNeg = coeffsel.lasso[,betaEst!=0] == 0
SEn = apply(tPos, 1, sum)/(apply(tPos, 1, sum)+ apply(fNeg, 1, sum))
mean(SEn)*100

SPe = apply(tNeg, 1, sum)/(apply(tNeg, 1, sum)+ apply(fPos, 1, sum))
mean(SPe)*100

## MSE of beta and sigma
beta.matrix=as.matrix(coeffsel.lasso)
dim(beta.matrix)
beta.mean=apply(beta.matrix,2,mean)
estimation.bias <- beta.mean - betaEst
beta.mse<- apply(beta.matrix,2,var)
beta.mse

## not correct
si=sigma(temp)
si.bias <- si - sigmae
sigma.mse=(si.bias^2)
sigma.mse

 

##SCAD ??????

library(ncvreg)
library(parallel)
library(lassoshooting)

XX = as.matrix(x)

lambdas <- 10^seq(1, -3, by = -.1)
coeffsel.scad=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.scad) = paste0("X",seq_len(p))

for (r in seq_len(Rep)) {
  cv.out=cv.ncvreg (X,y[,r],intercept=FALSE,lambda=lambdas_to_try,alpha =1,penalty=c("SCAD"))
  bestlam=cv.out$lambda.min
  scad.reg =ncvreg(X,y[,r],alpha=1,lambda=bestlam)
  coeffsel.scad[r,] = as.vector(unlist(coefficients(scad.reg)))[-1]
}
coeffsel.scad[,3] 

lasso_reg <- cv.glmnet(XX, y[,r], intercept=FALSE, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 10)
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  lasso_model <- glmnet(XX, y[,r], intercept=FALSE, alpha = 1, lambda = lambda_best, standardize = TRUE)
  
  coeffsel.lasso[r,] = as.vector(unlist(coefficients(lasso_model)))[-1]





model$beta












###scad

r=150
fit=ncvreg(X, y[,r], penalty="SCAD")
plot(fit)
cvfit <- cv.ncvreg(X,y[,r],penalty="SCAD")
plot(cvfit)
fitall <- cvfit$fit
beta <- fitall$beta[,cvfit$min]
beta

fit <- ncvreg(X, y[,r], penalty="SCAD")
fit$beta



ll <- log(fit$lambda)
op <- par(mfrow=c(2,1))
plot(ll, BIC(fit), type="l", xlim=rev(range(ll)))
lam <- fit$lambda[which.min(BIC(fit))]
b <- coef(fit, lambda=lam)
b[b!=0]
plot(fit)
abline(v=lam)
par(op)








r=2

cvfit <- cv.ncvreg(X, y[,r],  nfolds=5)
cvfit$lambda.min
lambdas_to_try
fit=ncvreg(X, y[,r], family="gaussian", penalty="SCAD", gamma=3.7,
           alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=300, lambda=lambdas_to_try, eps=1e-4)

fit$beta