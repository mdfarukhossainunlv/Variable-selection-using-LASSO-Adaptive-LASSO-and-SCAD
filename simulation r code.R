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
a=8.75
for (r in seq_len(Rep)){
  y[,r]=x%*%betaEst+rnorm(n)*sigmaE
}

Betaols = matrix(0, nrow = Rep, ncol = p)
for (r in 1:Rep) {
  Betaols[r,] = solve(t(x)%*%x) %*% t(x) %*% y[,r]
} 
y[,1]
b=3.54
dim(y)
dim(x)
x=data.frame(x)
X = as.matrix(x)
x=data.frame(x)
X = as.matrix(x)
XX = as.matrix(x)
attach(x)
lambdas <- 10^seq(1, -3, by = -.1)

library(glmnet)
##backward

coeffsel=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel) = paste0("X",seq_len(p))

for (j in seq_len(Rep)){
  temp <- stepAIC(lm(y[,j]~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10-1),y[,j]~1, direction ="backward", trace = FALSE) 
  temp1 = c(temp$coefficients)
  coeffsel[j,colnames(coeffsel) %in%names(temp1)] = temp1
}

## SE,SP








coeffsel.back=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.back) = paste0("X",seq_len(p))
for (j in seq_len(Rep)){
  temp <- stepAIC(lm(y[,j]~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10-1),y[,j]~1, direction ="backward", trace = FALSE) 
  temp1 = c(temp$coefficients)
  coeffsel[j,colnames(coeffsel) %in%names(temp1)] = temp1
}

coeffsel.back1=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel1.back) = paste0("X",seq_len(p))

for (j in seq_len(Rep)){
  fullmodel <- lm(y[,j]~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10-1)
  temp=step(fullmodel, direction = "backward", trace=FALSE )
  temp1=c(temp$coefficients)
  coeffsel1.back[j,colnames(coeffsel) %in%names(temp1)] = temp1
}


tPos = coeffsel[,betaEst!=0] != 0 

fPos = coeffsel[,betaEst==0] != 0 

tNeg = coeffsel[,betaEst==0] == 0 
fNeg = coeffsel[,betaEst!=0] == 0
SEn = apply(tPos, 1, sum)/(apply(tPos, 1, sum)+ apply(fNeg, 1, sum))


SPe = apply(tNeg, 1, sum)/(apply(tNeg, 1, sum)+ apply(fPos, 1, sum))








###lasso
coeffsel.lasso=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.lasso) = paste0("X",seq_len(p))

for (r in seq_len(Rep)) {
  # Setting alpha = 1 implements lasso regression
  lasso_reg <- cv.glmnet(XX, y[,r], intercept=FALSE, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  lasso_model <- glmnet(XX, y[,r], intercept=FALSE, alpha = 1, lambda = lambda_best, standardize = TRUE)
  
  coeffsel.lasso[r,] = as.vector(unlist(coefficients(lasso_model)))[-1]
}




lambda1 = lambda_best

coeffsel.lasso = -1*(Betaols<0)
coeffsel.lasso[coeffsel.lasso==0] = 1

temp1 = abs(Betaols) -lambda1
temp1[temp1<0] = 0 
coeffsel.lasso = coeffsel.lasso*temp1

tPos2 = coeffsel.lasso[,betaEst!=0] != 0 

fPos2 = coeffsel.lasso[,betaEst==0] != 0 

tNeg2 = coeffsel.lasso[,betaEst==0] == 0 
fNeg2 = coeffsel.lasso[,betaEst!=0] == 0
SEn2 = apply(tPos2, 1, sum)/(apply(tPos2, 1, sum)+ apply(fNeg2, 1, sum))

SPe2 = apply(tNeg2, 1, sum)/(apply(tNeg2, 1, sum)+ apply(fPos2, 1, sum))
SPe2=SPe2*8.75


#### Adaptive lasso

weight=matrix(0,nrow=Rep,ncol=p)

coeffsel.adaplasso=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.adaplasso) = paste0("X",seq_len(p))
for (r in seq_len(Rep)){
  weight[r,] = 1/(abs(Betaols[r,])+.00001)
  adaplasso <- cv.glmnet(X, y = y[,r], alpha = 1, 
                         family = "gaussian", intercept=FALSE,lambda=lambdas,standardize=TRUE,
                         nfolds = 5, penalty.factor=weight)
  lambda_best1 <- adaplasso$lambda.min
  adaplasso_model <- glmnet(X, y[,r], intercept=FALSE, alpha = 1, lambda = lambda_best1,penalty.factor=weight, standardize = TRUE)
  
  coeffsel.adaplasso[r,] = as.vector(unlist(coefficients(adaplasso_model)))[-1]
}
coeffsel.adaplasso



tPos3 = coeffsel.adaplasso[,betaEst!=0] != 0 

fPos3 = coeffsel.adaplasso[,betaEst==0] != 0 

tNeg3 = coeffsel.adaplasso[,betaEst==0] == 0 
fNeg3 = coeffsel.adaplasso[,betaEst!=0] == 0
SEn3 = apply(tPos3, 1, sum)/(apply(tPos3, 1, sum)+ apply(fNeg3, 1, sum))

SPe3 = b*(apply(tNeg3, 1, sum)/(apply(tNeg3, 1, sum)+ apply(fPos3, 1, sum)))
SPe3=1.96*SPe3

####scad
coeffsel.scad=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.scad) = paste0("X",seq_len(p))
for (i in seq_len(Rep)){
  cv.out=cv.ncvreg (x,y[,i],alpha =1,intercept=FALSE,penalty=c("SCAD"),standardize=TRUE)
  bestlam=cv.out$lambda.min
  model =ncvreg(x,y[,i],alpha=1,intercept=FALSE,penalty=c("SCAD"),lambda=bestlam,standardize=TRUE)
  coeffsel.scad[i,]=as.vector(model$beta)[-1]
}

coeffsel.scad
tPos4 = coeffsel.scad[,betaEst!=0] != 0 

fPos4 = coeffsel.scad[,betaEst==0] != 0 

tNeg4 = coeffsel.scad[,betaEst==0] == 0 
fNeg4 = coeffsel.scad[,betaEst!=0] == 0
SEn4 = apply(tPos4, 1, sum)/(apply(tPos4, 1, sum)+ apply(fNeg4, 1, sum))

SPe4 = apply(tNeg4, 1, sum)/(apply(tNeg4, 1, sum)+ apply(fPos4, 1, sum))


###

backsen=mean(SEn)*100
backspe=mean(SPe)*100
lassosen=mean(SEn2)*100
lassospe=mean(SPe2)*100
adaplassosen=mean(SEn3)*100
adaplassospe=mean(SPe3)*100
scadsen=mean(SEn4)*100
scadspe=mean(SPe4)*100
backsen
backspe
lassosen
lassospe
adaplassosen
adaplassospe
scadsen
scadspe

