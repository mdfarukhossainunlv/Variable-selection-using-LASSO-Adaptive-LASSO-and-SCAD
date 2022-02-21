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
attach(x)
###scad
coeffsel.scad=matrix(0,nrow=Rep,ncol=p)
colnames(coeffsel.scad) = paste0("X",seq_len(p))
lambdas <- 10^seq(1, -3, by = -.1)
for (i in seq_len(Rep)){
cv.out=cv.ncvreg (x,y[,i],alpha =1,penalty=c("SCAD"))
bestlam=cv.out$lambda.min
model =ncvreg(x,y[,i],alpha=1,penalty=c("SCAD"),lambda=bestlam)
coeffsel.scad[i,]=as.vector(model$beta)[-1]
}


coeffsel.scad[1,]

coeffsel.scad




