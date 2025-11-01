library(CppMatrix)
library(susieR)
library(devtools)
library(MASS)
library(RAMP)
library(glinternet)
load_all()
ARcov=function(p,rho){
s=c(1:p)
for(i in 1:p){
s[i]=rho^(i-1)
}
return(toeplitz(s))
}
CScov=function(p,rho){
diag(p)*(1-rho)+matrix(rho,p,p)
}
source("example/otherfunction.R")
n=1000
p=10
Main_TP=Main_TN=Int_TP=Int_TN=TIME=matrix(0,100,3)

for(iter in 1:100){
###############################################################################
R=kronecker(CScov(5,0.3),ARcov(2,0.5))
X=mvrnorm(n=n,mu=runif(p,0,1),R)
Xmean=colMeans(X)
Z=mvrnorm(n=n,mu=runif(15,0,1),diag(15))
colnames(Z)=paste0("UKBB",1:15)
alpha0=c(1:15)*0
alpha0[1:5]=rnorm(5,0,1/sqrt(5))
beta0=rep(0,p)
ind=sample(p,3)
ind=sort(ind)
beta0[ind]=0.5
eta=matrixVectorMultiply(Z,alpha0)+matrixVectorMultiply(X,beta0)+(X[,ind[1]]*X[,ind[2]]+X[,ind[1]]*Z[,1])/2
true_int_variable1=c(paste0("X",ind[1],"*X",ind[2]),paste0("X",ind[1],"*X",ind[3]))
true_int_variable2=c(paste0("X",ind[1],"X",ind[2]),paste0("X",ind[1],"X",ind[3]))
y=eta+rnorm(n,0,2)
Lmain=5;Linteraction=5;max.iter = 50;min.iter = 3;max.eps = 1e-5;susie.iter = 300;verbose=T;crossprodX=NULL;n_threads=10

t1=Sys.time()
fit_SuSiE4I=SuSiE4I(select_env=F,noint_env=11:15,X=X,Z=Z,y=y,L_main=5,L_int=5,n_threads=10,max_iter=15,verbose=T,susie_iter=1000,coverage_main=0.8)
t2=Sys.time()
t_SuSiE4I=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_RAMP=RAMP(X=cbind(X,Z),y,penalty="MCP",n.lambda=50)
t2=Sys.time()
t_RAMP=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_glinternet=glinternet.cv(X=cbind(X,Z),Y=y,numLevels=rep(1,ncol(X)+ncol(Z)),nLambda=50,numCores=4,verbose=F)
t2=Sys.time()
t_glinternet=difftime(t2, t1, units = "secs")

tptn_SuSiE4I=SuSiE4I_evaulate(true_main_index=ind,true_int_variable=true_int_variable1,fit_SuSiE4I)
tptn_RAMP=RAMP_evaulate(true_main_index=ind,true_int_variable=true_int_variable2,fit_RAMP)
tptn_glinternet=glinternet_evaulate(true_main_index=ind,true_int_variable=true_int_variable2,fit_glinternet)

TIME[iter,]=c(t_SuSiE4I,t_RAMP,t_glinternet)
Main_TP[iter,]=c(tptn_SuSiE4I$tp_main,tptn_RAMP$tp_main,tptn_glinternet$tp_main)
Main_TN[iter,]=c(tptn_SuSiE4I$tn_main,tptn_RAMP$tn_main,tptn_glinternet$tn_main)
Int_TP[iter,]=c(tptn_SuSiE4I$tp_int,tptn_RAMP$tp_int,tptn_glinternet$tp_int)
Int_TN[iter,]=c(tptn_SuSiE4I$tn_int,tptn_RAMP$tn_int,tptn_glinternet$tn_int)

if(iter%%10==0) print(iter)

}
