Run_GGE_Select <- function(X, Z, y,crossprodX=NULL, Lmain,Lenv, Lint, noint_env=NULL,max.iter, min.iter, max.eps, susie.iter,verbose=T,n_threads=1,coverage.main = 0.95, coverage.env=0.95, coverage.int = 0.85,returnModel=F,...) {

n <- length(y)

Xmean=colMeans(X)
etaX=etaZ=etaW=0
meanY=mean(y)
ZtZ = blockwise_crossprod(X=Z,n_threads=n_threads)
if(is.null(crossprodX)){
XtX = blockwise_crossprod(X=X,n_threads=n_threads)
}else{
XtX=crossprodX
}

g=c()
beta=0
alpha=0
gamma=0
meanY=mean(y)

for (iter in 1:max.iter) {
gamma_prev <- gamma
beta_prev <- beta
alpha_prev <- alpha

## --- update etaZ ---
rZ <- y - etaW - etaX - meanY
Zty <- crossprod(Z, rZ)
yty4Z=sum(rZ^2)
fitZ <- susie_suff_stat(XtX = ZtZ, Xty = Zty, yty = yty4Z, n = n, L = Lenv,
                    max_iter = susie.iter, estimate_prior_method = "EM",
                    X_colmeans=colMeans(Z),y_mean=mean(rZ),coverage=coverage.env,...)
alpha=coef(fitZ)[-1]
etaZ = matrixVectorMultiply(Z,alpha)
# Extract credible sets using summary information
CSdt <- summary(fitZ)$vars
cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
if(length(cs_indices) == 0) {
ZCS <- NULL
} else {
Alpha_filtered <- fitZ$alpha * 0
for(i in cs_indices) {
vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
if(length(vars_in_cs_i) > 0) {
Alpha_filtered[i, vars_in_cs_i] <- fitZ$alpha[i, vars_in_cs_i]
}
}
ZCS <- matrixMultiply(Z, t(as.matrix(Alpha_filtered)))
ZCS <- ZCS[, cs_indices, drop = FALSE]
if(is.null(dim(ZCS))) {
ZCS <- matrix(ZCS, ncol = 1)
}
colnames(ZCS) <- paste0("Env_CS", cs_indices)
ZCS <- as.matrix(ZCS)
}

## --- update etaX ---
rX <- y - etaW - etaZ - meanY
Xty <- as.vector(crossprod(X, rX))
yty4X <- sum(rX^2)
fitX <- susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty4X, n = n, L = Lmain,
                    max_iter = susie.iter, estimate_prior_method = "EM",
                    X_colmeans=Xmean,y_mean=mean(rX),coverage=coverage.main,...)
beta <- coef.susie(fitX)[-1]
etaX = matrixVectorMultiply(X,beta)#+coef.susie(fitX)[1]

## --- update etaW ---
rW <- y - etaX - etaZ - meanY
# Extract credible sets using summary information
CSdt <- summary(fitX)$vars
cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
if(length(cs_indices) == 0) {
warning("No credible set detected at iteration ", iter)
break
}
Alpha_filtered <- fitX$alpha * 0
for(i in cs_indices) {
vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
}
XCS <- matrixMultiply(X, t(as.matrix(Alpha_filtered)))
XCS <- XCS[, cs_indices, drop = FALSE]
if(is.null(dim(XCS))) {
XCS <- matrix(XCS, ncol = 1)
}
colnames(XCS) <- paste0("Main_CS", cs_indices)
XCS <- as.matrix(XCS)

G = cbind(1,ZCS,XCS)
W <- get_pairwise_interactions(XCS, ZCS)
#Wres=W-ProjectRes(A=W,B=G,inercept=F,n_threads=n_threads)
WtW <- blockwise_crossprod(X=W,n_threads=n_threads)
Wty <- as.vector(crossprod(W, rW))
yty4W <- sum(rW^2)
fitW <- susie_suff_stat(XtX = WtW, Xty = Wty, yty = yty4W, n = n, L = Lint,
                        max_iter = susie.iter, estimate_prior_method = "EM",
                        coverage=coverage.int,...)

CSdt <- summary(fitW)$vars
cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
if(length(cs_indices) == 0) {
WCS <- NULL
} else {
Alpha_filtered <- fitW$alpha * 0
for(i in cs_indices) {
vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
if(length(vars_in_cs_i) > 0) {
Alpha_filtered[i, vars_in_cs_i] <- fitW$alpha[i, vars_in_cs_i]
}
}
WCS <- matrixMultiply(W, t(as.matrix(Alpha_filtered)))
WCS <- WCS[, cs_indices, drop = FALSE]
if(is.null(dim(WCS))) {
WCS <- matrix(WCS, ncol = 1)
}
colnames(WCS) <- paste0("Int_CS", cs_indices)
WCS <- as.matrix(WCS)
}

if(is.null(ZCS)==F){
if(is.null(WCS)==F){
fit_final=lm(y~ZCS+XCS+WCS,model=F,x=F,y=F)
coefs=coef(fit_final)
meanY=coefs[1]
ZCSalpha=coefs[2:(ncol(ZCS)+1)]
XCSbeta=coefs[(ncol(ZCS)+2):(ncol(ZCS)+ncol(XCS)+1)]
WCSbeta=coefs[-c(1:(ncol(ZCS)+ncol(XCS)+1))]
etaZ=matrixVectorMultiply(ZCS,ZCSalpha)
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=matrixVectorMultiply(WCS,WCSbeta)
}else{
fit_final=lm(y~ZCS+XCS,model=F,x=F,y=F)
coefs=coef(fit_final)
meanY=coefs[1]
ZCSalpha=coefs[2:(ncol(ZCS)+1)]
XCSbeta=coefs[(ncol(ZCS)+2):(ncol(ZCS)+ncol(XCS)+1)]
etaZ=matrixVectorMultiply(ZCS,ZCSalpha)
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=0
}
}

if(is.null(ZCS)==T){
if(is.null(WCS)==F){
fit_final=lm(y~XCS+WCS,model=F,x=F,y=F)
coefs=coef(fit_final)
meanY=coefs[1]
XCSbeta=coefs[2:(ncol(XCS)+1)]
WCSbeta=coefs[-c(1:(ncol(XCS)+1))]
etaZ=0
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=matrixVectorMultiply(WCS,WCSbeta)
}else{
fit_final=lm(y~XCS,model=F,x=F,y=F)
coefs=coef(fit_final)
meanY=coefs[1]
XCSbeta=coefs[-1]
etaZ=0
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=0
}
}
## --- check convergence ---
errX <- sqrt(mean((beta - beta_prev)^2))
errZ <- sqrt(mean((alpha - alpha_prev)^2))
err = c(errX, errZ)
g[iter]=max(err)
if (verbose == TRUE){
cat(paste0("Iteration ",iter,"\n"))
}
if (max(err) < max.eps&iter>min.iter) {
break
}
}

if(is.null(ZCS)==F){
if(is.null(WCS)==F){
Dat=as.data.frame(cbind(y,ZCS,XCS,WCS))
fit_final=lm(y~.,data=Dat,model=F,x=F,y=F)
}else{
Dat=as.data.frame(cbind(y,ZCS,XCS))
fit_final=lm(y~.,data=Dat,model=F,x=F,y=F)
}
}

if(is.null(ZCS)==T){
if(is.null(WCS)==F){
Dat=as.data.frame(cbind(y,XCS,WCS))
fit_final=lm(y~.,data=Dat,model=F,x=F,y=F)
}else{
Dat=as.data.frame(cbind(y,XCS))
fit_final=lm(y~.,data=Dat,model=F,x=F,y=F)
}
}

## --- post-hoc: linear model of Z with offset ---
MainIndex=Identifying_MainEffect(fitX,colnames(X))
EnvIndex=Identifying_EnvEffect(fitZ,colnames(Z))
IntIndex=Identifying_IntEffect(fitW,colnames(W))

G=summary(fit_final)$coefficient
MainIndex <- safe_add_p(MainIndex, G)
EnvIndex  <- safe_add_p(EnvIndex,  G)
IntIndex  <- safe_add_p(IntIndex,  G)


if (verbose == TRUE) {
plot(g,
type = "o",
col = "black",
pch = 16,
xlab = "Iteration",
ylab = "Max Parameter Change",
main = "Convergence Trace of SuSiE4X (Max |Δ| in alpha and beta)")

for (i in seq_along(g)) {
text(x = i,
y = g[i],
labels = formatC(g[i], format = "e", digits = 1),
pos = 3,        #
cex = 0.7,      #
col = "red")
}
}
AA=list(iter=iter,
error=g,
fitX = fitX,
fitW = fitW,
fitZ = fitZ,
fitJoint = fit_final,
main_index=MainIndex,
interaction_index=IntIndex,
envrionment_index=EnvIndex,
JointCoef=G)
if(returnModel==T){
  AA$FinalModel=Dat
}
return(AA)
}
