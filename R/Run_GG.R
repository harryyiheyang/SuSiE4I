Run_GG <- function(X, y,crossprodX=NULL, Lmain, Lint, max.iter, min.iter, max.eps, susie.iter,verbose=T,n_threads=1,coverage.main = 0.95, coverage.int = 0.85,returnModel=F,...) {

n <- length(y)

Xmean=colMeans(X)
etaX=etaZ=etaW=0
meanY=mean(y)
if(is.null(crossprodX)){
XtX = blockwise_crossprod(X=X,n_threads=n_threads)
}else{
XtX=crossprodX
}

g=c()
beta=0
gamma=0

for (iter in 1:max.iter) {
gamma_prev <- gamma
beta_prev <- beta

## --- update etaX ---
rX <- y - etaW - meanY
Xty <- as.vector(crossprod(X, rX))
yty4X <- sum(rX^2)
fitX <- susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty4X, n = n, L = Lmain,
                  max_iter = susie.iter, estimate_prior_method = "EM",
                  X_colmeans=Xmean,y_mean=mean(rX),coverage=coverage.main,...)
beta <- coef.susie(fitX)[-1]
etaX = matrixVectorMultiply(X,beta)+coef.susie(fitX)[1]

## --- update etaW ---
rW <- y - etaX - meanY
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
G = cbind(1,XCS)
W <- get_pairwise_interactions(XCS)
#Wres=W-ProjectRes(A=W,B=G,inercept=F,n_threads=n_threads)
WtW <- blockwise_crossprod(X=W,n_threads=n_threads)
Wty <- as.vector(crossprod(W, rW))
yty4W <- sum(rW^2)
fitW <- susie_suff_stat(XtX = WtW, Xty = Wty, yty = yty4W, n = n, L = Lint,
                  max_iter = susie.iter, estimate_prior_method = "EM",
                  coverage=coverage.int,...)

## Update WCS
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

if(is.null(WCS)==F){
fit_final=lm(y~XCS+WCS,model=F,x=F,y=F)
coefs=coef(fit_final)
meanY=coefs[1]
XCSbeta=coefs[2:(ncol(XCS)+1)]
WCSbeta=coefs[-c(1:(ncol(XCS)+1))]
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=matrixVectorMultiply(WCS,WCSbeta)
}else{
fit_final=lm(y~XCS,model=F,x=F,y=F)
coefs=coef(fit_final)
meanY=coefs[1]
XCSbeta=coefs[-1]
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=0
}

## --- check convergence ---
err <- sqrt(mean((beta - beta_prev)^2))
g[iter]=max(err)
if (verbose == TRUE){
cat(paste0("Iteration ",iter,"\n"))
}
if (max(err) < max.eps&iter>min.iter) {
break
}
}

if(is.null(WCS)==F){
Dat=as.data.frame(cbind(y,XCS,WCS))
fit_final=lm(y~.,data=Dat,model=F,x=F,y=F)
}else{
Dat=as.data.frame(cbind(y,XCS))
fit_final=lm(y~.,data=Dat,model=F,x=F,y=F)
}

MainIndex=Identifying_MainEffect(fitX,colnames(X))
IntIndex=Identifying_IntEffect(fitW,colnames(W))
G=summary(fit_final)$coefficient
MainIndex <- safe_add_p(MainIndex, G)
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
fitJoint = fit_final,
main_index=MainIndex,
interaction_index=IntIndex,
JointCoef=G)
if(returnModel==T){
AA$FinalModel=Dat
}
return(AA)
}
