get_consecutive <- function(max, number) {
start <- sample(1:(max - number + 1), 1)
return(start:(start + number - 1))
}
###############################################################################
demean=function(X){
Xname=colnames(X)
p=ncol(X)
if(is.null(p)){
X=X-mean(X)
}else{
Xmean=colMeans(X)
X=sweep(X,2,Xmean,"-")
}
colnames(X)=Xname
return(X)
}
###############################################################################
get_pairwise_interactions <- function(W, Z = NULL, noint_env = NULL) {
W <- as.matrix(W)
n <- nrow(W); p <- ncol(W)
colnames_W <- colnames(W)
if (is.null(colnames_W)) colnames_W <- paste0("W", seq_len(p))
q <- 0L
Z_mat <- NULL
colnames_Z <- character(0)
if (!is.null(Z)) {
Z0 <- as.matrix(Z)
if (nrow(Z0) != n) stop("nrow(Z) must equal nrow(W).")
if (is.null(noint_env)) noint_env <- integer(0)
noint_env <- intersect(unique(as.integer(noint_env)), seq_len(ncol(Z0)))
indz <- setdiff(seq_len(ncol(Z0)), noint_env)
if (length(indz) > 0L) {
Z_mat <- Z0[, indz, drop = FALSE]
q <- ncol(Z_mat)
colnames_Z <- colnames(Z_mat)
if (is.null(colnames_Z)) colnames_Z <- paste0("Z", seq_len(q))
} else {
Z_mat <- NULL
q <- 0L
}
}
# --- W × W ---
ww_cols <- p * (p + 1L) / 2L
WW <- matrix(NA_real_, n, ww_cols)
colnames_WW <- character(ww_cols)
col_idx <- 1L
for (i in seq_len(p)) {
for (j in i:p) {
WW[, col_idx] <- W[, i] * W[, j]
colnames_WW[col_idx] <- paste0(colnames_W[i], "*", colnames_W[j])
col_idx <- col_idx + 1L
}
}
if (q > 0L) {
zw_cols <- q * p
ZW <- matrix(NA_real_, n, zw_cols)
colnames_ZW <- character(zw_cols)
col_idx <- 1L
for (i in seq_len(q)) {
for (j in seq_len(p)) {
ZW[, col_idx] <- Z_mat[, i] * W[, j]
colnames_ZW[col_idx] <- paste0(colnames_Z[i], "*", colnames_W[j])
col_idx <- col_idx + 1L
}
}
out <- cbind(WW, ZW)
colnames(out) <- c(colnames_WW, colnames_ZW)
} else {
out <- WW
colnames(out) <- colnames_WW
}
out
}

###############################################################################
get_active_indices <- function(fit) {
cs = tryCatch(summary(fit), error = function(e) NULL)
if (!is.null(cs) && length(cs$cs) > 0) {
active_idx = unique(unlist(cs$cs$cs))
}else{
active_idx=NULL
}
return(active_idx)
}
###############################################################################
Identifying_MainEffect=function(fit,nam){
summ=summary(fit)$vars
g=unique(summ$cs[which(summ$cs>0)])
bb=summary(fit)$cs
S=list()
for(i in g){
indi=which(summ$cs==i)
a=summ$variable[indi]
b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Main_CS",i),log10BF=bb$cs_log10bf[bb$cs==i],PIP=summ$variable_prob[indi])
S[[i]]=b
}
return(do.call(rbind,S))
}
###############################################################################
Identifying_EnvEffect=function(fit,nam){
  summ=summary(fit)$vars
  g=unique(summ$cs[which(summ$cs>0)])
  bb=summary(fit)$cs
  S=list()
  for(i in g){
    indi=which(summ$cs==i)
    a=summ$variable[indi]
    b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Env_CS",i),log10BF=bb$cs_log10bf[bb$cs==i],PIP=summ$variable_prob[indi])
    S[[i]]=b
  }
  return(do.call(rbind,S))
}
###############################################################################
Identifying_IntEffect=function(fitW,namW){
summ=summary(fitW)$vars
if(length(which(summ$cs>0))>0){
bb=summary(fitW)$cs
g=unique(summ$cs[which(summ$cs>0)])
S=list()
for(i in g){
indi=which(summ$cs==i)
a=summ$variable[indi]
b=data.frame(Index=a,Variable=namW[summ$variable[indi]],CS=paste0("Int_CS",i),log10BF=bb$cs_log10bf[bb$cs==i],PIP=summ$variable_prob[indi])
S[[i]]=b
}
return(do.call(rbind,S))
}else{
return(NULL)
}
}

ProjectRes=function(A,B,inercept=F,n_threads){
if(inercept==T){
B=cbind(1,B)
}
BtB = blockwise_crossprod(X=B,n_threads=n_threads)
BtA = blockwise_crossprod(B,A,n_threads)
ProjPart = matrixMultiply(B,(solve(BtB)%*%(BtA)))
return(ProjPart)
}

safe_add_p <- function(idx, Coefmat) {
  if (is.null(idx)) return(NULL)
  if (is.data.frame(idx) && nrow(idx) == 0) return(idx)
  if (!("CS" %in% names(idx))) return(idx)

  cs <- as.character(idx$CS)
  pos <- match(cs, rownames(Coefmat))  # 不存在的返回 NA
  p   <- rep(NA_real_, length(cs))
  p[!is.na(pos)] <- Coefmat[pos[!is.na(pos)], 4]

  idx$Pvalue <- p
  idx
}
