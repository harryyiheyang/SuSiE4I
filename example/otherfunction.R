glinternet_evaulate=function(true_main_index,true_int_variable,fit_glinternet){
i_1Std <- which(fit_glinternet$lambdaHat1Std == fit_glinternet$lambda)
coefs <- coef(fit_glinternet$glinternetFit)[[i_1Std]]
tp_main=ifelse(length(setdiff(true_main_index,coefs$mainEffects$cont))==0,1,0)
tn_main=ifelse(length(setdiff(coefs$mainEffects$cont,true_main_index))==0,1,0)
if(is.null(coefs$interactions$contcont)){
tp_int=ifelse(is.null(true_int_variable),1,0)
tn_int=1
}else{
p=nrow(coefs$interactions$contcont)
int_variable=c()
for(i in 1:p){
int_variable[i]=paste0("X",coefs$interactions$contcont[i,1],"X",coefs$interactions$contcont[i,2])
}
tp_int=ifelse(length(setdiff(true_int_variable,int_variable))==0,1,0)
tn_int=ifelse(length(setdiff(int_variable,true_int_variable))==0,1,0)
}
g=data.frame(tp_main=tp_main,tn_main=tn_main,tp_int=tp_int,tn_int=tn_int)
return(g)
}


RAMP_evaulate=function(true_main_index,true_int_variable,fit_RAMP){
tp_main=ifelse(length(setdiff(true_main_index,fit_RAMP$mainInd))==0,1,0)
tn_main=ifelse(length(setdiff(fit_RAMP$mainInd,true_main_index))==0,1,0)
tp_int=ifelse(length(setdiff(true_int_variable,fit_RAMP$interInd))==0,1,0)
tn_int=ifelse(length(setdiff(true_int_variable,fit_RAMP$interInd))==0,1,0)
g=data.frame(tp_main=tp_main,tn_main=tn_main,tp_int=tp_int,tn_int=tn_int)
return(g)
}

SuSiE4X_evaulate=function(true_main_index,true_int_variable,fit_SuSiE4X){
int=int_tptn(true_int_variable, fit_SuSiE4X$main_index, fit_SuSiE4X$interaction_index)
main=main_tptn(true_main_index,fit_SuSiE4X$main_index)
g=data.frame(tp_main=main$tp,tn_main=main$tn,tp_int=int$tp,tn_int=int$tn)
return(g)
}

################################################################################
int_tptn <- function(true_int_variable, main_index, int_index) {
if (is.null(int_index)) {
tn <- 1
tp <- ifelse(is.null(true_int_variable), 1, 0)
return(list(tp = tp, tn = tn))
}
# Step 1: mappings
cs_to_var <- split(main_index$Variable, main_index$CS)
var_to_cs <- setNames(main_index$CS, main_index$Variable)
# Step 2: Expand each interaction CS to variable-level pairs
expand_interaction_index <- function(pair_string) {
parts <- strsplit(pair_string, "\\*")[[1]]
part1 <- parts[1]
part2 <- parts[2]
vars1 <- if (grepl("^Main_CS", part1)) cs_to_var[[part1]] else part1
vars2 <- if (grepl("^Main_CS", part2)) cs_to_var[[part2]] else part2
if (is.null(vars1) || is.null(vars2)) return(character(0))
expand_grid <- expand.grid(a = vars1, b = vars2, stringsAsFactors = FALSE)
pairs <- unique(c(
paste0(expand_grid$a, "*", expand_grid$b),
paste0(expand_grid$b, "*", expand_grid$a)
))
return(pairs)
}
# Step 3: Expand all detected interaction pairs
all_detected <- unlist(lapply(int_index$Variable, expand_interaction_index))
# Step 4: Normalize all interactions as sorted(a,b)
normalize_pair <- function(x) paste(sort(strsplit(x, "\\*")[[1]]), collapse = "*")
detected_norm <- unique(sapply(all_detected, normalize_pair))
true_norm <- unique(sapply(true_int_variable, normalize_pair))
tp <- as.integer(all(true_norm %in% detected_norm))
# Step 5: TN — compute CS pairs corresponding to true interactions
get_cs_pair <- function(pair) {
parts <- strsplit(pair, "\\*")[[1]]
cs1 <- if (parts[1] %in% names(var_to_cs)) var_to_cs[parts[1]] else parts[1]
cs2 <- if (parts[2] %in% names(var_to_cs)) var_to_cs[parts[2]] else parts[2]
paste(sort(c(cs1, cs2)), collapse = "*")
}
true_cs_pairs <- unique(sapply(true_int_variable, get_cs_pair))
detected_cs_pairs <- unique(sapply(int_index$Variable, function(x) {
paste(sort(strsplit(x, "\\*")[[1]]), collapse = "*")
}))
tn <- as.integer(length(setdiff(detected_cs_pairs, true_cs_pairs)) == 0)
return(list(tp = tp, tn = tn))
}


main_tptn = function(true_main_index, main_index) {
if (is.null(true_main_index)) {
tp = tn = ifelse(is.null(main_index), 1, 0)
return(list(tp=tp, tn=tn))
}
# TP: all true main effects must be covered
tp = as.integer(all(true_main_index %in% main_index$Index))
# TN: check whether any extra credible sets (not involving true main effects) are identified
true_main_cs = unique(main_index$CS[main_index$Index %in% true_main_index])
all_cs = unique(main_index$CS)
fp_cs = setdiff(all_cs, true_main_cs)
tn = as.integer(length(fp_cs) == 0)
return(list(tp=tp, tn=tn))
}
