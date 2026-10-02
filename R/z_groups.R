###############################################################################
# Factor-coded Z columns (haplotypes).
#
# z_groups labels which Z indicator columns belong to the same factor. The
# interaction stage still runs susie_ss on single columns: every level times
# each main CS column, and every level-by-level product of two different
# factors, is its own candidate. Products within a factor are never formed,
# and a column is dropped when too few observations carry its level(s).
###############################################################################

normalize_z_groups <- function(z_groups, Z) {
if (is.null(z_groups)) return(NULL)
q <- ncol(Z)
nm <- colnames(Z)
out <- rep(NA_character_, q)
if (is.list(z_groups)) {
labels <- names(z_groups)
if (is.null(labels)) labels <- rep("", length(z_groups))
labels[!nzchar(labels)] <- paste0("G", which(!nzchar(labels)))
if (anyDuplicated(labels)) stop("z_groups names must be unique.")
for (k in seq_along(z_groups)) {
cols <- z_groups[[k]]
pos <- if (is.character(cols)) match(cols, nm) else as.integer(cols)
if (!length(pos) || any(is.na(pos)) || any(pos < 1L | pos > q)) {
stop("z_groups element '", labels[k], "' names columns not in Z.")
}
if (any(!is.na(out[pos]))) stop("A Z column belongs to more than one z_groups entry.")
out[pos] <- labels[k]
}
} else {
if (length(z_groups) != q) {
stop("z_groups must be a list of Z columns or a vector with one label per Z column.")
}
lab <- as.character(z_groups)
lab[lab %in% c("", "0")] <- NA_character_
out <- lab
}
if (all(is.na(out))) stop("z_groups does not assign any Z column to a group.")
out
}

z_group_support <- function(Z) {
Z <- as.matrix(Z)
out <- matrix(FALSE, nrow(Z), ncol(Z), dimnames = dimnames(Z))
for (j in seq_len(ncol(Z))) {
z <- Z[, j]
tab <- table(round(z, 8))
base <- as.numeric(names(tab)[which.max(tab)])
out[, j] <- abs(z - base) > 1e-8
}
out
}

check_z_group_coding <- function(Z_raw, z_groups) {
for (f in unique(stats::na.omit(z_groups))) {
cols <- which(z_groups == f)
if (length(cols) < 2L) next
s <- rowSums(as.matrix(Z_raw[, cols, drop = FALSE]))
if (stats::sd(s) < 1e-8) {
warning("The Z columns of group '", f, "' sum to a constant and are collinear with the intercept; drop one reference level.",
        call. = FALSE)
}
}
invisible(NULL)
}

get_grouped_pairwise_interactions <- function(W, Z, noint_env = NULL,
                                              include_x_squared = FALSE,
                                              z_groups,
                                              z_support = NULL,
                                              min_group_int_obs = 100L) {
Z <- as.matrix(Z)
n <- nrow(Z)
nmZ <- colnames(Z)
if (is.null(nmZ)) nmZ <- paste0("Z", seq_len(ncol(Z)))
if (length(z_groups) != ncol(Z)) stop("z_groups must have one entry per Z column.")
if (is.null(z_support)) z_support <- z_group_support(Z)
if (!identical(dim(z_support), dim(Z))) stop("z_support must have the same dimensions as Z.")
if (is.null(noint_env)) noint_env <- integer(0)
interacting <- setdiff(seq_len(ncol(Z)), intersect(unique(as.integer(noint_env)), seq_len(ncol(Z))))
ungrouped <- interacting[is.na(z_groups[interacting])]
grouped <- interacting[!is.na(z_groups[interacting])]
factors <- unique(z_groups[grouped])
n_support <- colSums(z_support)

base <- NULL
if (!is.null(W)) {
W <- as.matrix(W)
if (nrow(W) != n) stop("nrow(Z) must equal nrow(W).")
if (is.null(colnames(W))) colnames(W) <- paste0("W", seq_len(ncol(W)))
base <- get_pairwise_interactions(
W, Z = if (length(ungrouped)) Z[, ungrouped, drop = FALSE] else NULL,
include_x_squared = include_x_squared
)
}

blocks <- list()
add_block <- function(M) {
keep <- apply(M, 2L, function(v) {
s <- stats::sd(v)
is.finite(s) && s > 1e-8
})
if (any(keep)) blocks[[length(blocks) + 1L]] <<- M[, keep, drop = FALSE]
invisible(NULL)
}

if (!is.null(W)) {
for (f in factors) {
lv <- grouped[z_groups[grouped] == f]
lv <- lv[n_support[lv] >= min_group_int_obs]
if (!length(lv)) next
for (j in seq_len(ncol(W))) {
M <- Z[, lv, drop = FALSE] * W[, j]
colnames(M) <- paste0(nmZ[lv], "*", colnames(W)[j])
add_block(M)
}
}
}

if (length(factors) > 1L) {
for (a in seq_len(length(factors) - 1L)) {
for (b in (a + 1L):length(factors)) {
lva <- grouped[z_groups[grouped] == factors[a]]
lvb <- grouped[z_groups[grouped] == factors[b]]
cols <- list()
for (i in lva) {
for (k in lvb) {
if (sum(z_support[, i] & z_support[, k]) < min_group_int_obs) next
v <- Z[, i] * Z[, k]
cols[[paste0(nmZ[i], "*", nmZ[k])]] <- v
}
}
if (!length(cols)) next
M <- do.call(cbind, cols)
colnames(M) <- names(cols)
add_block(M)
}
}
}

if (is.null(base) && !length(blocks)) return(NULL)
do.call(cbind, c(if (is.null(base)) list() else list(base), blocks))
}

# Name the factor and level on each side of a selected interaction column,
# e.g. HA_a3*Main_CS1 -> HapA:HA_a3 x Main_CS1.
annotate_z_group_interactions <- function(IntIndex, z_names, z_groups) {
if (is.null(z_groups) || is.null(IntIndex) || !nrow(IntIndex)) return(IntIndex)
parts <- strsplit(as.character(IntIndex$Variable), "*", fixed = TRUE)
side <- function(k) vapply(parts, function(x) if (length(x) >= k) x[k] else NA_character_, character(1))
factor_of <- function(term) unname(z_groups[match(term, z_names)])
t1 <- side(1L)
t2 <- side(2L)
f1 <- factor_of(t1)
f2 <- factor_of(t2)
label <- function(f, t) ifelse(is.na(f), t, paste0(f, ":", t))
IntIndex$Factor1 <- f1
IntIndex$Level1 <- t1
IntIndex$Factor2 <- f2
IntIndex$Level2 <- t2
IntIndex$Pair <- paste(label(f1, t1), label(f2, t2), sep = " x ")
IntIndex
}
