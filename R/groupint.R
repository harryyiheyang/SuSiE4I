###############################################################################
# Pairwise interactions between groups of Z columns (e.g. haplotype dummies).
#
# groupint_ind lists groups of Z columns. For every pair of different groups,
# each column of one group times each column of the other is added to the
# interaction design as its own candidate; columns within a group are never
# paired. Z x main-CS interactions are unchanged and still follow noint_env.
###############################################################################

normalize_groupint_ind <- function(groupint_ind, Z) {
if (is.null(groupint_ind)) return(NULL)
if (!is.list(groupint_ind) || length(groupint_ind) < 2L) {
stop("groupint_ind must be a list of at least two groups of Z columns.")
}
q <- ncol(Z)
nm <- colnames(Z)
labels <- names(groupint_ind)
if (is.null(labels)) labels <- rep("", length(groupint_ind))
labels[!nzchar(labels)] <- paste0("G", which(!nzchar(labels)))
if (anyDuplicated(labels)) stop("groupint_ind names must be unique.")
out <- rep(NA_character_, q)
for (k in seq_along(groupint_ind)) {
cols <- groupint_ind[[k]]
pos <- if (is.character(cols)) match(cols, nm) else as.integer(cols)
if (!length(pos) || any(is.na(pos)) || any(pos < 1L | pos > q)) {
stop("groupint_ind element '", labels[k], "' names columns not in Z.")
}
if (anyDuplicated(pos) || any(!is.na(out[pos]))) {
stop("A Z column appears in more than one groupint_ind group.")
}
out[pos] <- labels[k]
}
out
}

get_groupint_interactions <- function(Z, groupint_ind, min_xtx = 1e-8) {
Z <- as.matrix(Z)
nmZ <- colnames(Z)
if (is.null(nmZ)) nmZ <- paste0("Z", seq_len(ncol(Z)))
groups <- unique(stats::na.omit(groupint_ind))
cols <- list()
for (a in seq_len(length(groups) - 1L)) {
for (b in (a + 1L):length(groups)) {
for (i in which(groupint_ind == groups[a])) {
for (j in which(groupint_ind == groups[b])) {
v <- Z[, i] * Z[, j]
if (sum(v^2) / length(v) < min_xtx) next
cols[[paste0(nmZ[i], "*", nmZ[j])]] <- v
}
}
}
}
if (!length(cols)) return(NULL)
out <- do.call(cbind, cols)
colnames(out) <- names(cols)
out
}

# Name the group and column on each side of a selected interaction column,
# e.g. HA_a1*HB_b2 -> HapA:HA_a1 x HapB:HB_b2.
annotate_groupint_interactions <- function(IntIndex, z_names, groupint_ind) {
if (is.null(groupint_ind) || is.null(IntIndex) || !nrow(IntIndex)) return(IntIndex)
parts <- strsplit(as.character(IntIndex$Variable), "*", fixed = TRUE)
side <- function(k) vapply(parts, function(x) if (length(x) >= k) x[k] else NA_character_, character(1))
group_of <- function(term) unname(groupint_ind[match(term, z_names)])
t1 <- side(1L)
t2 <- side(2L)
g1 <- group_of(t1)
g2 <- group_of(t2)
label <- function(g, t) ifelse(is.na(g), t, paste0(g, ":", t))
IntIndex$Group1 <- g1
IntIndex$Term1 <- t1
IntIndex$Group2 <- g2
IntIndex$Term2 <- t2
IntIndex$Pair <- paste(label(g1, t1), label(g2, t2), sep = " x ")
IntIndex
}

check_coverage_nonkilled <- function(x) {
if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x > 1) {
stop("coverage_nonkilled must be a number in (0, 1].")
}
invisible(TRUE)
}
