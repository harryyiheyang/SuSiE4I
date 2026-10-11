###############################################################################
# Factor groups in Z (groupint_ind).
#
# groupint_ind lists groups of Z columns, e.g. the indicator columns of one
# factor (baseline level dropped). Each Z column still interacts with the main
# CSs (following noint_env); the interaction stage treats the columns that pair
# one group with the same main CS as ONE group single effect (gsusie_ss).
# A group may be list(index = <columns>, baseline = <n-vector, 0/1>): its
# interaction columns are then built from orthonormal sum-to-zero contrasts of
# the full factor (index levels + baseline), so the group prior treats every
# level alike and no level acts as the reference.
###############################################################################

normalize_groupint_ind <- function(groupint_ind, Z) {
if (is.null(groupint_ind)) return(NULL)
if (!is.list(groupint_ind)) stop("groupint_ind must be a list of groups of Z columns.")
q <- ncol(Z)
nm <- colnames(Z)
labels <- names(groupint_ind)
if (is.null(labels)) labels <- rep("", length(groupint_ind))
labels[!nzchar(labels)] <- paste0("G", which(!nzchar(labels)))
if (anyDuplicated(labels)) stop("groupint_ind names must be unique.")
out <- rep(NA_character_, q)
baseline <- list()
for (k in seq_along(groupint_ind)) {
cols <- groupint_ind[[k]]
if (is.list(cols)) {
if (!is.null(cols$baseline)) {
if (!is.numeric(cols$baseline) || length(cols$baseline) != nrow(Z) || anyNA(cols$baseline)) {
stop("groupint_ind element '", labels[k], "': baseline must be a numeric 0/1 vector of length nrow(Z).")
}
baseline[[labels[k]]] <- as.numeric(cols$baseline)
}
cols <- cols$index
}
pos <- if (is.character(cols)) match(cols, nm) else as.integer(cols)
if (!length(pos) || any(is.na(pos)) || any(pos < 1L | pos > q)) {
stop("groupint_ind element '", labels[k], "' names columns not in Z.")
}
if (anyDuplicated(pos) || any(!is.na(out[pos]))) {
stop("A Z column appears in more than one groupint_ind group.")
}
out[pos] <- labels[k]
}
if (length(baseline)) attr(out, "baseline") <- baseline
out
}

# Z for interaction construction: for each group with a baseline, the index
# columns (dichotomised at their mean, so scaling does not matter) and the
# baseline form a one-hot factor D (K levels); the group's K - 1 columns are
# replaced by D %*% C, C orthonormal with columns summing to zero, centred.
# Main effects keep the original Z.
groupint_contrast_Z <- function(Z, groupint_ind) {
bl <- attr(groupint_ind, "baseline")
if (is.null(bl) || is.null(Z)) return(Z)
for (g in names(bl)) {
pos <- which(groupint_ind == g)
D <- cbind(vapply(pos, function(j) as.numeric(Z[, j] > mean(Z[, j])), numeric(nrow(Z))), as.numeric(bl[[g]] > 0))
if (any(rowSums(D) != 1)) stop("groupint_ind group '", g, "': index columns and baseline do not form a one-hot factor.")
H <- stats::contr.helmert(ncol(D))
Zc <- D %*% sweep(H, 2, sqrt(colSums(H^2)), "/")
Z[, pos] <- sweep(Zc, 2, colMeans(Zc))
}
Z
}

# Name the group and column on each side of a selected interaction column,
# e.g. dr_1*Main_CS1 -> Drink:dr_1 x Main_CS1.
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

# Group id of each interaction column for the group single-effect fit. Columns
# that pair the same two sides (a groupint_ind group with a main CS, e.g.
# Main_CS1 x Drink) form one group; every other column is a singleton.
groupint_column_groups <- function(namW, z_names, groupint_ind) {
if (is.null(groupint_ind) || is.null(namW)) return(NULL)
key <- vapply(strsplit(namW, "*", fixed = TRUE), function(x) {
g <- groupint_ind[match(x, z_names)]
paste(ifelse(is.na(g), x, g), collapse = "*")
}, character(1))
match(key, unique(key))
}
