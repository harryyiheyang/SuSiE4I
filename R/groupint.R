###############################################################################
# Factor groups in Z (groupint_ind).
#
# groupint_ind lists groups of Z columns, e.g. the indicator columns of one
# factor (baseline level dropped). Each Z column still interacts with the main
# CSs (following noint_env); the interaction stage treats the columns that pair
# one group with the same main CS as ONE group single effect (gsusie_ss).
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
