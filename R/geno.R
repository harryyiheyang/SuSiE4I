# geno_open(): a PLINK BED or PGEN genotype matrix kept in memory as 2-bit
# planes (n p / 4 bytes) instead of an n x p double matrix. The object behaves
# like a numeric matrix for dim(), dimnames() and [ , and gives the products
# fine-mapping needs through geno_crossprod(), geno_multiply() and
# geno_wcrossprod(). Values are A1 (BED) or ALT (PGEN) allele dosages, as in
# BEDMatrix; missing calls are filled with
# the variant's lower median (impute = "median") or its observed mean
# ("mean"); with scale = TRUE every product is that of the column-standardized
# matrix (mean 0, sd 1 with n - 1; constant columns are only centered).

#' Open a BED or PGEN file as a compact genotype matrix
#'
#' @param bedfile PLINK 1 BED path or prefix (`.bim` and `.fam` alongside).
#' @param pgenfile PLINK 2 PGEN path or prefix (hardcall-only, uncompressed
#'   `.pvar` and `.psam` alongside). Give either `bedfile` or `pgenfile`.
#' @param snp_vec Variants to keep, as IDs or 1-based file indices, in the
#'   column order wanted. `NULL` keeps all in file order.
#' @param sample_vec Samples to keep, as IIDs (or `FID_IID`) or 1-based file
#'   indices, in the row order wanted. `NULL` keeps all in file order.
#' @param impute `"median"` (default, the variant's lower median) or `"mean"`
#'   (the variant's observed mean; stores one more bit plane).
#' @param scale Whether products and dense blocks are of the standardized
#'   matrix.
#' @param threads Number of OpenMP threads.
#' @return An object of class `geno`.
#' @export
geno_open <- function(bedfile = NULL, pgenfile = NULL, snp_vec = NULL,
                      sample_vec = NULL, impute = c("median", "mean"),
                      scale = TRUE, threads = 4L) {
  impute <- match.arg(impute)
  threads <- .geno_threads(threads)
  if (is.null(bedfile) == is.null(pgenfile)) {
    stop("Give exactly one of bedfile and pgenfile.", call. = FALSE)
  }
  if (!is.null(bedfile)) {
    bed <- .geno_bed_path(bedfile, "bedfile")
    base <- sub("\\.bed$", "", bed)
    bim <- utils::read.table(paste0(base, ".bim"), colClasses = "character",
                             comment.char = "", quote = "")
    fam <- utils::read.table(paste0(base, ".fam"), colClasses = "character",
                             comment.char = "", quote = "")
    snp <- bim[[2L]]
    fid <- fam[[1L]]
    iid <- fam[[2L]]
  } else {
    files <- .geno_pgen_files(pgenfile)
    variants <- .geno_variants(files$pvar)
    snp <- variants$id
    ids <- strsplit(.geno_samples(files$psam), "\t", fixed = TRUE)
    fid <- vapply(ids, `[`, "", 1L)
    iid <- vapply(ids, `[`, "", 2L)
  }
  vidx <- .geno_index(snp_vec, snp, NULL, "snp_vec")
  sidx <- .geno_index(sample_vec, iid, paste(fid, iid, sep = "_"), "sample_vec")
  ptr <- if (!is.null(bedfile)) {
    geno_open_bed_cpp(bed, length(iid), length(snp), vidx - 1L, sidx - 1L,
                      impute == "mean", threads)
  } else {
    geno_open_pgen_cpp(files$pgen, length(iid), variants$allele_ct,
                       vidx - 1L, sidx - 1L, impute == "mean", threads)
  }
  info <- geno_info_cpp(ptr)
  structure(list(ptr = ptr, n = length(sidx), p = length(vidx),
                 snp = snp[vidx], sample = iid[sidx],
                 impute = impute, scale = scale,
                 center = info$mean, sd = info$sd, threads = threads),
            class = "geno")
}

.geno_index <- function(x, ids, alt_ids, arg) {
  if (is.null(x)) return(seq_along(ids))
  if (is.character(x)) {
    idx <- match(x, ids)
    if (!is.null(alt_ids) && anyNA(idx)) idx[is.na(idx)] <- match(x[is.na(idx)], alt_ids)
    if (anyNA(idx)) {
      stop(sum(is.na(idx)), " entries of ", arg, " are not in the file, e.g. ",
           x[which(is.na(idx))[1L]], call. = FALSE)
    }
    return(idx)
  }
  idx <- as.integer(x)
  if (anyNA(idx) || any(idx < 1L | idx > length(ids))) {
    stop(arg, " indices must lie between 1 and ", length(ids), ".", call. = FALSE)
  }
  idx
}

# Scale divisor: constant columns are centered but not divided (as scale()).
.geno_div <- function(X) ifelse(X$sd > 0, X$sd, 1)

#' @export
dim.geno <- function(x) c(x$n, x$p)

#' @export
dimnames.geno <- function(x) list(NULL, x$snp)

#' @export
print.geno <- function(x, ...) {
  cat(sprintf("geno matrix: %d samples x %d variants, impute = %s, scale = %s (%.2f GB in memory)\n",
              x$n, x$p, x$impute, x$scale,
              (if (x$impute == "mean") 2 else 1) * x$n * x$p / 4 / 1e9))
  invisible(x)
}

#' @export
`[.geno` <- function(x, i, j, drop = TRUE) {
  rows <- if (missing(i)) seq_len(x$n) else seq_len(x$n)[i]
  cols <- if (missing(j)) seq_len(x$p) else if (is.character(j)) match(j, x$snp) else seq_len(x$p)[j]
  if (anyNA(rows) || anyNA(cols)) stop("subscript out of bounds", call. = FALSE)
  G <- geno_dense_cpp(x$ptr, rows - 1L, cols - 1L, x$threads)
  if (x$scale) {
    G <- sweep(sweep(G, 2L, x$center[cols], "-"), 2L, .geno_div(x)[cols], "/")
  }
  colnames(G) <- x$snp[cols]
  if (drop) drop(G) else G
}

#' @export
as.matrix.geno <- function(x, ...) x[, , drop = FALSE]

#' Cross products of a geno matrix
#'
#' `geno_crossprod(X)` is `crossprod(X)` by exact popcount; `geno_crossprod(X, M)`
#' is `crossprod(X, M)` for a dense `M` with one row per sample.
#'
#' @param X A `geno` object from [geno_open()].
#' @param M Optional dense matrix (or vector) with one row per sample.
#' @param threads Number of OpenMP threads.
#' @return A p by p matrix when `M` is `NULL`, otherwise a p by `ncol(M)` matrix.
#' @export
geno_crossprod <- function(X, M = NULL, threads = X$threads) {
  if (is.null(M)) {
    out <- geno_xtx_cpp(X$ptr, threads)
    if (X$scale) {
      d <- .geno_div(X)
      out <- (out - X$n * tcrossprod(X$center)) / tcrossprod(d)
    }
    dimnames(out) <- list(X$snp, X$snp)
    return(out)
  }
  M <- as.matrix(M) + 0
  out <- geno_xtm_cpp(X$ptr, M, threads)
  if (X$scale) out <- (out - tcrossprod(X$center, colSums(M))) / .geno_div(X)
  rownames(out) <- X$snp
  out
}

#' Multiply a geno matrix by a dense matrix or vector
#'
#' `geno_multiply(X, B)` is `X %*% B` for `B` with one row per variant; a
#' vector `B` gives a vector.
#'
#' @param X A `geno` object from [geno_open()].
#' @param B Dense matrix with one row per variant, or a numeric vector.
#' @param threads Number of OpenMP threads.
#' @return An n by `ncol(B)` matrix, or a numeric vector when `B` is a vector.
#' @export
geno_multiply <- function(X, B, threads = X$threads) {
  vec <- is.null(dim(B))
  B <- as.matrix(B) + 0
  if (X$scale) {
    B <- B / .geno_div(X)
    out <- geno_mx_cpp(X$ptr, B, threads)
    out <- sweep(out, 2L, as.numeric(crossprod(X$center, B)), "-")
  } else {
    out <- geno_mx_cpp(X$ptr, B, threads)
  }
  if (vec) as.numeric(out) else out
}

#' Row-chunked weighted cross products of a geno matrix
#'
#' Returns `crossprod(X, X * w)` and `crossprod(X, M)` from dense row blocks
#' of at most `block_size` rows (capped so a block stays near 512 MB).
#'
#' @param X A `geno` object from [geno_open()].
#' @param w Numeric vector of n row weights.
#' @param M Optional dense matrix with one row per sample.
#' @param block_size Maximum number of rows per dense block.
#' @return A list with `XtWX` (p by p) and `XtM` (p by `ncol(M)`).
#' @export
geno_wcrossprod <- function(X, w, M = NULL, block_size = 10000L) {
  w <- as.numeric(w)
  M <- if (is.null(M)) matrix(0, X$n, 0L) else as.matrix(M) + 0
  rows <- max(1L, min(as.integer(block_size), 2^26 %/% X$p))
  XtWX <- matrix(0, X$p, X$p)
  XtM <- matrix(0, X$p, ncol(M))
  for (start in seq.int(1L, X$n, by = rows)) {
    idx <- start:min(X$n, start + rows - 1L)
    G <- X[idx, , drop = FALSE]
    XtWX <- XtWX + crossprod(G, G * w[idx])
    if (ncol(M)) XtM <- XtM + crossprod(G, M[idx, , drop = FALSE])
  }
  dimnames(XtWX) <- list(X$snp, X$snp)
  rownames(XtM) <- X$snp
  list(XtWX = XtWX, XtM = XtM)
}

.geno_bed_path <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(arg, " must be a single BED file path or prefix.", call. = FALSE)
  }
  if (!grepl("\\.bed$", x)) x <- paste0(x, ".bed")
  if (!file.exists(x)) stop("BED file not found: ", x, call. = FALSE)
  normalizePath(x, winslash = "/")
}

.geno_pgen_files <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("pgenfile must be a single PGEN file path or prefix.", call. = FALSE)
  }
  prefix <- sub("\\.pgen$", "", x, ignore.case = TRUE)
  pgen <- paste0(prefix, ".pgen")
  pvar <- paste0(prefix, ".pvar")
  psam <- paste0(prefix, ".psam")
  if (!file.exists(pvar) && file.exists(paste0(pvar, ".zst"))) {
    stop("Decompress the .pvar.zst first, e.g. plink2 --pfile ", prefix,
         " vzs --make-just-pvar --out ", prefix, call. = FALSE)
  }
  for (path in c(pgen, pvar, psam)) {
    if (!file.exists(path)) stop("Required PGEN companion file is missing: ",
                                 path, call. = FALSE)
  }
  list(pgen = normalizePath(pgen, mustWork = TRUE),
       pvar = normalizePath(pvar, mustWork = TRUE),
       psam = normalizePath(psam, mustWork = TRUE))
}

# First line that does not start with "##" and the number of lines up to and
# including it.
.geno_header <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return(list(line = NULL, skip = skip))
    skip <- skip + 1L
    if (!startsWith(line, "##")) return(list(line = line, skip = skip))
  }
}

# Only the named columns of a header-described text file, as character.
.geno_columns <- function(path, header, cols, keep, sep) {
  classes <- rep("NULL", length(cols))
  classes[cols %in% keep] <- "character"
  utils::read.table(path, sep = sep, skip = header$skip, header = FALSE,
                    colClasses = classes, col.names = cols,
                    comment.char = "", quote = "",
                    na.strings = character(0), check.names = FALSE)
}

.geno_samples <- function(path) {
  header <- .geno_header(path)
  if (is.null(header$line)) {
    stop("PSAM must contain a header and at least one sample: ", path,
         call. = FALSE)
  }
  cols <- strsplit(trimws(header$line), "[[:space:]]+")[[1L]]
  iid_col <- if ("#IID" %in% cols) "#IID" else "IID"
  if (!iid_col %in% cols) {
    stop("PSAM is missing its IID column: ", path, call. = FALSE)
  }
  fid_col <- if ("#FID" %in% cols) "#FID" else "FID"
  tab <- .geno_columns(path, header, cols, c(iid_col, fid_col), "")
  if (!nrow(tab)) {
    stop("PSAM must contain a header and at least one sample: ", path,
         call. = FALSE)
  }
  fid <- if (fid_col %in% cols) tab[[fid_col]] else rep("0", nrow(tab))
  paste(fid, tab[[iid_col]], sep = "\t")
}

# Variant IDs and allele counts (REF + ALTs) of a plain-text .pvar, in file
# order. "##" lines are skipped; the "#CHROM" header line names the columns.
.geno_variants <- function(path) {
  header <- .geno_header(path)
  if (is.null(header$line) || !startsWith(header$line, "#CHROM")) {
    stop("PVAR must have a #CHROM header line: ", path, call. = FALSE)
  }
  cols <- strsplit(sub("^#", "", header$line), "\t", fixed = TRUE)[[1L]]
  if (!all(c("ID", "ALT") %in% cols)) {
    stop("PVAR header is missing its ID or ALT column: ", path, call. = FALSE)
  }
  tab <- .geno_columns(path, header, cols, c("ID", "ALT"), "\t")
  if (!nrow(tab)) stop("PVAR file has no variants: ", path, call. = FALSE)
  alt <- tab[["ALT"]]
  n_alt <- nchar(alt) - nchar(gsub(",", "", alt, fixed = TRUE)) + 1L
  list(id = tab[["ID"]], allele_ct = as.integer(1L + n_alt))
}

.geno_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      !is.finite(threads) || threads < 1 || threads != round(threads) ||
      threads > .Machine$integer.max) {
    stop("threads must be a single positive integer.", call. = FALSE)
  }
  as.integer(threads)
}
