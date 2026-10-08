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
#' Keeps the genotypes of a PLINK BED or PGEN file in 2-bit form, so SuSiE4I
#' never holds an n by p double matrix in R. Values are A1 (BED) or ALT (PGEN)
#' allele counts as in BEDMatrix; multiallelic PGEN variants count all non-REF
#' alleles.
#'
#' @param bedfile PLINK 1 BED path or prefix (`.bim` and `.fam` alongside).
#' @param pgenfile PLINK 2 PGEN path or prefix (hardcall-only, uncompressed
#'   `.pvar` and `.psam` alongside). Give either `bedfile` or `pgenfile`.
#' @param snp_vec Variants to keep, as IDs or 1-based file indices, in the
#'   column order wanted. `NULL` keeps all in file order. Selection by ID takes
#'   the first match, so use indices when IDs are duplicated (e.g. ".").
#' @param sample_vec Samples to keep, as IIDs (or `FID_IID`) or 1-based file
#'   indices, in the row order wanted. `NULL` keeps all in file order. The rows
#'   of `y`, `Z` and `status` passed to `SuSiE4I()` must be in this order (or in
#'   `.fam`/`.psam` order when `NULL`). IIDs repeated across families need the
#'   `FID_IID` form.
#' @param impute `"median"` (default, the variant's lower median) or `"mean"`
#'   (the variant's observed mean; stores one more bit plane).
#' @param scale Whether products and dense blocks are of the standardized
#'   matrix.
#' @param threads Number of OpenMP threads.
#' @return An object of class `geno` with fields `n`, `p`, `snp`, `chr`, `pos`,
#'   `a1`, `a2` (aligned with `snp`; `a1` is the counted allele: bim column 5
#'   for BED, ALT for PGEN, and `a2` the other allele: bim column 6 or REF),
#'   `sample` (IIDs in row order), `impute`, `scale`, `center` and `sd`. It holds an
#'   external pointer, so it is valid only in the current R session (it cannot
#'   be saved with `saveRDS` or sent to parallel workers).
#' @export
geno_open <- function(bedfile = NULL, pgenfile = NULL, snp_vec = NULL,
                      sample_vec = NULL, impute = c("median", "mean"),
                      scale = TRUE, threads = 4L) {
  impute <- match.arg(impute)
  threads <- as.integer(threads)
  if (is.null(bedfile) == is.null(pgenfile)) {
    stop("Give exactly one of bedfile and pgenfile.", call. = FALSE)
  }
  if (!is.null(bedfile)) {
    bed <- .geno_bed_path(bedfile)
    base <- sub("\\.bed$", "", bed)
    bim <- data.table::fread(paste0(base, ".bim"), header = FALSE, colClasses = "character")
    fam <- data.table::fread(paste0(base, ".fam"), header = FALSE, colClasses = "character")
    snp <- bim[[2L]]
    chr <- bim[[1L]]
    pos <- bim[[4L]]
    a1 <- bim[[5L]]
    a2 <- bim[[6L]]
    fid <- fam[[1L]]
    iid <- fam[[2L]]
  } else {
    files <- .geno_pgen_files(pgenfile)
    pv <- data.table::fread(files$pvar, skip = "#CHROM", colClasses = "character")
    snp <- pv[["ID"]]
    chr <- pv[["#CHROM"]]
    pos <- pv[["POS"]]
    a1 <- pv[["ALT"]]
    a2 <- pv[["REF"]]
    allele_ct <- 1L + nchar(pv$ALT) - nchar(gsub(",", "", pv$ALT, fixed = TRUE))
    ps <- data.table::fread(files$psam, colClasses = "character")
    iid <- ps[[if ("#IID" %in% names(ps)) "#IID" else "IID"]]
    fid_col <- intersect(c("#FID", "FID"), names(ps))
    fid <- if (length(fid_col)) ps[[fid_col[1L]]] else rep("0", length(iid))
  }
  vidx <- .geno_index(snp_vec, snp, NULL, "snp_vec")
  sidx <- .geno_index(sample_vec, iid, paste(fid, iid, sep = "_"), "sample_vec")
  ptr <- if (!is.null(bedfile)) {
    geno_open_bed_cpp(bed, length(iid), length(snp), vidx - 1L, sidx - 1L,
                      impute == "mean", threads)
  } else {
    geno_open_pgen_cpp(files$pgen, length(iid), allele_ct,
                       vidx - 1L, sidx - 1L, impute == "mean", threads)
  }
  info <- geno_info_cpp(ptr)
  structure(list(ptr = ptr, n = length(sidx), p = length(vidx),
                 snp = snp[vidx], chr = chr[vidx], pos = pos[vidx],
                 a1 = a1[vidx], a2 = a2[vidx], sample = iid[sidx],
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

# crossprod(X) by exact popcount, or crossprod(X, M) for a dense M with one row per sample.
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

# X %*% B for B with one row per variant; a vector B gives a vector.
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

# crossprod(X, X * w) and crossprod(X, M) from dense row blocks of at most block_size rows.
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

.geno_bed_path <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("bedfile must be one BED file path.", call. = FALSE)
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
  for (path in c(pgen, pvar, psam)) {
    if (!file.exists(path)) stop("Required PGEN companion file is missing: ",
                                 path, call. = FALSE)
  }
  list(pgen = normalizePath(pgen, mustWork = TRUE),
       pvar = normalizePath(pvar, mustWork = TRUE),
       psam = normalizePath(psam, mustWork = TRUE))
}

# Adds CHR/POS/A1/A2 of the reported variants to the discovery tables of a
# SuSiE4I() fit on a geno object. Main rows are located by column index and
# get Sign (+1/-1: direction of the variant's effect on A1; the refit column of
# a credible set is sign-aligned, so its coefficient times Sign is the A1
# effect). Interaction terms are built from credible-set columns ("Main_CSk"),
# so each such factor is reported through the credible set's lead variant
# (highest PIP) and Sign_k; a factor that is a Z column has no variant (NA).
.add_alleles <- function(res, X) {
  fit_sign <- function(l, j) {
    m <- res$fitX
    m <- if (is.null(m)) NULL else if (is.null(m$group)) m$mu else m$unit_mu
    if (is.null(m) || is.na(l) || j > ncol(m)) return(NA_real_)
    sign(m[l, j])
  }
  main <- res$main_discoveries
  if (is.data.frame(main) && nrow(main) > 0L && "Index" %in% names(main)) {
    i <- main$Index
    l <- suppressWarnings(as.integer(sub("^Main_CS", "", main$CS)))
    main$CHR <- X$chr[i]; main$POS <- X$pos[i]
    main$A1 <- X$a1[i]; main$A2 <- X$a2[i]
    main$Sign <- mapply(fit_sign, l, i)
    res$main_discoveries <- main
    lead <- main[!duplicated(main$CS), , drop = FALSE]
  }
  int <- res$interaction_discoveries
  if (is.data.frame(int) && nrow(int) > 0L && "Variable" %in% names(int)) {
    parts <- strsplit(as.character(int$Variable), "*", fixed = TRUE)
    pick <- function(v, k) {
      v <- v[grepl("^Main_CS[0-9]+$", v)]
      if (length(v) < k) return(NA_integer_)
      match(v[k], lead$CS)
    }
    for (k in 1:2) {
      r <- vapply(parts, pick, 1L, k = k)
      i <- lead$Index[r]
      int[[paste0("CHR_", k)]] <- X$chr[i]; int[[paste0("POS_", k)]] <- X$pos[i]
      int[[paste0("A1_", k)]] <- X$a1[i]; int[[paste0("A2_", k)]] <- X$a2[i]
      int[[paste0("Sign_", k)]] <- lead$Sign[r]
    }
    res$interaction_discoveries <- int
  }
  res
}
