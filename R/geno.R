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
#' allele counts as in BEDMatrix. Multiallelic PGEN variants are counted as REF
#' versus all non-REF alleles. To choose the counted direction, set the target
#' allele as REF upstream (`plink2 --ref-allele force target.txt 1 2
#' --make-pgen`); the counted allele is then everything else.
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
#' @return An object of class `geno` with fields `n`, `p`, `snp`,
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
    a1 <- bim[[5L]]
    a2 <- bim[[6L]]
    fid <- fam[[1L]]
    iid <- fam[[2L]]
  } else {
    files <- .geno_pgen_files(pgenfile)
    pv <- data.table::fread(files$pvar, skip = "#CHROM", colClasses = "character")
    snp <- pv[["ID"]]
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
                 snp = snp[vidx],
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

# Adds A1/A2 and Effect/Effect_SE to the discovery tables of a
# SuSiE4I() fit on a geno object. CS columns are oriented like the lead
# (highest-PIP) member, so the refit coefficient of a CS is per lead A1 copy up
# to the unit change below, which is exact for a single-variant CS and uses the
# lead's sd otherwise. Effect/Effect_SE are filled on the lead row of each CS:
# coefficient / sd(lead) (main), / (sd(lead1) sd(lead2)) (G x G) or / sd(lead)
# (E x G); sd = 1 when X is not scaled.
.add_alleles <- function(res, X) {
  sdx <- if (isTRUE(X$scale)) .geno_div(X) else rep(1, X$p)
  f <- res$fitJoint
  tab <- if (inherits(f, "coxph")) cox_coef_table(f)
         else if (inherits(f, "gam")) summary(f)$p.table
         else ocat_coef_table(f)
  eff <- function(cs, div) {
    r <- match(cs, rownames(tab))
    if (is.na(r)) return(c(NA_real_, NA_real_))
    c(tab[r, 1L], tab[r, 2L]) / div
  }
  main <- res$main_discoveries
  lead <- NULL
  if (is.data.frame(main) && nrow(main) > 0L && "Index" %in% names(main)) {
    i <- main$Index
    main$A1 <- X$a1[i]; main$A2 <- X$a2[i]
    is_lead <- !duplicated(main$CS)
    main$Effect <- NA_real_; main$Effect_SE <- NA_real_
    for (r in which(is_lead)) {
      e <- eff(main$CS[r], sdx[i[r]])
      main$Effect[r] <- e[1L]; main$Effect_SE[r] <- e[2L]
    }
    res$main_discoveries <- main
    lead <- main[is_lead, , drop = FALSE]
  }
  int <- res$interaction_discoveries
  if (is.data.frame(int) && nrow(int) > 0L && "Variable" %in% names(int) && !is.null(lead)) {
    parts <- strsplit(as.character(int$Variable), "*", fixed = TRUE)
    rr <- lapply(parts, function(v) {
      v <- v[grepl("^Main_CS[0-9]+$", v)]
      match(v, lead$CS)[seq_len(min(2L, length(v)))]
    })
    int$Effect <- NA_real_; int$Effect_SE <- NA_real_
    has_group <- if ("Group1" %in% names(int)) !is.na(int$Group1) | !is.na(int$Group2) else FALSE
    for (r in which(!duplicated(int$CS) & !has_group)) {
      d <- prod(sdx[lead$Index[rr[[r]]]])
      e <- eff(int$CS[r], d)
      int$Effect[r] <- e[1L]; int$Effect_SE[r] <- e[2L]
    }
    res$interaction_discoveries <- int
  }
  res
}
