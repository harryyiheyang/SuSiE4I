// geno: a PLINK BED / PGEN genotype matrix held in memory as 2-bit planes
// (n p / 4 bytes) with the products SuSiE-style fine-mapping needs: X'X by
// exact popcount, X'M and X B by plane scans, and dense sub-blocks on demand.
// Values are A1 (BED) / ALT (PGEN) dosages as in BEDMatrix, missing calls filled with the variant's
// lower median, or with its observed mean (impute = "mean"), in which case
// X = X0 + M diag(mu0) with M the missing plane, kept as a second Planes.
// Centering and scaling are applied on the R side from mu and sd.

#include <Rcpp.h>
#include <fstream>
#include <string>
#include "geno_kernel.h"
#include "pgen_file.h"

namespace {

constexpr std::size_t kGroup = 64;  // variants per read group

struct Geno {
  std::size_t n = 0, p = 0;
  bool mean = false;
  Planes G;                     // lo = [x == 1], hi = [x == 2]
  Planes M;                     // mean imputation only: lo = missing, hi = 0
  std::vector<double> mu, sd;   // column mean and sd (n - 1) of the imputed matrix
};

inline int omp_thread() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline int resolve_threads(int threads) {
#ifdef _OPENMP
  return threads <= 0 ? omp_get_max_threads() : threads;
#else
  return 1;
#endif
}

// Packed column of the selected samples, in their requested order.
inline void gather_packed(const unsigned char* column, const std::vector<int>& rows,
                          unsigned char* buf) {
  std::fill(buf, buf + (rows.size() + 3) / 4, 0);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const unsigned int r = static_cast<unsigned int>(rows[i]);
    const unsigned int code = (column[r >> 2] >> (2 * (r & 3u))) & 3u;
    buf[i >> 2] |= static_cast<unsigned char>(code << (2 * (i & 3u)));
  }
}

inline void encode_column(Geno& g, const NibbleTable& table, const unsigned char* buf,
                          std::size_t j) {
  const std::size_t W = g.G.words;
  u64* lo = g.G.bits.data() + 2 * j * W;
  u64* miss = g.mean ? g.M.bits.data() + 2 * j * W : nullptr;
  encode_packed(table, buf, g.n, (g.n + 3) / 4, W, lo, lo + W, g.G.s[j], g.G.v[j], miss);
}

void init_geno(Geno& g, std::size_t n, std::size_t p, bool mean) {
  if (n < 2) Rcpp::stop("The genotype matrix needs at least two samples.");
  if (!p) Rcpp::stop("No variants selected.");
  g.n = n;
  g.p = p;
  g.mean = mean;
  reset_planes(g.G, n, p);
  if (mean) reset_planes(g.M, n, p);
}

void finish_stats(Geno& g) {
  const double n = static_cast<double>(g.n);
  g.mu.assign(g.p, 0.0);
  g.sd.assign(g.p, 0.0);
  for (std::size_t j = 0; j < g.p; ++j) {
    const double s = static_cast<double>(g.G.s[j]);
    const double v = static_cast<double>(g.G.v[j]);
    double mean, sumsq;
    if (g.mean) {
      std::int64_t cm = 0;
      const u64* miss = g.M.lo(j);
      for (std::size_t w = 0; w < g.M.words; ++w) cm += popcount64(miss[w]);
      const double obs = n - static_cast<double>(cm);
      mean = obs > 0 ? s / obs : 0.0;
      sumsq = v + static_cast<double>(cm) * mean * mean;
      g.M.s[j] = cm;
    } else {
      mean = s / n;
      sumsq = (v + s * s) / n;
    }
    g.mu[j] = mean;
    const double var = (sumsq - n * mean * mean) / (n - 1);
    g.sd[j] = var > 0 ? std::sqrt(var) : 0.0;
  }
}

std::vector<int> to_rows(const Rcpp::IntegerVector& x) {
  return std::vector<int>(x.begin(), x.end());
}

} // namespace

// [[Rcpp::export]]
SEXP geno_open_bed_cpp(const std::string& path, int n_total, int m_total,
                       const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                       bool mean, int threads) {
  threads = resolve_threads(threads);
  const std::size_t stride = (static_cast<std::size_t>(n_total) + 3) / 4;
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) Rcpp::stop("Cannot open BED file: " + path);
  const std::streamoff actual = input.tellg();
  if (actual < 0 || static_cast<std::uint64_t>(actual) != 3 + static_cast<std::uint64_t>(m_total) * stride)
    Rcpp::stop("BED size does not match its BIM and FAM files: " + path);
  input.seekg(0);
  unsigned char header[3];
  input.read(reinterpret_cast<char*>(header), 3);
  if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("BED must be a PLINK SNP-major BED file: " + path);

  Rcpp::XPtr<Geno> g(new Geno, true);
  const std::vector<int> rows = to_rows(samples);
  init_geno(*g, rows.size(), variants.size(), mean);
  static const NibbleTable table({2, -1, 1, 0});  // A1 dosage: 00:2, 01:NA, 10:1, 11:0
  const std::size_t sel_stride = (rows.size() + 3) / 4;
  std::vector<unsigned char> buffer(kGroup * stride);
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  const std::size_t p = g->p;
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    for (std::size_t k = 0; k < count; ++k) {
      input.seekg(static_cast<std::streamoff>(3 + static_cast<std::uint64_t>(variants[g0 + k]) * stride));
      input.read(reinterpret_cast<char*>(buffer.data() + k * stride), static_cast<std::streamsize>(stride));
      if (!input) Rcpp::stop("Could not read BED file: " + path);
    }
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      unsigned char* buf = packed.data() + static_cast<std::size_t>(omp_thread()) * sel_stride;
      gather_packed(buffer.data() + k * stride, rows, buf);
      encode_column(*g, table, buf, g0 + k);
    }
    Rcpp::checkUserInterrupt();
  }
  finish_stats(*g);
  return g;
}

// [[Rcpp::export]]
SEXP geno_open_pgen_cpp(const std::string& path, int n_total, const Rcpp::IntegerVector& allele_ct,
                        const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                        bool mean, int threads) {
  threads = resolve_threads(threads);
  PgenFile file(path, static_cast<std::size_t>(n_total), allele_ct, threads);
  Rcpp::XPtr<Geno> g(new Geno, true);
  const std::vector<int> rows = to_rows(samples);
  init_geno(*g, rows.size(), variants.size(), mean);
  const NibbleTable& table = PgenFile::table();
  const std::size_t sel_stride = (rows.size() + 3) / 4;
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  const std::size_t p = g->p;
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    int failed = 0;
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = omp_thread();
      const unsigned char* column = file.read(static_cast<std::size_t>(variants[g0 + k]), t);
      if (!column) {
        #pragma omp atomic write
        failed = 1;
        continue;
      }
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      gather_packed(column, rows, buf);
      encode_column(*g, table, buf, g0 + k);
    }
    if (failed) Rcpp::stop("Could not read PGEN file: " + path);
    Rcpp::checkUserInterrupt();
  }
  finish_stats(*g);
  return g;
}

// [[Rcpp::export]]
Rcpp::List geno_info_cpp(SEXP ptr) {
  Rcpp::XPtr<Geno> g(ptr);
  return Rcpp::List::create(Rcpp::Named("n") = static_cast<double>(g->n),
                            Rcpp::Named("p") = static_cast<double>(g->p),
                            Rcpp::Named("mean") = g->mu,
                            Rcpp::Named("sd") = g->sd);
}

// X'X of the imputed, unstandardized matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtx_cpp(SEXP ptr, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const TileKernel kernel = select_kernel();
  const std::size_t p = g->p;
  const double n = static_cast<double>(g->n);
  Rcpp::NumericMatrix out(p, p);
  cor_block(g->G, 0, p, g->G, 0, p, true, n, out.begin(), p, threads, kernel, true);
  if (g->mean) {
    const std::vector<double>& mu = g->mu;
    Rcpp::NumericMatrix t(p, p);
    cor_block(g->G, 0, p, g->M, 0, p, false, n, t.begin(), p, threads, kernel, true);
    for (std::size_t b = 0; b < p; ++b)
      for (std::size_t a = 0; a < p; ++a)
        out[a + b * p] += t[a + b * p] * mu[b] + t[b + a * p] * mu[a];
    std::fill(t.begin(), t.end(), 0.0);
    cor_block(g->M, 0, p, g->M, 0, p, true, n, t.begin(), p, threads, kernel, true);
    for (std::size_t b = 0; b < p; ++b)
      for (std::size_t a = 0; a < p; ++a)
        out[a + b * p] += mu[a] * mu[b] * t[a + b * p];
  }
  return out;
}

// X'M for a dense n x k matrix M.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtm_cpp(SEXP ptr, const Rcpp::NumericMatrix& M, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p, W = g->G.words;
  const std::size_t k = M.ncol();
  if (static_cast<std::size_t>(M.nrow()) != n) Rcpp::stop("M must have one row per sample.");
  Rcpp::NumericMatrix out(p, k);
  const double* m = M.begin();
  double* o = out.begin();
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    std::vector<double> s1(k, 0.0), s2(k, 0.0), sm(k, 0.0);
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    const u64* miss = g->mean ? g->M.lo(j) : nullptr;
    for (std::size_t w = 0; w < W; ++w) {
      for (u64 bits = lo[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) s1[c] += m[i + c * n];
      }
      for (u64 bits = hi[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) s2[c] += m[i + c * n];
      }
      if (miss)
        for (u64 bits = miss[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) sm[c] += m[i + c * n];
        }
    }
    for (std::size_t c = 0; c < k; ++c)
      o[j + c * p] = s1[c] + 2.0 * s2[c] + (miss ? g->mu[j] * sm[c] : 0.0);
  }
  return out;
}

// X B for a dense p x k matrix B; all-zero rows of B are skipped.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_mx_cpp(SEXP ptr, const Rcpp::NumericMatrix& B, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p, W = g->G.words;
  const std::size_t k = B.ncol();
  if (static_cast<std::size_t>(B.nrow()) != p) Rcpp::stop("B must have one row per variant.");
  Rcpp::NumericMatrix out(n, k);
  const double* b = B.begin();
  double* o = out.begin();
  std::vector<char> nonzero(p, 0);
  for (std::size_t j = 0; j < p; ++j)
    for (std::size_t c = 0; c < k; ++c)
      if (b[j + c * p] != 0.0) { nonzero[j] = 1; break; }
  const std::size_t used_words = (n + 63) / 64;
  #pragma omp parallel num_threads(threads)
  {
    const int nt = resolve_threads(threads), t = omp_thread();
    const std::size_t w0 = used_words * t / nt, w1 = used_words * (t + 1) / nt;
    std::vector<double> bj(k);
    for (std::size_t j = 0; j < p; ++j) {
      if (!nonzero[j]) continue;
      for (std::size_t c = 0; c < k; ++c) bj[c] = b[j + c * p];
      const u64* lo = g->G.lo(j);
      const u64* hi = g->G.hi(j);
      const u64* miss = g->mean ? g->M.lo(j) : nullptr;
      for (std::size_t w = w0; w < w1; ++w) {
        for (u64 bits = lo[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) o[i + c * n] += bj[c];
        }
        for (u64 bits = hi[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) o[i + c * n] += 2.0 * bj[c];
        }
        if (miss) {
          const double mj = g->mu[j];
          for (u64 bits = miss[w]; bits; bits &= bits - 1) {
            const std::size_t i = 64 * w + __builtin_ctzll(bits);
            for (std::size_t c = 0; c < k; ++c) o[i + c * n] += mj * bj[c];
          }
        }
      }
    }
  }
  return out;
}

// Dense block X[rows, cols] of the imputed, unstandardized matrix (0-based).
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_dense_cpp(SEXP ptr, const Rcpp::IntegerVector& rows,
                                   const Rcpp::IntegerVector& cols, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t nr = rows.size(), nc = cols.size();
  Rcpp::NumericMatrix out(nr, nc);
  double* o = out.begin();
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(nc); ++cc) {
    const std::size_t j = static_cast<std::size_t>(cols[cc]);
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    const u64* miss = g->mean ? g->M.lo(j) : nullptr;
    const double mj = g->mu[j];
    double* col = o + static_cast<std::size_t>(cc) * nr;
    for (std::size_t r = 0; r < nr; ++r) {
      const std::size_t i = static_cast<std::size_t>(rows[r]);
      const std::size_t w = i / 64;
      const u64 bit = u64(1) << (i % 64);
      col[r] = (lo[w] & bit) ? 1.0 : (hi[w] & bit) ? 2.0 : (miss && (miss[w] & bit)) ? mj : 0.0;
    }
  }
  return out;
}
