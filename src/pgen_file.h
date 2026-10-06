// PgenFile: hardcall-only PGEN reader on the bundled pgenlib, one PgenReader
// per thread, used by the geno object.
#ifndef SUSIE4I_PGEN_FILE_H
#define SUSIE4I_PGEN_FILE_H

#include <Rcpp.h>
#include <memory>
#include <string>
#include "geno_kernel.h"
#include "pgenlib/pgenlib_read.h"

namespace {

struct PgenHeader {
  plink2::PgenFileInfo info;
  unsigned char* alloc = nullptr;
  std::vector<uintptr_t> allele_idx_offsets, nonref_flags;
  PgenHeader() { plink2::PreinitPgfi(&info); }
  PgenHeader(const PgenHeader&) = delete;
  PgenHeader& operator=(const PgenHeader&) = delete;
  ~PgenHeader() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(&info, &e);
    if (alloc) plink2::aligned_free(alloc);
  }
};

struct ThreadReader {
  plink2::PgenReader pgr;
  plink2::PgrSampleSubsetIndex pssi;
  unsigned char* alloc = nullptr;
  unsigned char* genovec = nullptr;
  ThreadReader() { plink2::PreinitPgr(&pgr); }
  ThreadReader(const ThreadReader&) = delete;
  ThreadReader& operator=(const ThreadReader&) = delete;
  ~ThreadReader() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgr(&pgr, &e);
    if (alloc) plink2::aligned_free(alloc);
    if (genovec) plink2::aligned_free(genovec);
  }
};

class PgenFile {
public:
  PgenFile(const std::string& path, std::size_t n, const Rcpp::IntegerVector& allele_ct,
           int threads)
    : path_(path) {
    plink2::PgenFileInfo& info = header_.info;
    char err[plink2::kPglErrstrBufBlen];
    err[0] = '\0';
    plink2::PgenHeaderCtrl header_ctrl;
    uintptr_t pgfi_cachelines;
    if (plink2::PgfiInitPhase1(path.c_str(), nullptr, UINT32_MAX, UINT32_MAX, &header_ctrl,
                               &info, &pgfi_cachelines, err) != plink2::kPglRetSuccess)
      fail(err);
    if (info.gflags & plink2::kfPgenGlobalDosagePresent) fail_dosage();
    if (info.raw_sample_ct != n)
      Rcpp::stop("PGEN sample count does not match its PSAM file: " + path);
    const std::size_t m = info.raw_variant_ct;
    if (static_cast<std::size_t>(allele_ct.size()) != m)
      Rcpp::stop("PGEN variant count does not match its PVAR file: " + path);

    // Allele counts come from the .pvar, as in pgenlibr with a pvar; pgenlib
    // checks them against the .pgen when it stores them.
    uint32_t max_allele_ct = 2;
    for (int a : allele_ct) max_allele_ct = std::max<uint32_t>(max_allele_ct, a);
    info.max_allele_ct = 2;
    if (max_allele_ct > 2 || (header_ctrl & 0x30)) {
      std::vector<uintptr_t>& offsets = header_.allele_idx_offsets;
      offsets.resize(m + 1);
      offsets[0] = 0;
      for (std::size_t j = 0; j < m; ++j) offsets[j + 1] = offsets[j] + allele_ct[j];
      info.allele_idx_offsets = offsets.data();
      info.max_allele_ct = max_allele_ct;
    }
    if ((header_ctrl & 0xc0) == 0xc0) {
      header_.nonref_flags.resize(plink2::DivUp(m, plink2::kBitsPerWord) + 1);
      info.nonref_flags = header_.nonref_flags.data();
    }
    if (plink2::cachealigned_malloc(
          std::max<uintptr_t>(pgfi_cachelines, 1) * plink2::kCacheline, &header_.alloc))
      Rcpp::stop("Out of memory opening PGEN file: " + path);
    uint32_t max_vrec_width;
    uintptr_t pgr_cachelines;
    if (plink2::PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, info.raw_variant_ct, &max_vrec_width,
                               &info, header_.alloc, &pgr_cachelines, err) !=
        plink2::kPglRetSuccess)
      fail(err);
    if (info.gflags & plink2::kfPgenGlobalDosagePresent) fail_dosage();

    const uintptr_t genovec_bytes = plink2::NypCtToVecCt(n) * plink2::kBytesPerVec;
    for (int t = 0; t < threads; ++t) {
      readers_.emplace_back(new ThreadReader);
      ThreadReader& r = *readers_.back();
      plink2::PgrSetFreadBuf(nullptr, &r.pgr);
      if (plink2::cachealigned_malloc(
            std::max<uintptr_t>(pgr_cachelines, 1) * plink2::kCacheline, &r.alloc) ||
          plink2::cachealigned_malloc(genovec_bytes, &r.genovec))
        Rcpp::stop("Out of memory opening PGEN file: " + path);
      const plink2::PglErr e =
        plink2::PgrInit(path.c_str(), max_vrec_width, &info, &r.pgr, r.alloc);
      if (e != plink2::kPglRetSuccess)
        Rcpp::stop("Could not open PGEN file (pgenlib error " +
                   std::to_string(static_cast<int>(e)) + "): " + path);
      plink2::PgrClearSampleSubsetIndex(&r.pgr, &r.pssi);
    }
  }

  std::size_t samples() const { return header_.info.raw_sample_ct; }
  // Packed genovec of one variant (all samples, file order) read by one thread.
  const unsigned char* read(std::size_t variant, int thread) {
    ThreadReader& r = *readers_[thread];
    const plink2::PglErr e =
      plink2::PgrGet(nullptr, r.pssi, static_cast<uint32_t>(samples()),
                     static_cast<uint32_t>(variant), &r.pgr,
                     reinterpret_cast<uintptr_t*>(r.genovec));
    if (e != plink2::kPglRetSuccess) return nullptr;
    return r.genovec;
  }

  // ALT dosage for the geno object: genovec codes 0 hom REF, 1 het, 2 two
  // non-REF alleles, 3 missing. This matches the BED table of geno (A1
  // dosage, as BEDMatrix), since plink2 imports BED with A1 as ALT.
  static const NibbleTable& table() {
    static const NibbleTable t({0, 1, 2, -1});
    return t;
  }

private:
  [[noreturn]] void fail(const char* err) const {
    std::string msg(err);
    if (msg.compare(0, 7, "Error: ") == 0) msg.erase(0, 7);
    while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) msg.pop_back();
    Rcpp::stop("Could not open PGEN file " + path_ + (msg.empty() ? "" : ": " + msg));
  }

  [[noreturn]] void fail_dosage() const {
    Rcpp::stop("PGEN files with dosages are not supported; use hardcall PGEN files: " + path_);
  }

  std::string path_;
  PgenHeader header_;
  std::vector<std::unique_ptr<ThreadReader>> readers_;
};

} // namespace

#endif
