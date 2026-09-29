// 6R.284 live capture. Compiles GKL 0.8.8 AVX float and double kernels and
// runs the IntelPairHmm.cc branch unchanged:
//   result_float = g_compute_full_prob_float
//   if (result_float < MIN_ACCEPTED) double kernel
//   else (double)(log10f(result_float) - LOG10_INITIAL_CONSTANT)
// No arithmetic in this file other than that branch.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gnu/libc-version.h>
#include <immintrin.h>

#include "avx_impl.h"
#include "pairhmm_common.h"
#include "Context.h"

namespace {

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : line) {
        if (ch == '\t') {
            out.push_back(cur);
            cur.clear();
        } else if (ch != '\r') {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    return out;
}

std::vector<char> parse_bytes(const std::string& csv) {
    std::vector<char> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) {
            continue;
        }
        int v = std::atoi(item.c_str());
        out.push_back(static_cast<char>(v));
    }
    return out;
}

struct ReadRec {
    int ri;
    int skipped;
    std::string bases;
    std::vector<char> bq;
    std::vector<char> iq;
    std::vector<char> dq;
    std::vector<char> gcp;
};

void print_bits_f32(const char* key, float x) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &x, 4);
    std::printf("META\t%s\t0x%08x\t%.9g\n", key, bits, static_cast<double>(x));
}

void print_bits_f64(const char* key, double x) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, 8);
    std::printf("META\t%s\t0x%016llx\t%.17g\n", key, static_cast<unsigned long long>(bits), x);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: 6r284_live_capture inputs.tsv\n");
        return 2;
    }
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    ConvertChar::init();
    Context<float> ctxf;
    Context<double> ctxd;

    std::printf("META\tuname_machine\tx86_64\n");
    std::printf("META\tglibc\t%s\n", gnu_get_libc_version());
    std::printf("META\tavx\t%d\n", __builtin_cpu_supports("avx") ? 1 : 0);
    std::printf("META\tavx2\t%d\n", __builtin_cpu_supports("avx2") ? 1 : 0);
    std::printf("META\tavx512f\t%d\n", __builtin_cpu_supports("avx512f") ? 1 : 0);
    std::printf("META\tsimd_path\tAVX compute_fp_avxs / compute_fp_avxd\n");
    std::printf("META\tuse_double_precision\tfalse\n");
    std::printf("META\tftz\tON\n");
    std::printf("META\tdaz\tunchanged\n");
    print_bits_f32("log10_initial_f32", Context<float>::LOG10_INITIAL_CONSTANT);
    print_bits_f64("log10_initial_f64", Context<double>::LOG10_INITIAL_CONSTANT);
    (void)ctxf;
    (void)ctxd;

    std::vector<std::string> haps;
    std::vector<ReadRec> reads;
    std::ifstream in(argv[1]);
    if (!in) {
        std::perror(argv[1]);
        return 1;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        auto f = split_tab(line);
        if (f[0] == "HAP") {
            haps.push_back(f[3]);
        } else if (f[0] == "READ") {
            ReadRec r;
            r.ri = std::atoi(f[1].c_str());
            r.skipped = std::atoi(f[4].c_str());
            r.bases = f[6];
            r.bq = parse_bytes(f[7]);
            r.iq = parse_bytes(f[8]);
            r.dq = parse_bytes(f[9]);
            r.gcp = parse_bytes(f[10]);
            if (r.bq.size() != r.bases.size() || r.iq.size() != r.bases.size() ||
                r.dq.size() != r.bases.size() || r.gcp.size() != r.bases.size()) {
                std::fprintf(stderr, "length mismatch read %d\n", r.ri);
                return 1;
            }
            reads.push_back(std::move(r));
        }
    }
    if (haps.size() != 5) {
        std::fprintf(stderr, "expected 5 haplotypes, got %zu\n", haps.size());
        return 1;
    }

    const float min_accepted = MIN_ACCEPTED;
    for (const auto& r : reads) {
        for (int k = 0; k < 5; k++) {
            testcase tc;
            tc.rs = r.bases.data();
            tc.rslen = static_cast<int>(r.bases.size());
            tc.hap = haps[static_cast<std::size_t>(k)].data();
            tc.haplen = static_cast<int>(haps[static_cast<std::size_t>(k)].size());
            tc.q = r.bq.data();
            tc.i = r.iq.data();
            tc.d = r.dq.data();
            tc.c = r.gcp.data();

            float result_float = compute_fp_avxs(&tc);
            const bool below = result_float < min_accepted;
            double result_final;
            const char* branch;
            if (below) {
                double result_double = compute_fp_avxd(&tc);
                result_final = std::log10(result_double) - Context<double>::LOG10_INITIAL_CONSTANT;
                branch = "double";
            } else {
                result_final = static_cast<double>(::log10f(result_float) -
                                                   Context<float>::LOG10_INITIAL_CONSTANT);
                branch = "float";
            }
            std::uint32_t fb = 0;
            std::uint64_t db = 0;
            std::memcpy(&fb, &result_float, 4);
            std::memcpy(&db, &result_final, 8);
            std::printf(
                "CELL\t%d\t%d\t%d\t%s\t0x%08x\t%.9g\t0x%016llx\t%.17g\n",
                r.ri,
                k,
                r.skipped,
                branch,
                fb,
                static_cast<double>(result_float),
                static_cast<unsigned long long>(db),
                result_final);
        }
    }
    return 0;
}
