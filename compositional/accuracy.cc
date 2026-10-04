// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2026 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sergiu.deitsch@gmail.com (Sergiu Deitsch)
//
// Measures the accuracy of the variadic AccurateNorm and AccurateRNorm against
// correctly rounded references computed using MPFR. The program prints one CSV
// row per floating-point type, number of arguments, magnitude distribution,
// and function.

#include <mpfr.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>

#include "ceres/accurate_norm.h"

namespace {

// Precision that represents the sum of squares of all tested arguments exactly.
constexpr mpfr_prec_t kReferencePrecision = 5000;
constexpr long kSamples = 200000;

// Constants for deriving distinct random seeds from the configuration.
constexpr std::size_t kCountSeedFactor = 1000003;
constexpr int kCenterSeedFactor = 7919;

constexpr double kPercent = 100;

// Returns the correctly rounded norm and reciprocal norm of the values.
template <typename T, std::size_t N>
std::pair<T, T> Reference(const std::array<T, N>& values) {
  mpfr_t sum;
  mpfr_t square;
  mpfr_t result;
  mpfr_init2(sum, kReferencePrecision);
  mpfr_init2(square, kReferencePrecision);
  mpfr_init2(result, std::numeric_limits<T>::digits);
  mpfr_set_zero(sum, 1);

  for (const T value : values) {
    mpfr_set_ld(square, static_cast<long double>(value), MPFR_RNDN);
    mpfr_sqr(square, square, MPFR_RNDN);
    mpfr_add(sum, sum, square, MPFR_RNDN);
  }

  mpfr_sqrt(result, sum, MPFR_RNDN);
  const T norm = static_cast<T>(mpfr_get_ld(result, MPFR_RNDN));
  mpfr_rec_sqrt(result, sum, MPFR_RNDN);
  const T rnorm = static_cast<T>(mpfr_get_ld(result, MPFR_RNDN));
  mpfr_clears(sum, square, result, static_cast<mpfr_ptr>(nullptr));
  return {norm, rnorm};
}

// Counts the floating-point values between two non-negative values using their
// binary representations, which are ordered like the values themselves.
template <typename T>
long UlpDistance(T a, T b) {
  using Integer = std::conditional_t<sizeof(T) == sizeof(std::int64_t),
                                     std::int64_t,
                                     std::int32_t>;
  static_assert(sizeof(Integer) == sizeof(T));

  Integer x;
  Integer y;
  std::memcpy(&x, &a, sizeof(T));
  std::memcpy(&y, &b, sizeof(T));
  return static_cast<long>(x > y ? x - y : y - x);
}

template <typename T, std::size_t N, std::size_t... I>
T Norm(const std::array<T, N>& values, std::index_sequence<I...>) {
  return ceres::AccurateNorm(values[I]...);
}

template <typename T, std::size_t N, std::size_t... I>
T RNorm(const std::array<T, N>& values, std::index_sequence<I...>) {
  return ceres::AccurateRNorm(values[I]...);
}

struct Statistics {
  long count = 0;
  long correctly_rounded = 0;
  long above_one_ulp = 0;
  long maximum = 0;
  double total = 0;

  void Add(long distance) {
    ++count;
    if (distance == 0) {
      ++correctly_rounded;
    }
    if (distance > 1) {
      ++above_one_ulp;
    }
    if (distance > maximum) {
      maximum = distance;
    }
    total += static_cast<double>(distance);
  }
};

void Print(const char* type,
           std::size_t count,
           const char* distribution,
           const char* function,
           const Statistics& statistics) {
  std::printf("%s,%zu,%s,%s,%.4f,%.4f,%ld,%.4f\n",
              type,
              count,
              distribution,
              function,
              kPercent * statistics.correctly_rounded / statistics.count,
              statistics.total / statistics.count,
              statistics.maximum,
              kPercent * statistics.above_one_ulp / statistics.count);
}

// Draws values m * 2^e with a uniformly distributed mantissa m in [1, 2), a
// uniformly distributed exponent e in [center - spread, center + spread], and a
// random sign.
template <typename T, std::size_t N>
void Run(const char* type, const char* distribution, int center, int spread) {
  std::mt19937_64 generator(N * kCountSeedFactor + center * kCenterSeedFactor +
                            spread);
  std::uniform_real_distribution<double> mantissa(1, 2);
  std::uniform_int_distribution<int> exponent(center - spread, center + spread);
  std::bernoulli_distribution negative(0.5);
  Statistics norm;
  Statistics rnorm;

  for (long sample = 0; sample < kSamples; ++sample) {
    std::array<T, N> values;

    // Draw the random numbers in a fixed order since the evaluation order of
    // function arguments is unspecified.
    for (T& value : values) {
      const double m = mantissa(generator);
      const int e = exponent(generator);
      const bool sign = negative(generator);
      value = static_cast<T>(std::ldexp(sign ? -m : m, e));
    }

    const auto [expected_norm, expected_rnorm] = Reference(values);
    norm.Add(UlpDistance(Norm(values, std::make_index_sequence<N>{}),
                         expected_norm));
    rnorm.Add(UlpDistance(RNorm(values, std::make_index_sequence<N>{}),
                          expected_rnorm));
  }

  Print(type, N, distribution, "norm", norm);
  Print(type, N, distribution, "rnorm", rnorm);
}

template <typename T>
void RunAll(const char* type, int wide_spread, int extreme_exponent) {
  constexpr int kNarrowSpread = 4;
  constexpr int kSubnormalSpread = 12;
  // Centers the subnormal distribution such that most values are subnormal.
  constexpr int kSubnormalCenter = std::numeric_limits<T>::min_exponent -
                                   std::numeric_limits<T>::digits / 2 - 1;
  // Keeps the largest drawn magnitude finite.
  constexpr int kFullSpread = std::numeric_limits<T>::max_exponent - 4;

  const auto run = [type, wide_spread, extreme_exponent](auto count) {
    constexpr std::size_t kCount = decltype(count)::value;
    Run<T, kCount>(type, "same-binade", 0, 0);
    Run<T, kCount>(type, "spread-4", 0, kNarrowSpread);
    Run<T, kCount>(type, "spread-wide", 0, wide_spread);
    Run<T, kCount>(type, "near-overflow", extreme_exponent, kNarrowSpread);
    Run<T, kCount>(type, "near-underflow", -extreme_exponent, kNarrowSpread);
    Run<T, kCount>(type, "subnormal", kSubnormalCenter, kSubnormalSpread);
    Run<T, kCount>(type, "full-range", 0, kFullSpread);
  };

  run(std::integral_constant<std::size_t, 3>{});
  run(std::integral_constant<std::size_t, 4>{});
  run(std::integral_constant<std::size_t, 5>{});
  run(std::integral_constant<std::size_t, 8>{});
  run(std::integral_constant<std::size_t, 16>{});
}

}  // namespace

int main() {
  constexpr int kDoubleWideSpread = 100;
  constexpr int kDoubleExtremeExponent = 1018;
  constexpr int kFloatWideSpread = 20;
  constexpr int kFloatExtremeExponent = 122;

  std::printf(
      "type,n,distribution,function,correct_pct,mean_ulp,max_ulp,"
      "above_1ulp_pct\n");
  RunAll<double>("double", kDoubleWideSpread, kDoubleExtremeExponent);
  RunAll<float>("float", kFloatWideSpread, kFloatExtremeExponent);
}
