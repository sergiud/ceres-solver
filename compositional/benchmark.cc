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
// Measures the throughput of the variadic AccurateNorm and AccurateRNorm. The
// program prints one CSV row per floating-point type, number of arguments, and
// magnitude range with the best time per call in nanoseconds.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

#include "ceres/accurate_norm.h"

namespace {

constexpr std::size_t kSets = 1 << 14;
constexpr int kRepetitions = 7;
constexpr int kPasses = 40;
constexpr double kMagnitude = 10;

template <typename T, std::size_t N, std::size_t... I>
T Norm(const std::array<T, N>& values, std::index_sequence<I...>) {
  return ceres::AccurateNorm(values[I]...);
}

template <typename T, std::size_t N, std::size_t... I>
T RNorm(const std::array<T, N>& values, std::index_sequence<I...>) {
  return ceres::AccurateRNorm(values[I]...);
}

// Accumulates the results to prevent the compiler from eliminating the calls.
double sink = 0;

// Returns the best time per call in nanoseconds.
template <typename T, std::size_t N, typename Function>
double Time(const std::vector<std::array<T, N>>& inputs, Function function) {
  double best = std::numeric_limits<double>::infinity();

  for (int repetition = 0; repetition < kRepetitions; ++repetition) {
    T sum = 0;
    const auto start = std::chrono::steady_clock::now();

    for (int pass = 0; pass < kPasses; ++pass) {
      for (const auto& values : inputs) {
        sum += function(values);
      }
    }

    const std::chrono::duration<double, std::nano> elapsed =
        std::chrono::steady_clock::now() - start;
    sink += static_cast<double>(sum);
    best = std::min(best,
                    elapsed.count() / (static_cast<double>(kSets) * kPasses));
  }

  return best;
}

// Draws values uniformly from [-10, 10] scaled by 2^exponent.
template <typename T, std::size_t N>
void Run(const char* type, const char* range, int exponent) {
  std::mt19937_64 generator(N);
  std::uniform_real_distribution<double> uniform(-kMagnitude, kMagnitude);
  std::vector<std::array<T, N>> inputs(kSets);

  for (auto& values : inputs) {
    for (T& value : values) {
      value = static_cast<T>(std::ldexp(uniform(generator), exponent));
    }
  }

  const double norm = Time<T, N>(inputs, [](const std::array<T, N>& values) {
    return Norm(values, std::make_index_sequence<N>{});
  });
  const double rnorm = Time<T, N>(inputs, [](const std::array<T, N>& values) {
    return RNorm(values, std::make_index_sequence<N>{});
  });
  std::printf("%s,%zu,%s,%.3f,%.3f\n", type, N, range, norm, rnorm);
}

template <typename T>
void RunAll(const char* type) {
  // Forces the rescaling of double arguments.
  constexpr int kLargeExponent = 600;

  const auto run = [type](auto count) {
    constexpr std::size_t kCount = decltype(count)::value;
    Run<T, kCount>(type, "typical", 0);
    if constexpr (std::is_same_v<T, double>) {
      Run<T, kCount>(type, "large", kLargeExponent);
    }
  };

  run(std::integral_constant<std::size_t, 3>{});
  run(std::integral_constant<std::size_t, 4>{});
  run(std::integral_constant<std::size_t, 8>{});
  run(std::integral_constant<std::size_t, 16>{});
}

}  // namespace

int main() {
  std::printf("type,n,range,norm_ns,rnorm_ns\n");
  RunAll<double>("double");
  RunAll<float>("float");
  std::fprintf(stderr, "%g\n", sink);
}
