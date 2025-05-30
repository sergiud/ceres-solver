// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2025 Google Inc. All rights reserved.
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

#include "ceres/accurate_norm.h"

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "gtest/gtest.h"

namespace {

#if 0
template <typename T>
auto KahanSum1(T a, T b)
    -> std::enable_if_t<std::is_floating_point_v<T>, std::pair<T, T>> {
  using std::fmax;
  using std::fmin;

  const T x = fmax(a, b);
  const T y = fmin(a, b);

  return Fast2Sum(x, y);
}

template <typename T>
auto KahanSum1(T a, const std::pair<T, T>& st)
    -> std::enable_if_t<std::is_floating_point_v<T>, std::pair<T, T>> {
  const auto& [s, t] = st;
  return Fast2Sum(a - t, s);
}

template<typename T>
auto KahanSum1(T a) -> std::enable_if_t<std::is_floating_point_v<T>, std::pair<T, T>>
{
    return std::make_pair(a, T(0));
}

template<typename T, typename ...Ts>
auto KahanSum(T a, T b, Ts&& ...args) -> std::enable_if_t<(std::is_same_v<T, std::decay_t<Ts>> && ... && true), T>
{
    return KahanSum1(a, KahanSum1(b,  std::forward<Ts>(args)...)).first;
}
#endif

template <typename T>
constexpr auto FloatDistance(T a, T b)
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  using std::copysign;
  using std::fabs;
  using std::fmin;
  using std::fpclassify;
  using std::ilogb;
  using std::isgreater;
  using std::isless;
  using std::scalbn;
  using std::signbit;

  if (isgreater(a, b)) {
    return -FloatDistance(b, a);
  }

  if (fpclassify(a - b) == FP_ZERO) {
    return T{0};
  }

  if (fpclassify(a) == FP_ZERO) {
    return T{1} + fabs(FloatDistance(
                      copysign(std::numeric_limits<T>::denorm_min(), b), b));
  }

  if (fpclassify(b) == FP_ZERO) {
    return T{1} + fabs(FloatDistance(
                      copysign(std::numeric_limits<T>::denorm_min(), a), a));
  }

  if (signbit(a) != signbit(b)) {
    return T{2} +
           fabs(FloatDistance(copysign(std::numeric_limits<T>::denorm_min(), b),
                              b)) +
           fabs(FloatDistance(copysign(std::numeric_limits<T>::denorm_min(), a),
                              a));
  }

  // a, b are either both positive or both negative. Above we already ensure a <
  // b.
  if (isless(a, 0)) {
    return FloatDistance(-b, -a);
  }

  int e1 = ilogb(a) + 1;
  const T upper1 = scalbn(T{1}, e1);

  T result{0};

  if (isgreater(b, upper1)) {
    const int e2 = ilogb(b);
    const T upper2 = scalbn(T{1}, e2);

    result = FloatDistance(upper2, b) +
             scalbn(e2 - e1, std::numeric_limits<T>::digits - 1);
  }

  e1 = std::numeric_limits<T>::digits - e1;

  const T mb = -fmin(upper1, b);

  using ceres::internal::Fast2Sum;

  const auto [s, t] = Fast2Sum(a, mb);
  const T x = s;
  const T y = Fast2Sum(b, mb - t).first;

  return result + scalbn(fabs(x), e1) + scalbn(fabs(y), e1);
}

}  // namespace

TEST(AccurateNorm, Promote) {
  static_assert(std::is_same_v<ceres::internal::Promote_t<int>, double>,
                "Promotion of an int must be a double");
  static_assert(std::is_same_v<ceres::internal::Promote_t<int, int>, double>,
                "Promotion of multuple ints must be a double");
  static_assert(std::is_same_v<ceres::internal::Promote_t<unsigned>, double>,
                "Promotion of an unsigned int must be double");
  static_assert(std::is_same_v<ceres::internal::Promote_t<long>, double>,
                "Promotion of a long must be double");
  static_assert(
      std::is_same_v<ceres::internal::Promote_t<int, long, float>, double>,
      "Promotion of arithmetic types must be double");
}

#if GTEST_HAS_TYPED_TEST

template <typename T>
class AccurateNormTest : public testing::Test {
 public:
  static constexpr auto kTiny = std::numeric_limits<T>::min();
  static constexpr auto kHuge = std::numeric_limits<T>::max();
};

using Types = testing::Types<float, double, long double>;

TYPED_TEST_SUITE(AccurateNormTest, Types);

TYPED_TEST(AccurateNormTest, Ulp) {
  using Scalar = TypeParam;

  EXPECT_TRUE(std::isnan(
      ceres::internal::Ulp(std::numeric_limits<Scalar>::quiet_NaN())));
  EXPECT_EQ(std::fpclassify(
                ceres::internal::Ulp(+std::numeric_limits<Scalar>::infinity())),
            FP_INFINITE);
  EXPECT_EQ(std::fpclassify(
                ceres::internal::Ulp(-std::numeric_limits<Scalar>::infinity())),
            FP_INFINITE);
  EXPECT_EQ(ceres::internal::Ulp(Scalar{0}),
            std::numeric_limits<Scalar>::min());

  EXPECT_EQ(ceres::internal::Ulp(Scalar{+1}),
            std::numeric_limits<Scalar>::epsilon());
  EXPECT_EQ(ceres::internal::Ulp(Scalar{-1}),
            std::numeric_limits<Scalar>::epsilon());
}

TYPED_TEST(AccurateNormTest, Norm) {
  using Scalar = TypeParam;

  EXPECT_EQ(ceres::AccurateNorm(this->kTiny, Scalar{0}), this->kTiny);
  EXPECT_EQ(ceres::AccurateNorm(Scalar{0}, this->kTiny), this->kTiny);

  EXPECT_EQ(ceres::AccurateNorm(this->kHuge, Scalar{0}), this->kHuge);
  EXPECT_EQ(ceres::AccurateNorm(Scalar{0}, this->kHuge), this->kHuge);

  EXPECT_EQ(ceres::AccurateNorm(this->kTiny, this->kTiny),
            this->kTiny * std::sqrt(Scalar{2}));

  EXPECT_EQ(ceres::AccurateNorm(Scalar{0}, Scalar{0}), 0);

  EXPECT_TRUE(std::isinf(
      ceres::AccurateNorm(+std::numeric_limits<Scalar>::infinity(), 0)));
  EXPECT_TRUE(std::isinf(
      ceres::AccurateNorm(-std::numeric_limits<Scalar>::infinity(), 0)));

  EXPECT_TRUE(std::isinf(
      ceres::AccurateNorm(0, +std::numeric_limits<Scalar>::infinity())));
  EXPECT_TRUE(std::isinf(
      ceres::AccurateNorm(0, -std::numeric_limits<Scalar>::infinity())));

  EXPECT_TRUE(std::isnan(
      ceres::AccurateNorm(std::numeric_limits<Scalar>::quiet_NaN(), 0)));
  EXPECT_TRUE(std::isnan(
      ceres::AccurateNorm(0, std::numeric_limits<Scalar>::quiet_NaN())));
}

TYPED_TEST(AccurateNormTest, RNorm) {
  using Scalar = TypeParam;

  EXPECT_EQ(ceres::AccurateRNorm(this->kTiny, Scalar{0}), 1 / this->kTiny);
  EXPECT_EQ(ceres::AccurateRNorm(Scalar{0}, this->kTiny), 1 / this->kTiny);

  EXPECT_EQ(ceres::AccurateRNorm(this->kHuge, Scalar{0}), 1 / this->kHuge);
  EXPECT_EQ(ceres::AccurateRNorm(Scalar{0}, this->kHuge), 1 / this->kHuge);

  EXPECT_TRUE(std::isnan(ceres::AccurateRNorm(0, 0)));

  const auto large = std::sqrt(this->kHuge / 2);
  const auto a = ceres::AccurateRNorm(large, large);
  const auto expected = 1 / std::sqrt(this->kHuge);

  const auto d = FloatDistance(a, expected);
  EXPECT_LE(d, 1);

  EXPECT_EQ(ceres::AccurateRNorm(+std::numeric_limits<Scalar>::infinity(), 0),
            0);
  EXPECT_EQ(ceres::AccurateRNorm(-std::numeric_limits<Scalar>::infinity(), 0),
            0);

  EXPECT_EQ(ceres::AccurateRNorm(0, +std::numeric_limits<Scalar>::infinity()),
            0);
  EXPECT_EQ(ceres::AccurateRNorm(0, -std::numeric_limits<Scalar>::infinity()),
            0);

  EXPECT_TRUE(std::isnan(
      ceres::AccurateRNorm(std::numeric_limits<Scalar>::quiet_NaN(), 0)));
  EXPECT_TRUE(std::isnan(
      ceres::AccurateRNorm(0, std::numeric_limits<Scalar>::quiet_NaN())));
}

#endif
