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
//
// [1] Borges, C. F. (2021). Fast Compensated Algorithms for the Reciprocal
//     Square Root, the Reciprocal Hypotenuse, and Givens Rotations.
//     http://arxiv.org/abs/2103.08694
//
// [2] Borges, C. F. (2021). Algorithm 1014: An Improved Algorithm for
//     hypot(x,y). ACM Transactions on Mathematical Software, 47(1), 1–12.
//     https://doi.org/10.1145/3428446

#ifndef CERES_PUBLIC_ACCURATE_NORM_
#define CERES_PUBLIC_ACCURATE_NORM_

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

namespace ceres {

namespace internal {

template <typename T, typename Enable = void>
struct Promote {};

template <typename T>
struct Promote<T, std::enable_if_t<std::is_integral_v<T>>> {
  using type = double;
};

template <typename T>
struct Promote<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  using type = T;
};

template <typename... Ts>
using Promote_t = decltype((typename Promote<Ts>::type(0) + ... + 0));

// Determines the unit in the last place (ulp) of a value x. ulp is the spacing
// between consecutive floating-point numbers. For instance, ulp(1) = b^(1-p) =
// ε corresponds to the machine epsilon for a floating-point type with radix b
// and precision p.
template <typename T>
constexpr auto Ulp(T x) -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  using std::fpclassify;
  using std::ilogb;
  using std::scalbn;

  const int cls = fpclassify(x);

  if (cls == FP_NAN) {
    return x;
  }

  if (cls == FP_INFINITE) {
    // Return ∞ for ±∞ since ulp(-x) = ulp(x)
    return std::numeric_limits<T>::infinity();
  }

  if (cls == FP_NORMAL) {
    // Compute b^(e-p+1) = ε·b^e where e = ⌊log_b x⌋ is the logarithm to base
    // b of x, b is the radix and p is the floating-point type precision.
    return scalbn(std::numeric_limits<T>::epsilon(), ilogb(x));
  }

  // Subnormal or zero
  return std::numeric_limits<T>::min();
}

// Compute the ulp while promoting the argument type to floating-point, e.g., if
// the function is called as Ulp(0) where 0 is integer literal.
template <typename T>
constexpr auto Ulp(T x)
    -> std::enable_if_t<!std::is_floating_point_v<T>, Promote_t<T>> {
  return Ulp(static_cast<Promote_t<T>>(x));
}

// To ensure compile-time evaluation we need std::sqrt to be constexpr which is
// the case since C++26
#if defined(__cpp_lib_constexpr_cmath) && (__cpp_lib_constexpr_cmath >= 202306L)
#define CERES_HAS_CONSTEXPR_CMATH26
#endif

// GCC is non-conforming in regard to constexpr support as the compiler supports
// constexpr cmath using builtin functions without raising the C++ standard to
// C++26. The constexpr cmath support is available unless the code is built
// using -fno-builtin.
#if (defined(__GNUG__) && !defined(__clang__)) && defined(__has_builtin)
#if __has_builtin(sqrt) && __has_builtin(scalbn) && __has_builtin(ilogb)
#define CERES_HAS_CONSTEXPR_FOR_ACCURATENORMTRAITS
#endif
#endif

#if defined(CERES_HAS_CONSTEXPR_CMATH26) && \
    !defined(CERES_HAS_CONSTEXPR_FOR_ACCURATENORMTRAITS)
#define CERES_HAS_CONSTEXPR_FOR_ACCURATENORMTRAITS
#endif

#if defined(CERES_HAS_CONSTEXPR_FOR_ACCURATENORMTRAITS)
#define CERES_ACCURATENORM_CONSTEXPR constexpr
#else
#define CERES_ACCURATENORM_CONSTEXPR const
#endif

template <typename T, typename Enable = void>
struct AccurateNormTraits {
  static constexpr T Varying() noexcept {
    using std::sqrt;
    return sqrt(std::numeric_limits<T>::epsilon() / 2);
  }

  static constexpr T Huge() noexcept {
    using std::sqrt;
    return sqrt(std::numeric_limits<T>::max() / 2);
  }

  static constexpr T Tiny() noexcept {
    using std::sqrt;
    return sqrt(std::numeric_limits<T>::min());
  }

  static constexpr T Scale() noexcept { return Ulp(Tiny()); }
};

// Unless Ceres is compiled with extended constexpr support for cmath introduced
// in C++26, provide a specialization for IEEE-754 double arithmetic that avoids
// computing the thresholds at runtime.
#if !defined(CERES_HAS_CONSTEXPR_FOR_ACCURATENORMTRAITS)
template <>
struct AccurateNormTraits<
    double,
    std::enable_if<std::numeric_limits<double>::is_iec559>> {
  // ulp(√ε/2)
  static constexpr double Varying() noexcept { return 0x1.6a09e667f3bcdp-27; }
  // √F_max/2
  static constexpr double Huge() noexcept { return 0x1.6a09e667f3bccp+511; }
  // √F_min
  static constexpr double Tiny() noexcept { return 0x1p-511; }
  // ulp(√F_min)
  static constexpr double Scale() noexcept { return 0x1p-563; }
};
#endif  // !defined(CERES_HAS_CONSTEXPR_FOR_ACCURATENORMTRAITS)

// Compute two values s, t that satisfy s + t = x + y exactly where s is the sum
// nearest to x + y and t is the round-off error. The algorithm assumes the
// round-to-nearest (RN) mode which is the default.
template <typename T>
constexpr auto Fast2Sum(T x, T y)
    -> std::enable_if_t<std::is_floating_point_v<T>, std::pair<T, T>> {
  using std::fabs;
  using std::isgreaterequal;
  assert(isgreaterequal(fabs(x), fabs(y)));
  const T s = x + y;
  const T z = s - x;
  const T t = y - z;
  return std::make_pair(s, t);
}

template <typename T>
constexpr auto UnscaledAccurateNormWithError(T x, T y)
    -> std::enable_if_t<std::is_floating_point_v<T>, std::pair<T, T>> {
  using std::fma;
  using std::sqrt;

  const T x_sq = x * x;
  const T y_sq = y * y;
  // Recover the rounding error of the floating-point addition of both squares.
  const auto [sigma, sigma_e] = Fast2Sum(x_sq, y_sq);
  // Use the 2MultFMA algorithm to recover the rounding error due to squaring x
  // and y.
  return std::make_pair(sigma, sigma_e + fma(y, y, -y_sq) + fma(x, x, -x_sq));
}

// Computes the hypotenuse of x and y without checking the arguments and
// ensuring invariants. Not intended to be invoked by users.
//
// The functions assumes the arguments to be finite, scaled correctly to avoid
// an under-/overflow and passed in the correct order ensuring |x| ≥ |y|.
template <typename T>
constexpr auto UnscaledAccurateNorm(T x, T y)
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  using std::fma;
  using std::sqrt;

  const auto [sigma, sigma_e] = UnscaledAccurateNormWithError(x, y);
  const T h = sqrt(sigma);
  const T tau = sigma_e + fma(-h, h, sigma);
  return fma(tau / h, T(0.5), h);
}

template <typename T>
constexpr auto UnscaledAccurateRNorm(T x, T y)
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  using std::fma;
  using std::sqrt;

  auto [sigma, sigma_e] = UnscaledAccurateNormWithError(x, y);
  const T r = T(1) / sigma;
  sigma = fma(-r, sigma_e, fma(-r, sigma, T(1)));
  const T rho = sqrt(r);
  const T tau = fma(-rho, rho, r);
  const T nu = fma(sigma, tau, sigma) / 2;
  return fma(rho, nu, rho);
}

}  // namespace internal

template <typename T>
constexpr auto AccurateNorm(T a, T b)
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  using std::fabs;
  using std::fmax;
  using std::fmin;
  using std::isfinite;

  a = fabs(a);
  b = fabs(b);

  if (!isfinite(a)) {
    return a;
  }

  if (!isfinite(b)) {
    return b;
  }

  // Ensure |x| ≥ |y|. While we could have used std::minmax, we want to avoid
  // floating point exceptions.
  const T x = fmax(a, b);
  const T y = fmin(a, b);

  using internal::AccurateNormTraits;

  if (y <= x * AccurateNormTraits<T>::Varying()) {
    return x;
  }

  using internal::UnscaledAccurateNorm;

  CERES_ACCURATENORM_CONSTEXPR T scale = AccurateNormTraits<T>::Scale();

  if (x > AccurateNormTraits<T>::Huge()) {
    // Scale x to prevent an overflow
    return UnscaledAccurateNorm(x * scale, y * scale) / scale;
  }

  if (y < AccurateNormTraits<T>::Tiny()) {
    // Scale y to prevent an underflow
    return UnscaledAccurateNorm(x / scale, y / scale) * scale;
  }

  // Avoid rounding errors due to unnecessary scaling
  return UnscaledAccurateNorm(x, y);
}

template <typename T, typename... Args>
constexpr auto AccurateNorm(T a, T b, Args&&... args)
    -> std::enable_if_t<(sizeof...(Args) > 0 &&
                         (std::is_same_v<T, std::decay_t<Args>> && ...)),
                        T> {
  return AccurateNorm(a, AccurateNorm(b, std::forward<Args>(args)...));
}

template <typename T, typename U, typename... Args>
constexpr internal::Promote_t<T, U, Args...> AccurateNorm(T a,
                                                          U b,
                                                          Args&&... args) {
  using Type = internal::Promote_t<T, U, Args...>;
  return AccurateNorm(Type(a), Type(b), Type(std::forward<Args>(args))...);
}

template <typename T>
constexpr auto AccurateRNorm(T a, T b)
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  using std::fabs;
  using std::fmax;
  using std::fmin;
  using std::fpclassify;
  using std::isfinite;
  using std::isnan;

  a = fabs(a);
  b = fabs(b);

  // Ensure |x| ≥ |y|. While we could have used std::minmax, we want to avoid
  // floating point exceptions.
  const T x = fmax(a, b);
  const T y = fmin(a, b);

  const int cls = fpclassify(x);

  if (cls == FP_NAN) {
    // FIXME when y not +/-oo
    return x;
  }

  if (cls == FP_INFINITE) {
    return T(0);  // +/-oo, y
  }

  if (cls == FP_ZERO) {
    return std::numeric_limits<T>::quiet_NaN();
  }

  using internal::AccurateNormTraits;

  if (y <= x * AccurateNormTraits<T>::Varying()) {
    return T(1) / x;
  }

  using internal::UnscaledAccurateRNorm;

  CERES_ACCURATENORM_CONSTEXPR T scale = AccurateNormTraits<T>::Scale();

  // The rescaling differs from the one used in AccurateNorm because scaling the
  // arguments x and y of a reciprocal hypotenuse yields
  //
  //     1/sqrt(x'^2+y'^2)
  // <=> 1/sqrt((x*s)^2+(y*s)^2)
  // <=> 1/(s*sqrt(x^2+y^2))
  //
  // i.e., to cancel the scale, we need reapply it to the result.

  if (x > AccurateNormTraits<T>::Huge()) {
    // Scale x to prevent an overflow
    return UnscaledAccurateRNorm(x * scale, y * scale) * scale;
  }

  if (y < AccurateNormTraits<T>::Tiny()) {
    // Scale y to prevent an underflow
    return UnscaledAccurateRNorm(x / scale, y / scale) / scale;
  }

  // Avoid rounding errors due to unnecessary scaling
  return UnscaledAccurateRNorm(x, y);
}

template <typename T, typename... Args>
constexpr auto AccurateRNorm(T a, T b, Args&&... args)
    -> std::enable_if_t<(sizeof...(Args) > 0 &&
                         (std::is_same_v<T, std::decay_t<Args>> && ...)),
                        T> {
  // Note that we compose the reciprocal hypotenuse with the non-reciprocal one as this
  // is the convention of the arguments. Additionally, this avoids division by
  // zero in cases such AccurateRNorm(x, 0, 0).
  return AccurateRNorm(a, AccurateNorm(b, std::forward<Args>(args)...));
}

template <typename T, typename U, typename... Args>
constexpr internal::Promote_t<T, U, Args...> AccurateRNorm(T a,
                                                           U b,
                                                           Args&&... args) {
  using Type = internal::Promote_t<T, U, Args...>;
  return AccurateRNorm(Type(a), Type(b), Type(std::forward<Args>(args))...);
}

}  // namespace ceres

#endif  // CERES_PUBLIC_ACCURATE_NORM_
