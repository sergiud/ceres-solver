# Compositional Euclidean norms

This directory compares the compensated variadic `AccurateNorm` and
`AccurateRNorm` against variants that compose the 2-argument functions, as
suggested during the review of the change introducing them.

## Background

Novaković [1] computes the Frobenius norm of large arrays by composing a
hypotenuse function. Two partial norms of disjoint subarrays are combined into
the norm of their concatenation using a single `hypot` call. The paper
considers two compositions:

* The sequential composition of Eq. 5 combines the norm of all preceding
  elements with the next element. Its relative error bound grows linearly with
  the number of elements.
* The recursive composition of Listing 1 splits the array into two contiguous
  halves, the first one being at most one element longer, computes the norms
  of both halves recursively, and combines them. Its relative error bound
  grows logarithmically with the number of elements.

Both compositions were applied to the variadic overloads using the existing
2-argument `AccurateNorm` as the hypotenuse function. The variadic
`AccurateRNorm` composes the norms in the same way but performs the final
combination using the 2-argument `AccurateRNorm`. The 2-argument and Jet
overloads are unchanged.

## Variants

| Variant | Revision | Description |
|---|---|---|
| `baseline` | `74fadb55` | Compensated accumulation of all squares, the base of both compositions |
| `recursive` | `44268ea1` | Recursive composition of 2-argument norms |
| `sequential` | `048d4d69` | Sequential composition of 2-argument norms |

## Method

The accuracy is measured by `accuracy.cc` using 200,000 random samples per
configuration for `float` and `double` with 3, 4, 5, 8, and 16 arguments. The
reference is the sum of squares computed exactly using MPFR whose square root
and reciprocal square root are rounded once to the target type. The error is
the distance in ULP between the result and the reference. Every argument is
drawn as m·2^e with a uniform mantissa m in [1, 2), a uniform exponent e, and a
random sign. The distributions differ in the range of the exponent:

| Distribution | Exponent range of `double` | Exponent range of `float` |
|---|---|---|
| `same-binade` | 0 | 0 |
| `spread-4` | [-4, 4] | [-4, 4] |
| `spread-wide` | [-100, 100] | [-20, 20] |
| `near-overflow` | [1014, 1022] | [118, 126] |
| `near-underflow` | [-1022, -1014] | [-126, -118] |
| `subnormal` | [-1060, -1036] | [-150, -126] |
| `full-range` | [-1020, 1020] | [-124, 124] |

The runtime is measured by `benchmark.cc` as the best time per call over two
rounds of seven repetitions, each evaluating 16,384 argument sets 40 times. The
arguments are drawn uniformly from [-10, 10]. For `double`, the `large` range
scales them by 2^600 to force rescaling.

The results below were obtained on an AMD Ryzen 9 5950X running Linux 7.2.8
using GCC 16.2.1 and Clang 23.1.1 with `-std=c++17 -O3 -march=native`, MPFR
4.2.2, and GMP 6.3.0. Other processes used about two and a half cores during
the runtime measurements. Since every round measures all variants one after
another, the ratios between the variants are more reliable than the absolute
times.

## Findings

* The compensated `AccurateNorm` and `AccurateRNorm` are correctly rounded in
  all samples of all configurations, including subnormal arguments and
  arguments near overflow.
* Both compositions lose correct rounding because they round every partial
  norm. For 3 arguments of similar magnitude, both are correctly rounded in
  83% to 94% of the samples of the norm and in 80% to 94% of the samples of
  the reciprocal norm. For 16 arguments, the share drops to 61% to 73% for the
  recursive and to 47% to 58% for the sequential composition of the norm, and
  to 58% to 90% and 50% to 81% of the reciprocal norm, respectively. The
  maximum error of the norm is 3 ULP for the recursive and 4 ULP for the
  sequential composition.
* The maximum error of the composed reciprocal norm is 3 ULP for the recursive
  and 7 ULP for the sequential composition. Subnormal `float` arguments
  increase it to 13 and 17 ULP, respectively, because the partial norms are
  rounded to the few significant bits of subnormal numbers. The reciprocal
  norm of subnormal `double` arguments overflows for all variants, which makes
  these configurations trivially exact.
* The compositions fail in simple cases. `AccurateNorm(1.0, 1.0, 1.0)` is 1 ULP
  away from √3. For four equal arguments x, the exact result 2x is
  representable, but both compositions return the value 1 ULP below it.
* The composed `AccurateRNorm` returns zero if a partial norm overflows even
  though the reciprocal of the norm is finite. For eight arguments equal to
  1.875·2^1023, both compositions return zero instead of the subnormal
  0x0.1822cb17ff2ecp-1022, which the compensated implementation returns
  exactly. This causes the large errors of the `near-overflow` reciprocal norm
  with 8 and 16 arguments.
* The compositions are not more accurate than the compensated implementation
  in any configuration. For three `double` arguments in the `spread-wide`
  distribution, for instance, the composed reciprocal norm is correctly
  rounded in 96.6% of the samples compared to all samples.
* Both compositions are slower in all configurations. The recursive
  composition is 1.8 to 3.5 times slower using GCC and 1.7 to 3.3 times slower
  using Clang. The sequential composition is 1.8 to 5.9 times slower using GCC
  and 1.6 to 5.5 times slower using Clang.

The composition targets large arrays where the logarithmic growth of the error
bound and vectorization pay off. For the small number of arguments relevant to
Ceres, the compensated accumulation is correctly rounded, avoids the overflow
of partial norms, and is faster.

## Unit tests

Both composition revisions fail the tests that expect correctly rounded
variadic results, namely `ReciprocalNormIsCorrectlyRounded`,
`VariadicNormAccuracy`, `VariadicNormIsAccurateAroundUnscaledRange`,
`VariadicReciprocalNormOfMixedSigns`,
`VariadicRNormIsAccurateAroundUnscaledRange`, and
`VariadicRNormIsAccurateForExtremeMagnitudes`. All other tests pass.

## Reproduction

The revisions listed above must be present in the repository.

```sh
cd compositional
make -j4 accuracy
make benchmark  # run on an otherwise idle machine
make report
```

## References

[1] Novaković, V. (2025). Recursive vectorized computation of the vector
    p-norm. https://arxiv.org/abs/2509.06220

## Results

### Accuracy

Each cell shows the share of correctly rounded results, the mean error, and the maximum error in ULP.

#### double norm

| n | distribution | baseline | recursive | sequential |
|---|---|---|---|---|
| 3 | full-range | 100.00% / 0.000 / 0 | 99.99% / 0.000 / 1 | 99.99% / 0.000 / 1 |
| 3 | near-overflow | 100.00% / 0.000 / 0 | 83.27% / 0.167 / 1 | 83.27% / 0.167 / 1 |
| 3 | near-underflow | 100.00% / 0.000 / 0 | 83.42% / 0.166 / 1 | 83.42% / 0.166 / 1 |
| 3 | same-binade | 100.00% / 0.000 / 0 | 82.59% / 0.174 / 1 | 82.59% / 0.174 / 1 |
| 3 | spread-4 | 100.00% / 0.000 / 0 | 83.27% / 0.167 / 1 | 83.27% / 0.167 / 1 |
| 3 | spread-wide | 100.00% / 0.000 / 0 | 99.15% / 0.009 / 1 | 99.15% / 0.009 / 1 |
| 3 | subnormal | 100.00% / 0.000 / 0 | 84.58% / 0.154 / 1 | 84.58% / 0.154 / 1 |
| 4 | full-range | 100.00% / 0.000 / 0 | 99.98% / 0.000 / 1 | 99.98% / 0.000 / 1 |
| 4 | near-overflow | 100.00% / 0.000 / 0 | 77.87% / 0.221 / 1 | 77.42% / 0.226 / 1 |
| 4 | near-underflow | 100.00% / 0.000 / 0 | 77.77% / 0.222 / 1 | 77.45% / 0.226 / 1 |
| 4 | same-binade | 100.00% / 0.000 / 0 | 78.75% / 0.212 / 1 | 74.59% / 0.254 / 1 |
| 4 | spread-4 | 100.00% / 0.000 / 0 | 77.85% / 0.222 / 1 | 77.36% / 0.226 / 1 |
| 4 | spread-wide | 100.00% / 0.000 / 0 | 98.36% / 0.016 / 1 | 98.36% / 0.016 / 1 |
| 4 | subnormal | 100.00% / 0.000 / 0 | 77.31% / 0.227 / 1 | 78.06% / 0.219 / 1 |
| 5 | full-range | 100.00% / 0.000 / 0 | 99.97% / 0.000 / 1 | 99.97% / 0.000 / 1 |
| 5 | near-overflow | 100.00% / 0.000 / 0 | 75.03% / 0.250 / 1 | 72.93% / 0.272 / 2 |
| 5 | near-underflow | 100.00% / 0.000 / 0 | 74.95% / 0.251 / 1 | 72.97% / 0.271 / 2 |
| 5 | same-binade | 100.00% / 0.000 / 0 | 73.99% / 0.260 / 1 | 69.36% / 0.307 / 2 |
| 5 | spread-4 | 100.00% / 0.000 / 0 | 74.95% / 0.251 / 1 | 73.03% / 0.271 / 2 |
| 5 | spread-wide | 100.00% / 0.000 / 0 | 97.45% / 0.025 / 1 | 97.44% / 0.026 / 1 |
| 5 | subnormal | 100.00% / 0.000 / 0 | 73.46% / 0.265 / 1 | 73.13% / 0.270 / 2 |
| 8 | full-range | 100.00% / 0.000 / 0 | 99.92% / 0.001 / 1 | 99.92% / 0.001 / 1 |
| 8 | near-overflow | 100.00% / 0.000 / 0 | 72.86% / 0.271 / 2 | 64.15% / 0.370 / 3 |
| 8 | near-underflow | 100.00% / 0.000 / 0 | 72.98% / 0.270 / 2 | 64.23% / 0.368 / 3 |
| 8 | same-binade | 100.00% / 0.000 / 0 | 82.21% / 0.178 / 2 | 70.88% / 0.293 / 2 |
| 8 | spread-4 | 100.00% / 0.000 / 0 | 72.83% / 0.272 / 2 | 63.80% / 0.373 / 3 |
| 8 | spread-wide | 100.00% / 0.000 / 0 | 94.25% / 0.058 / 1 | 94.27% / 0.057 / 2 |
| 8 | subnormal | 100.00% / 0.000 / 0 | 67.94% / 0.321 / 2 | 63.75% / 0.374 / 3 |
| 16 | full-range | 100.00% / 0.000 / 0 | 99.65% / 0.004 / 1 | 99.66% / 0.003 / 1 |
| 16 | near-overflow | 100.00% / 0.000 / 0 | 71.09% / 0.290 / 2 | 51.01% / 0.546 / 4 |
| 16 | near-underflow | 100.00% / 0.000 / 0 | 70.93% / 0.292 / 2 | 50.84% / 0.548 / 4 |
| 16 | same-binade | 100.00% / 0.000 / 0 | 72.72% / 0.273 / 2 | 47.15% / 0.586 / 3 |
| 16 | spread-4 | 100.00% / 0.000 / 0 | 70.90% / 0.292 / 2 | 51.03% / 0.546 / 4 |
| 16 | spread-wide | 100.00% / 0.000 / 0 | 85.39% / 0.146 / 2 | 85.53% / 0.145 / 2 |
| 16 | subnormal | 100.00% / 0.000 / 0 | 61.41% / 0.391 / 2 | 50.21% / 0.562 / 4 |

#### double rnorm

| n | distribution | baseline | recursive | sequential |
|---|---|---|---|---|
| 3 | full-range | 100.00% / 0.000 / 0 | 99.66% / 0.003 / 1 | 99.66% / 0.003 / 1 |
| 3 | near-overflow | 100.00% / 0.000 / 0 | 86.37% / 0.136 / 1 | 86.37% / 0.136 / 1 |
| 3 | near-underflow | 100.00% / 0.000 / 0 | 83.19% / 0.168 / 1 | 83.19% / 0.168 / 1 |
| 3 | same-binade | 100.00% / 0.000 / 0 | 80.34% / 0.197 / 1 | 80.34% / 0.197 / 1 |
| 3 | spread-4 | 100.00% / 0.000 / 0 | 83.24% / 0.168 / 1 | 83.24% / 0.168 / 1 |
| 3 | spread-wide | 100.00% / 0.000 / 0 | 96.59% / 0.034 / 1 | 96.59% / 0.034 / 1 |
| 3 | subnormal | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 |
| 4 | full-range | 100.00% / 0.000 / 0 | 99.55% / 0.004 / 1 | 99.34% / 0.007 / 2 |
| 4 | near-overflow | 100.00% / 0.000 / 0 | 83.15% / 0.169 / 1 | 82.90% / 0.172 / 2 |
| 4 | near-underflow | 100.00% / 0.000 / 0 | 77.67% / 0.223 / 1 | 77.27% / 0.230 / 2 |
| 4 | same-binade | 100.00% / 0.000 / 0 | 81.65% / 0.183 / 1 | 77.61% / 0.224 / 1 |
| 4 | spread-4 | 100.00% / 0.000 / 0 | 77.72% / 0.223 / 1 | 77.25% / 0.230 / 2 |
| 4 | spread-wide | 100.00% / 0.000 / 0 | 95.41% / 0.046 / 1 | 93.72% / 0.063 / 2 |
| 4 | subnormal | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 |
| 5 | full-range | 100.00% / 0.000 / 0 | 99.35% / 0.006 / 1 | 98.99% / 0.010 / 1 |
| 5 | near-overflow | 100.00% / 0.000 / 0 | 82.59% / 0.175 / 2 | 80.99% / 0.193 / 3 |
| 5 | near-underflow | 100.00% / 0.000 / 0 | 75.03% / 0.251 / 2 | 72.97% / 0.278 / 3 |
| 5 | same-binade | 100.00% / 0.000 / 0 | 81.60% / 0.184 / 1 | 78.14% / 0.219 / 2 |
| 5 | spread-4 | 100.00% / 0.000 / 0 | 75.03% / 0.251 / 2 | 72.90% / 0.278 / 3 |
| 5 | spread-wide | 100.00% / 0.000 / 0 | 93.50% / 0.065 / 2 | 91.05% / 0.090 / 3 |
| 5 | subnormal | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 |
| 8 | full-range | 100.00% / 0.000 / 0 | 98.81% / 0.012 / 2 | 97.93% / 0.021 / 2 |
| 8 | near-overflow | 100.00% / 0.000 / 0 | 84.61% / 0.154 / 2 | 79.84% / 5250866989.401 / 1050173397838403 |
| 8 | near-underflow | 100.00% / 0.000 / 0 | 72.33% / 0.279 / 2 | 63.55% / 0.390 / 5 |
| 8 | same-binade | 100.00% / 0.000 / 0 | 74.98% / 0.250 / 2 | 59.53% / 0.417 / 3 |
| 8 | spread-4 | 100.00% / 0.000 / 0 | 72.35% / 0.279 / 2 | 63.38% / 0.392 / 4 |
| 8 | spread-wide | 100.00% / 0.000 / 0 | 89.81% / 0.102 / 2 | 84.77% / 0.154 / 3 |
| 8 | subnormal | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 |
| 16 | full-range | 100.00% / 0.000 / 0 | 97.52% / 0.025 / 2 | 95.52% / 0.045 / 2 |
| 16 | near-overflow | 100.00% / 0.000 / 0 | 89.56% / 35889577468.955 / 1100266561135794 | 80.44% / 2519219614094.072 / 1125635333639789 |
| 16 | near-underflow | 100.00% / 0.000 / 0 | 69.91% / 0.304 / 3 | 49.99% / 0.587 / 6 |
| 16 | same-binade | 100.00% / 0.000 / 0 | 76.48% / 0.235 / 1 | 52.68% / 0.503 / 3 |
| 16 | spread-4 | 100.00% / 0.000 / 0 | 69.91% / 0.303 / 3 | 50.14% / 0.586 / 6 |
| 16 | spread-wide | 100.00% / 0.000 / 0 | 81.92% / 0.182 / 2 | 74.66% / 0.261 / 4 |
| 16 | subnormal | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 | 100.00% / 0.000 / 0 |

#### float norm

| n | distribution | baseline | recursive | sequential |
|---|---|---|---|---|
| 3 | full-range | 100.00% / 0.000 / 0 | 99.87% / 0.001 / 1 | 99.87% / 0.001 / 1 |
| 3 | near-overflow | 100.00% / 0.000 / 0 | 83.22% / 0.168 / 1 | 83.22% / 0.168 / 1 |
| 3 | near-underflow | 100.00% / 0.000 / 0 | 83.32% / 0.167 / 1 | 83.32% / 0.167 / 1 |
| 3 | same-binade | 100.00% / 0.000 / 0 | 82.58% / 0.174 / 1 | 82.58% / 0.174 / 1 |
| 3 | spread-4 | 100.00% / 0.000 / 0 | 83.31% / 0.167 / 1 | 83.31% / 0.167 / 1 |
| 3 | spread-wide | 100.00% / 0.000 / 0 | 95.84% / 0.042 / 1 | 95.84% / 0.042 / 1 |
| 3 | subnormal | 100.00% / 0.000 / 0 | 94.02% / 0.060 / 2 | 94.02% / 0.060 / 2 |
| 4 | full-range | 100.00% / 0.000 / 0 | 99.75% / 0.003 / 1 | 99.73% / 0.003 / 1 |
| 4 | near-overflow | 100.00% / 0.000 / 0 | 77.70% / 0.223 / 1 | 77.44% / 0.226 / 1 |
| 4 | near-underflow | 100.00% / 0.000 / 0 | 77.77% / 0.222 / 1 | 77.41% / 0.226 / 1 |
| 4 | same-binade | 100.00% / 0.000 / 0 | 78.68% / 0.213 / 1 | 74.50% / 0.255 / 1 |
| 4 | spread-4 | 100.00% / 0.000 / 0 | 77.70% / 0.223 / 1 | 77.54% / 0.225 / 1 |
| 4 | spread-wide | 100.00% / 0.000 / 0 | 92.75% / 0.072 / 1 | 92.88% / 0.071 / 1 |
| 4 | subnormal | 100.00% / 0.000 / 0 | 88.00% / 0.121 / 2 | 88.12% / 0.120 / 2 |
| 5 | full-range | 100.00% / 0.000 / 0 | 99.55% / 0.004 / 1 | 99.55% / 0.004 / 1 |
| 5 | near-overflow | 100.00% / 0.000 / 0 | 75.17% / 0.248 / 1 | 72.91% / 0.272 / 2 |
| 5 | near-underflow | 100.00% / 0.000 / 0 | 75.10% / 0.249 / 2 | 72.92% / 0.272 / 2 |
| 5 | same-binade | 100.00% / 0.000 / 0 | 73.95% / 0.261 / 1 | 69.53% / 0.305 / 2 |
| 5 | spread-4 | 100.00% / 0.000 / 0 | 74.82% / 0.252 / 1 | 72.93% / 0.272 / 2 |
| 5 | spread-wide | 100.00% / 0.000 / 0 | 89.83% / 0.102 / 1 | 89.84% / 0.102 / 2 |
| 5 | subnormal | 100.00% / 0.000 / 0 | 82.93% / 0.172 / 2 | 82.82% / 0.175 / 2 |
| 8 | full-range | 100.00% / 0.000 / 0 | 98.83% / 0.012 / 1 | 98.84% / 0.012 / 1 |
| 8 | near-overflow | 100.00% / 0.000 / 0 | 72.70% / 0.273 / 2 | 64.05% / 0.371 / 3 |
| 8 | near-underflow | 100.00% / 0.000 / 0 | 72.82% / 0.272 / 2 | 64.14% / 0.370 / 3 |
| 8 | same-binade | 100.00% / 0.000 / 0 | 82.21% / 0.178 / 2 | 70.73% / 0.294 / 2 |
| 8 | spread-4 | 100.00% / 0.000 / 0 | 73.00% / 0.270 / 2 | 64.29% / 0.369 / 3 |
| 8 | spread-wide | 100.00% / 0.000 / 0 | 82.61% / 0.174 / 2 | 82.33% / 0.177 / 3 |
| 8 | subnormal | 100.00% / 0.000 / 0 | 73.16% / 0.274 / 2 | 71.98% / 0.291 / 3 |
| 16 | full-range | 100.00% / 0.000 / 0 | 96.20% / 0.038 / 1 | 96.18% / 0.038 / 2 |
| 16 | near-overflow | 100.00% / 0.000 / 0 | 71.23% / 0.288 / 2 | 50.90% / 0.547 / 4 |
| 16 | near-underflow | 100.00% / 0.000 / 0 | 71.00% / 0.291 / 2 | 50.88% / 0.547 / 4 |
| 16 | same-binade | 100.00% / 0.000 / 0 | 72.87% / 0.272 / 2 | 47.37% / 0.582 / 3 |
| 16 | spread-4 | 100.00% / 0.000 / 0 | 71.15% / 0.289 / 2 | 50.98% / 0.547 / 4 |
| 16 | spread-wide | 100.00% / 0.000 / 0 | 72.99% / 0.271 / 2 | 69.90% / 0.308 / 3 |
| 16 | subnormal | 100.00% / 0.000 / 0 | 65.05% / 0.365 / 3 | 58.20% / 0.465 / 4 |

#### float rnorm

| n | distribution | baseline | recursive | sequential |
|---|---|---|---|---|
| 3 | full-range | 100.00% / 0.000 / 0 | 98.64% / 0.014 / 1 | 98.64% / 0.014 / 1 |
| 3 | near-overflow | 100.00% / 0.000 / 0 | 86.23% / 0.138 / 1 | 86.23% / 0.138 / 1 |
| 3 | near-underflow | 100.00% / 0.000 / 0 | 83.36% / 0.166 / 1 | 83.36% / 0.166 / 1 |
| 3 | same-binade | 100.00% / 0.000 / 0 | 80.28% / 0.197 / 1 | 80.28% / 0.197 / 1 |
| 3 | spread-4 | 100.00% / 0.000 / 0 | 83.34% / 0.167 / 1 | 83.34% / 0.167 / 1 |
| 3 | spread-wide | 100.00% / 0.000 / 0 | 92.02% / 0.080 / 1 | 92.02% / 0.080 / 1 |
| 3 | subnormal | 100.00% / 0.000 / 0 | 94.02% / 0.075 / 5 | 94.02% / 0.075 / 5 |
| 4 | full-range | 100.00% / 0.000 / 0 | 98.17% / 0.018 / 1 | 97.39% / 0.026 / 2 |
| 4 | near-overflow | 100.00% / 0.000 / 0 | 83.25% / 0.168 / 1 | 82.97% / 0.172 / 2 |
| 4 | near-underflow | 100.00% / 0.000 / 0 | 77.96% / 0.220 / 1 | 77.52% / 0.227 / 2 |
| 4 | same-binade | 100.00% / 0.000 / 0 | 81.74% / 0.183 / 1 | 77.62% / 0.224 / 1 |
| 4 | spread-4 | 100.00% / 0.000 / 0 | 77.73% / 0.223 / 1 | 77.28% / 0.230 / 2 |
| 4 | spread-wide | 100.00% / 0.000 / 0 | 89.26% / 0.107 / 1 | 86.87% / 0.132 / 2 |
| 4 | subnormal | 100.00% / 0.000 / 0 | 88.89% / 0.140 / 5 | 87.21% / 0.168 / 8 |
| 5 | full-range | 100.00% / 0.000 / 0 | 97.37% / 0.026 / 2 | 96.17% / 0.038 / 2 |
| 5 | near-overflow | 100.00% / 0.000 / 0 | 82.47% / 0.176 / 2 | 81.00% / 0.193 / 3 |
| 5 | near-underflow | 100.00% / 0.000 / 0 | 74.77% / 0.254 / 2 | 72.67% / 0.281 / 3 |
| 5 | same-binade | 100.00% / 0.000 / 0 | 81.43% / 0.186 / 1 | 78.04% / 0.220 / 2 |
| 5 | spread-4 | 100.00% / 0.000 / 0 | 75.14% / 0.250 / 2 | 72.88% / 0.278 / 3 |
| 5 | spread-wide | 100.00% / 0.000 / 0 | 86.00% / 0.140 / 2 | 82.73% / 0.175 / 3 |
| 5 | subnormal | 100.00% / 0.000 / 0 | 82.97% / 0.221 / 8 | 81.09% / 0.259 / 10 |
| 8 | full-range | 100.00% / 0.000 / 0 | 95.62% / 0.044 / 2 | 92.94% / 0.071 / 2 |
| 8 | near-overflow | 100.00% / 0.000 / 0 | 84.74% / 0.153 / 2 | 79.65% / 10.671 / 2092187 |
| 8 | near-underflow | 100.00% / 0.000 / 0 | 72.23% / 0.280 / 2 | 63.22% / 0.394 / 4 |
| 8 | same-binade | 100.00% / 0.000 / 0 | 75.02% / 0.250 / 2 | 59.51% / 0.417 / 3 |
| 8 | spread-4 | 100.00% / 0.000 / 0 | 72.47% / 0.277 / 2 | 63.60% / 0.390 / 4 |
| 8 | spread-wide | 100.00% / 0.000 / 0 | 80.42% / 0.197 / 2 | 74.94% / 0.257 / 4 |
| 8 | subnormal | 100.00% / 0.000 / 0 | 70.97% / 0.389 / 9 | 67.97% / 0.475 / 14 |
| 16 | full-range | 100.00% / 0.000 / 0 | 91.22% / 0.088 / 2 | 85.71% / 0.144 / 3 |
| 16 | near-overflow | 100.00% / 0.000 / 0 | 89.58% / 103.583 / 2071575 | 80.50% / 4370.347 / 2096309 |
| 16 | near-underflow | 100.00% / 0.000 / 0 | 70.00% / 0.303 / 2 | 50.12% / 0.587 / 6 |
| 16 | same-binade | 100.00% / 0.000 / 0 | 76.52% / 0.235 / 2 | 52.51% / 0.505 / 3 |
| 16 | spread-4 | 100.00% / 0.000 / 0 | 70.07% / 0.302 / 2 | 49.98% / 0.589 / 7 |
| 16 | spread-wide | 100.00% / 0.000 / 0 | 72.76% / 0.277 / 3 | 64.69% / 0.377 / 4 |
| 16 | subnormal | 100.00% / 0.000 / 0 | 58.34% / 0.559 / 13 | 50.88% / 0.769 / 17 |

### Runtime

Each cell shows the best time per call in nanoseconds.

#### clang

| type | n | range | function | baseline | recursive | sequential |
|---|---|---|---|---|---|---|
| double | 3 | large | norm | 4.96 | 9.23 | 9.18 |
| double | 4 | large | norm | 6.12 | 12.76 | 15.75 |
| double | 8 | large | norm | 12.91 | 39.48 | 60.20 |
| double | 16 | large | norm | 28.91 | 95.32 | 157.69 |
| double | 3 | typical | norm | 4.67 | 7.80 | 7.58 |
| double | 4 | typical | norm | 5.63 | 11.57 | 14.53 |
| double | 8 | typical | norm | 10.99 | 35.31 | 55.14 |
| double | 16 | typical | norm | 30.06 | 83.70 | 144.67 |
| double | 3 | large | rnorm | 5.51 | 10.21 | 9.90 |
| double | 4 | large | rnorm | 7.94 | 14.19 | 16.23 |
| double | 8 | large | rnorm | 12.62 | 39.16 | 62.79 |
| double | 16 | large | rnorm | 31.49 | 94.37 | 159.92 |
| double | 3 | typical | rnorm | 5.01 | 9.23 | 9.00 |
| double | 4 | typical | rnorm | 6.05 | 12.71 | 14.98 |
| double | 8 | typical | rnorm | 10.85 | 34.02 | 57.66 |
| double | 16 | typical | rnorm | 31.33 | 80.32 | 146.02 |
| float | 3 | typical | norm | 3.87 | 6.99 | 7.05 |
| float | 4 | typical | norm | 4.91 | 9.97 | 12.29 |
| float | 8 | typical | norm | 10.62 | 26.13 | 46.68 |
| float | 16 | typical | norm | 29.70 | 69.90 | 120.72 |
| float | 3 | typical | rnorm | 4.51 | 7.70 | 7.23 |
| float | 4 | typical | rnorm | 5.70 | 10.63 | 12.60 |
| float | 8 | typical | rnorm | 10.69 | 27.88 | 47.43 |
| float | 16 | typical | rnorm | 34.03 | 69.23 | 121.33 |

#### gcc

| type | n | range | function | baseline | recursive | sequential |
|---|---|---|---|---|---|---|
| double | 3 | large | norm | 4.90 | 9.82 | 9.89 |
| double | 4 | large | norm | 6.45 | 14.80 | 15.86 |
| double | 8 | large | norm | 12.96 | 43.65 | 58.35 |
| double | 16 | large | norm | 45.38 | 99.19 | 156.36 |
| double | 3 | typical | norm | 4.54 | 8.91 | 8.77 |
| double | 4 | typical | norm | 5.92 | 12.54 | 14.82 |
| double | 8 | typical | norm | 12.12 | 38.98 | 54.16 |
| double | 16 | typical | norm | 42.75 | 88.13 | 142.99 |
| double | 3 | large | rnorm | 5.31 | 10.80 | 10.61 |
| double | 4 | large | rnorm | 7.04 | 15.76 | 16.66 |
| double | 8 | large | rnorm | 13.45 | 41.92 | 63.35 |
| double | 16 | large | rnorm | 28.01 | 96.83 | 161.55 |
| double | 3 | typical | rnorm | 4.96 | 9.74 | 9.56 |
| double | 4 | typical | rnorm | 6.43 | 14.14 | 15.29 |
| double | 8 | typical | rnorm | 12.43 | 37.56 | 57.73 |
| double | 16 | typical | rnorm | 24.80 | 85.79 | 146.96 |
| float | 3 | typical | norm | 3.94 | 6.99 | 7.21 |
| float | 4 | typical | norm | 5.18 | 10.15 | 12.34 |
| float | 8 | typical | norm | 11.78 | 31.49 | 46.10 |
| float | 16 | typical | norm | 40.40 | 73.86 | 119.35 |
| float | 3 | typical | rnorm | 4.38 | 8.39 | 8.07 |
| float | 4 | typical | rnorm | 5.67 | 11.73 | 12.87 |
| float | 8 | typical | rnorm | 11.70 | 30.43 | 48.07 |
| float | 16 | typical | rnorm | 24.08 | 70.05 | 121.88 |
