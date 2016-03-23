#ifndef CERES_INTERNAL_PARFOR_H
#define CERES_INTERNAL_PARFOR_H

#if defined(CERES_USE_OPENMP)
#include <omp.h>

namespace ceres {
namespace internal {

template <typename Var, typename Body>
void parfor(Var first, Var last, Var step, int n_threads, Body body) {
#pragma omp parallel for num_threads(n_threads)
  for (Var i = first; i < last; i += step) body(omp_get_thread_num(), i);
}

template <typename Var, typename Body>
void parfor_dynamic(Var first, Var last, Var step, int n_threads, Body body) {
#pragma omp parallel for num_threads(n_threads) schedule(dynamic)
  for (Var i = first; i < last; i += step) body(omp_get_thread_num(), i);
}

}  // namespace internal
}  // namespace ceres

#elif defined(CERES_USE_TBB)
#include <tbb/tbb.h>
#include <cassert>

namespace ceres {
namespace internal {

template <typename Var, typename Body>
void parfor(Var first, Var last, Var step, int n_threads, Body body) {
  assert(n_threads > 0 && "parfor n_threads must be positive");
  assert(first <= last && step > 0 && "parfor step must be positive");
  tbb::parallel_for(0, n_threads, 1, [&](int thread_id) {
    for (Var i = first + thread_id; i < last; i += step * n_threads)
      body(thread_id, i);
  });
}

template <typename Var, typename Body>
void parfor_dynamic(Var first, Var last, Var step, int threads, Body body) {
  parfor(first, last, step, threads, body);
}

}  // namespace internal
}  // namespace ceres

#else  // Sequential

#define CERES_PARFOR_THREAD_ID() 0

namespace ceres {
namespace internal {

template <typename Var, typename Body>
void parfor(Var first, Var last, Var step, int /*threads*/, Body body) {
  for (Var i = first; i < last; i += step) body(0, i);
}

template <typename Var, typename Body>
void parfor_dynamic(Var first, Var last, Var step, int n_threads, Body body) {
  parfor(first, last, step, n_threads, body);
}

}  // namespace internal
}  // namespace ceres
#endif

#endif  // CERES_INTERNAL_PARFOR_H
