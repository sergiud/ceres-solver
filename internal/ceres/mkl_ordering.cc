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

#include "ceres/mkl_ordering.h"

#ifndef CERES_NO_MKL

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/event_logger.h"
#include "ceres/mkl_diagnostics.h"
#include "ceres/mkl_normal_matrix.h"
#include "ceres/mkl_pardiso.h"
#include "ceres/mkl_sparse_matrix.h"
#include "ceres/types.h"
#include "mkl.h"

namespace ceres::internal {

bool MklComputeOrdering(const CompressedRowSparseMatrix& matrix,
                        const LinearSolverOrderingType ordering_type,
                        const int max_num_threads,
                        int* ordering,
                        std::string* message) {
  EventLogger event_logger("MklComputeOrdering");
  std::unique_ptr<CompressedRowSparseMatrix> normal_matrix;
  if (!ComputeAtAUsingMkl(matrix, max_num_threads, &normal_matrix, message)) {
    return false;
  }
  event_logger.AddEvent("Form J'J");
  const bool success = ComputePardisoOrdering(
      *normal_matrix, ordering_type, max_num_threads, ordering, message);
  event_logger.AddEvent("PARDISO Ordering");
  return success;
}

bool MklComputeSchurOrdering(const CompressedRowSparseMatrix& e_matrix,
                             const CompressedRowSparseMatrix& f_matrix,
                             const LinearSolverOrderingType ordering_type,
                             const int max_num_threads,
                             int* ordering,
                             std::string* message) {
  EventLogger event_logger("MklComputeSchurOrdering");
  const MklThreadScope thread_scope(max_num_threads);
  MklCsrMatrix e_handle;
  MklCsrMatrix f_handle;
  if (!e_handle.Create(e_matrix, message) ||
      !f_handle.Create(f_matrix, message)) {
    return false;
  }
  matrix_descr descriptor{};
  descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
  descriptor.mode = SPARSE_FILL_MODE_FULL;
  descriptor.diag = SPARSE_DIAG_NON_UNIT;
  sparse_matrix_t raw_etf_handle = nullptr;
  if (!CheckMklStatus(mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE,
                                      descriptor,
                                      e_handle.get(),
                                      SPARSE_OPERATION_NON_TRANSPOSE,
                                      descriptor,
                                      f_handle.get(),
                                      SPARSE_STAGE_FULL_MULT,
                                      &raw_etf_handle),
                      "sparse Schur product",
                      message)) {
    return false;
  }
  MklSparseHandle etf_handle(raw_etf_handle);
  event_logger.AddEvent("Compute E'F");

  std::unique_ptr<CompressedRowSparseMatrix> etf;
  if (!ExportMklCsr(etf_handle.get(),
                    CompressedRowSparseMatrix::StorageType::UNSYMMETRIC,
                    &etf,
                    message)) {
    return false;
  }
  std::fill_n(etf->mutable_values(), etf->num_nonzeros(), 1.0);

  std::unique_ptr<CompressedRowSparseMatrix> etf_normal;
  std::unique_ptr<CompressedRowSparseMatrix> f_normal;
  if (!ComputeAtAUsingMkl(*etf, max_num_threads, &etf_normal, message) ||
      !ComputeAtAUsingMkl(f_matrix, max_num_threads, &f_normal, message)) {
    return false;
  }
  event_logger.AddEvent("Form Schur pattern");

  const int dimension = f_matrix.num_cols();
  CompressedRowSparseMatrix combined(
      dimension,
      dimension,
      etf_normal->num_nonzeros() + f_normal->num_nonzeros());
  int nonzero = 0;
  combined.mutable_rows()[0] = 0;
  std::vector<int> columns;
  for (int row = 0; row < dimension; ++row) {
    const int etf_row_size =
        etf_normal->rows()[row + 1] - etf_normal->rows()[row];
    const int f_row_size = f_normal->rows()[row + 1] - f_normal->rows()[row];
    columns.clear();
    columns.reserve(static_cast<std::size_t>(etf_row_size + f_row_size));
    for (int index = etf_normal->rows()[row];
         index < etf_normal->rows()[row + 1];
         ++index) {
      if (etf_normal->cols()[index] >= row) {
        columns.push_back(etf_normal->cols()[index]);
      }
    }
    for (int index = f_normal->rows()[row]; index < f_normal->rows()[row + 1];
         ++index) {
      if (f_normal->cols()[index] >= row) {
        columns.push_back(f_normal->cols()[index]);
      }
    }
    std::sort(columns.begin(), columns.end());
    columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
    for (const int column : columns) {
      combined.mutable_cols()[nonzero] = column;
      combined.mutable_values()[nonzero] = 1.0;
      ++nonzero;
    }
    combined.mutable_rows()[row + 1] = nonzero;
  }
  combined.set_storage_type(
      CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);

  const bool success = ComputePardisoOrdering(
      combined, ordering_type, max_num_threads, ordering, message);
  event_logger.AddEvent("PARDISO Ordering");
  return success;
}

}  // namespace ceres::internal

#endif  // CERES_NO_MKL
