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

#include "ceres/mkl_pardiso.h"

#ifndef CERES_NO_MKL

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/log/vlog_is_on.h"
#include "absl/strings/str_format.h"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/event_logger.h"
#include "ceres/mkl_pardiso_diagnostics.h"
#include "ceres/mkl_sparse_matrix.h"
#include "ceres/types.h"
#include "mkl.h"

namespace ceres::internal {

namespace {

std::array<MKL_INT, kPardisoParameterCount> CopyPardisoIparm(
    const MKL_INT* iparm) {
  std::array<MKL_INT, kPardisoParameterCount> copy{};
  std::copy_n(iparm, kPardisoParameterCount, copy.begin());
  return copy;
}

}  // namespace

// PARDISO-compatible CSR matrix and the map back to the Ceres values.
struct PardisoMatrixState {
  bool defined = false;
  std::vector<MKL_INT> rows;
  std::vector<MKL_INT> columns;
  std::vector<int> source_indices;
  std::vector<double> values;
  // Dimensions of the source matrix that source_indices refers to.
  int input_num_rows = 0;
  int input_num_nonzeros = 0;
};

class CERES_NO_EXPORT PardisoSolver final {
 public:
  PardisoSolver(OrderingType ordering_type,
                int max_num_threads,
                bool use_two_level_factorization)
      : ordering_type_(ordering_type),
        max_num_threads_(max_num_threads),
        use_two_level_factorization_(use_two_level_factorization) {}

  ~PardisoSolver() {
    if (initialized_) {
      Call(kPardisoRelease, nullptr, nullptr, nullptr, nullptr);
    }
  }

  // Convert Ceres upper triangular storage to PARDISO scalar CSR.
  bool DefineStructure(const CompressedRowSparseMatrix& matrix,
                       std::string* message) {
    CHECK_EQ(matrix.storage_type(),
             CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
    if (matrix_.defined) {
      // SparseCholesky requires the sparsity structure to stay fixed after
      // the first factorization. Only the dimensions are verified because
      // RefreshValues depends on them.
      if (matrix.num_rows() != matrix_.input_num_rows ||
          matrix.num_nonzeros() != matrix_.input_num_nonzeros) {
        *message = absl::StrFormat(
            "PARDISO matrix structure changed after symbolic analysis: got "
            "%d rows and %d nonzeros, expected %d rows and %d nonzeros.",
            matrix.num_rows(),
            matrix.num_nonzeros(),
            matrix_.input_num_rows,
            matrix_.input_num_nonzeros);
        return false;
      }
      return true;
    }

    std::vector<MKL_INT> rows(matrix.num_rows() + 1);
    std::vector<MKL_INT> columns;
    std::vector<int> source_indices;
    columns.reserve(matrix.num_nonzeros() + matrix.num_rows());
    source_indices.reserve(matrix.num_nonzeros() + matrix.num_rows());
    std::vector<std::pair<int, int>> entries;
    for (int row = 0; row < matrix.num_rows(); ++row) {
      rows[row] = static_cast<MKL_INT>(columns.size());
      entries.clear();
      for (int index = matrix.rows()[row]; index < matrix.rows()[row + 1];
           ++index) {
        // Block upper triangular matrices store their diagonal blocks in full,
        // so entries below the diagonal are skipped.
        if (matrix.cols()[index] >= row) {
          entries.emplace_back(matrix.cols()[index], index);
        }
      }
      // PARDISO requires sorted columns and an explicit diagonal, which is
      // the first entry of an upper triangular row.
      std::sort(entries.begin(), entries.end());
      if (entries.empty() || entries.front().first != row) {
        entries.insert(entries.begin(), {row, -1});
      }
      for (const auto& [column, index] : entries) {
        columns.push_back(column);
        source_indices.push_back(index);
      }
    }
    rows.back() = static_cast<MKL_INT>(columns.size());

    Initialize();
    matrix_.rows = std::move(rows);
    matrix_.columns = std::move(columns);
    matrix_.source_indices = std::move(source_indices);
    matrix_.values.resize(matrix_.columns.size());
    matrix_.input_num_rows = matrix.num_rows();
    matrix_.input_num_nonzeros = matrix.num_nonzeros();
    if (ordering_type_ == OrderingType::NATURAL) {
      permutation_.resize(matrix.num_rows());
      std::iota(permutation_.begin(), permutation_.end(), 0);
    }
    matrix_.defined = true;
    return true;
  }

  // Refresh the PARDISO value array from the source matrix. Fill-in entries
  // introduced by DefineStructure stay zero.
  void RefreshValues(const CompressedRowSparseMatrix& matrix) {
    for (int index = 0; index < static_cast<int>(matrix_.source_indices.size());
         ++index) {
      if (matrix_.source_indices[index] >= 0) {
        matrix_.values[index] = matrix.values()[matrix_.source_indices[index]];
      }
    }
  }

  LinearSolverTerminationType AnalyzeStructure(std::string* message) {
    if (analyzed_) {
      return LinearSolverTerminationType::SUCCESS;
    }
    if (!matrix_.defined) {
      *message =
          "PARDISO symbolic analysis requires a defined matrix "
          "structure.";
      return LinearSolverTerminationType::FATAL_ERROR;
    }

    const MKL_INT status = Call(
        kPardisoAnalysis,
        matrix_.values.data(),
        ordering_type_ == OrderingType::NATURAL ? permutation_.data() : nullptr,
        nullptr,
        nullptr);
    CheckPardisoStatus(status, "symbolic analysis", message);
    analyzed_ = status == PARDISO_NO_ERROR;
    return PardisoErrorToTerminationType(status);
  }

  LinearSolverTerminationType Factorize(std::string* message) {
    if (!analyzed_) {
      *message =
          "PARDISO factorization requires symbolic analysis of the matrix "
          "structure first.";
      return LinearSolverTerminationType::FATAL_ERROR;
    }
    const MKL_INT status = Call(
        kPardisoNumericalFactorization,
        matrix_.values.data(),
        ordering_type_ == OrderingType::NATURAL ? permutation_.data() : nullptr,
        nullptr,
        nullptr);
    CheckPardisoStatus(status, "factorization", message);
    return PardisoErrorToTerminationType(status);
  }

  LinearSolverTerminationType Solve(const double* rhs,
                                    double* solution,
                                    std::string* message) {
    const MKL_INT status = Call(kPardisoSolve, nullptr, nullptr, rhs, solution);
    CheckPardisoStatus(status, "solve", message);
    return PardisoErrorToTerminationType(status);
  }

  bool ComputeOrdering(const CompressedRowSparseMatrix& matrix,
                       int* ordering,
                       std::string* message) {
    if (ordering_type_ == OrderingType::NATURAL) {
      std::iota(ordering, ordering + matrix.num_cols(), 0);
      return true;
    }
    if (!DefineStructure(matrix, message)) {
      return false;
    }
    iparm_[kPardisoUserPermutationParameter] = kPardisoReturnPermutation;
    std::vector<MKL_INT> permutation(matrix.num_cols());
    const MKL_INT status =
        Call(kPardisoAnalysis, nullptr, permutation.data(), nullptr, nullptr);
    if (!CheckPardisoStatus(status, "ordering", message)) {
      return false;
    }
    std::vector<bool> seen(matrix.num_cols(), false);
    for (int index = 0; index < matrix.num_cols(); ++index) {
      if (permutation[index] < 0 || permutation[index] >= matrix.num_cols()) {
        *message = absl::StrFormat(
            "PARDISO ordering at index %d is %d, expected a value in the "
            "range [0, %d).",
            index,
            permutation[index],
            matrix.num_cols());
        return false;
      }
      if (seen[permutation[index]]) {
        *message = absl::StrFormat(
            "PARDISO ordering contains duplicate value %d at index %d, "
            "expected a permutation of [0, %d).",
            permutation[index],
            index,
            matrix.num_cols());
        return false;
      }
      seen[permutation[index]] = true;
      ordering[index] = static_cast<int>(permutation[index]);
    }
    return true;
  }

 private:
  void Initialize() {
    constexpr MKL_INT matrix_type = kPardisoPositiveDefinite;
    // PARDISOINIT sets documented defaults before overrides. See:
    // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/intel-onemkl-pardiso-parameters-in-tabular-form.html
    PARDISOINIT(pparam_, &matrix_type, iparm_);
    iparm_[kPardisoUseDefaultValuesParameter] = kPardisoUserParameters;
    iparm_[kPardisoIndexBaseParameter] = kPardisoZeroBasedIndexing;
    iparm_[kPardisoIterativeRefinementParameter] =
        kPardisoAutomaticIterativeRefinement;
    iparm_[kPardisoParallelFactorizationParameter] =
        use_two_level_factorization_ ? kPardisoTwoLevelFactorization
                                     : kPardisoClassicFactorization;
    iparm_[kPardisoParallelSolveParameter] = kPardisoParallelSolve;
#ifndef NDEBUG
    iparm_[kPardisoMatrixCheckerParameter] = kPardisoEnableMatrixChecker;
#endif
    if (ordering_type_ == OrderingType::NESDIS) {
      iparm_[kPardisoOrderingParameter] = kPardisoUseNestedDissectionOrdering;
    } else if (ordering_type_ == OrderingType::AMD) {
      iparm_[kPardisoOrderingParameter] = kPardisoUseMinimumDegreeOrdering;
    } else if (ordering_type_ == OrderingType::NATURAL) {
      iparm_[kPardisoUserPermutationParameter] = kPardisoUseUserPermutation;
    }
    if (VLOG_IS_ON(3)) {
      iparm_[kPardisoFactorNonzerosParameter] = kPardisoDebugStatistics;
      iparm_[kPardisoFactorOperationsParameter] = kPardisoDebugStatistics;
    }
    initialized_ = true;
    if (VLOG_IS_ON(3)) {
      previous_iparm_ = CopyPardisoIparm(iparm_);
      has_iparm_snapshot_ = true;
      const std::vector<PardisoIparmEntry>& entries = GetPardisoIparmEntries();
      std::vector<int> indices;
      indices.reserve(entries.size());
      for (const PardisoIparmEntry& entry : entries) {
        indices.push_back(entry.index);
      }
      VLOG(3) << FormatPardisoIparmTable(
          "PARDISO iparm initialization", indices, previous_iparm_, nullptr);
    }
  }

  MKL_INT Call(MKL_INT phase,
               const double* values,
               MKL_INT* permutation,
               const double* rhs,
               double* solution) {
    const MklThreadScope thread_scope(max_num_threads_);
    MKL_INT max_factors = 1;
    MKL_INT factor_number = 1;
    MKL_INT rhs_count = 1;
    MKL_INT error = 0;
    constexpr MKL_INT matrix_type = kPardisoPositiveDefinite;
    const MKL_INT size = matrix_.rows.empty() ? 0 : matrix_.rows.size() - 1;
    pardiso(pparam_,
            &max_factors,
            &factor_number,
            &matrix_type,
            &phase,
            &size,
            values,
            matrix_.rows.data(),
            matrix_.columns.data(),
            permutation,
            &rhs_count,
            iparm_,
            &message_level_,
            const_cast<double*>(rhs),
            solution,
            &error);
    if (VLOG_IS_ON(3)) {
      const std::array<MKL_INT, kPardisoParameterCount> current_iparm =
          CopyPardisoIparm(iparm_);
      if (phase != kPardisoRelease && has_iparm_snapshot_) {
        const std::vector<int> changed_indices =
            GetChangedPardisoIparmIndices(previous_iparm_, current_iparm);
        if (!changed_indices.empty()) {
          VLOG(3) << FormatPardisoIparmTable(
              absl::StrFormat("PARDISO iparm changes after %s (%d)",
                              PardisoPhaseToString(phase),
                              phase),
              changed_indices,
              current_iparm,
              &previous_iparm_);
        }
        previous_iparm_ = current_iparm;
      }
      VLOG(3) << absl::StrFormat(
          "PARDISO %s (%d) returned %s (%d) for %d x %d matrix, requested "
          "threads %d, active MKL threads %d.",
          PardisoPhaseToString(phase),
          phase,
          PardisoErrorToString(error),
          error,
          size,
          size,
          max_num_threads_,
          mkl_get_max_threads());
    }
    return error;
  }

  const OrderingType ordering_type_;
  const int max_num_threads_;
  const bool use_two_level_factorization_;
  MKL_INT iparm_[kPardisoParameterCount] = {};
  void* pparam_[kPardisoParameterCount] = {};
  std::array<MKL_INT, kPardisoParameterCount> previous_iparm_{};
  MKL_INT message_level_ = 0;
  bool initialized_ = false;
  bool has_iparm_snapshot_ = false;
  bool analyzed_ = false;
  PardisoMatrixState matrix_;
  std::vector<MKL_INT> permutation_;
};

bool ComputePardisoOrdering(const CompressedRowSparseMatrix& matrix,
                            const LinearSolverOrderingType ordering_type,
                            const int max_num_threads,
                            int* ordering,
                            std::string* message) {
  PardisoSolver solver(
      ordering_type == AMD ? OrderingType::AMD : OrderingType::NESDIS,
      max_num_threads,
      /*use_two_level_factorization=*/false);
  return solver.ComputeOrdering(matrix, ordering, message);
}

MklSparseCholesky::MklSparseCholesky(const OrderingType ordering_type,
                                     const int max_num_threads,
                                     const bool use_two_level_factorization)
    : solver_(std::make_unique<PardisoSolver>(
          ordering_type, max_num_threads, use_two_level_factorization)) {}

MklSparseCholesky::~MklSparseCholesky() = default;

std::unique_ptr<MklSparseCholesky> MklSparseCholesky::Create(
    const OrderingType ordering_type,
    const int max_num_threads,
    const bool use_two_level_factorization) {
  return std::unique_ptr<MklSparseCholesky>(new MklSparseCholesky(
      ordering_type, max_num_threads, use_two_level_factorization));
}

CompressedRowSparseMatrix::StorageType MklSparseCholesky::StorageType() const {
  return CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR;
}

LinearSolverTerminationType MklSparseCholesky::Factorize(
    CompressedRowSparseMatrix* lhs, std::string* message) {
  EventLogger event_logger("MklSparseCholesky::Factorize");
  if (!solver_->DefineStructure(*lhs, message)) {
    return LinearSolverTerminationType::FATAL_ERROR;
  }
  event_logger.AddEvent("Define structure");
  solver_->RefreshValues(*lhs);
  const LinearSolverTerminationType analysis_termination_type =
      solver_->AnalyzeStructure(message);
  if (analysis_termination_type != LinearSolverTerminationType::SUCCESS) {
    return analysis_termination_type;
  }
  event_logger.AddEvent("Analyze structure");
  const LinearSolverTerminationType termination_type =
      solver_->Factorize(message);
  event_logger.AddEvent("Factorize");
  return termination_type;
}

LinearSolverTerminationType MklSparseCholesky::Solve(const double* rhs,
                                                     double* solution,
                                                     std::string* message) {
  return solver_->Solve(rhs, solution, message);
}

}  // namespace ceres::internal

#endif  // CERES_NO_MKL
