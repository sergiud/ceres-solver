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

#include "ceres/mkl_pardiso_diagnostics.h"

#ifndef CERES_NO_MKL

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <vector>

#include "absl/strings/str_format.h"

namespace ceres::internal {

namespace {

constexpr MKL_INT kPardisoInputInconsistent = -1;
constexpr MKL_INT kPardisoNotEnoughMemory = -2;
constexpr MKL_INT kPardisoReorderingProblem = -3;
constexpr MKL_INT kPardisoNumericalFactorizationProblem = -4;
constexpr MKL_INT kPardisoInternalProblem = -5;
constexpr MKL_INT kPardisoNonsymmetricReorderingFailed = -6;
constexpr MKL_INT kPardisoSingularDiagonalMatrix = -7;
constexpr MKL_INT kPardisoIntegerOverflow = -8;
constexpr MKL_INT kPardisoNotEnoughOutOfCoreMemory = -9;
constexpr MKL_INT kPardisoOutOfCoreFileOpenError = -10;
constexpr MKL_INT kPardisoOutOfCoreFileAccessError = -11;

constexpr int kPardisoIparmColumnWidth = sizeof("iparm[]") - 1;

std::string_view PardisoIparmValueDescription(const int index,
                                              const MKL_INT value) {
  switch (index) {
    case kPardisoUseDefaultValuesParameter:
      return value == kPardisoDisabled ? "use PARDISO defaults"
                                       : "use supplied values";
    case kPardisoOrderingParameter:
      switch (value) {
        case kPardisoUseMinimumDegreeOrdering:
          return "minimum degree";
        case kPardisoMetisOrdering:
          return "METIS nested dissection";
        case kPardisoUseNestedDissectionOrdering:
          return "parallel nested dissection";
        default:
          return "unknown ordering";
      }
    case kPardisoPreconditionedCgParameter:
      return value == kPardisoDisabled ? "direct factorization and solve"
                                       : "iterative mode";
    case kPardisoIterativeRefinementParameter:
      return value == kPardisoDisabled
                 ? "automatic refinement on perturbed pivots"
                 : "maximum configured steps";
    case kPardisoUserPermutationParameter:
      switch (value) {
        case kPardisoIgnoreUserPermutation:
          return "ignore user permutation";
        case kPardisoUseUserPermutation:
          return "use supplied permutation";
        case kPardisoReturnPermutation:
          return "return computed permutation";
        default:
          return "unknown permutation mode";
      }
    case kPardisoParallelFactorizationParameter:
      switch (value) {
        case kPardisoClassicFactorization:
          return "classic factorization";
        case kPardisoTwoLevelFactorization:
          return "two-level factorization";
        default:
          return "unknown factorization mode";
      }
    case kPardisoParallelSolveParameter:
      switch (value) {
        case kPardisoDisabled:
          return "parallel by matrix or right-hand sides";
        case kPardisoSequentialSolve:
          return "sequential solve";
        case kPardisoParallelSolve:
          return "parallel matrix partitioning";
        default:
          return "unknown solve mode";
      }
    case kPardisoMatrixCheckerParameter:
      return value == kPardisoDisabled ? "disabled" : "enabled";
    case kPardisoIndexBaseParameter:
      switch (value) {
        case kPardisoOneBasedIndexing:
          return "one-based";
        case kPardisoZeroBasedIndexing:
          return "zero-based";
        default:
          return "unknown index base";
      }
    case kPardisoFactorNonzerosParameter:
    case kPardisoFactorOperationsParameter:
      return value < 0 ? "report enabled" : std::string_view{};
    default:
      return {};
  }
}

std::string FormatPardisoIparmValue(const int index, const MKL_INT value) {
  const std::string_view description =
      PardisoIparmValueDescription(index, value);
  if (description.empty()) {
    return absl::StrFormat("%d", value);
  }
  return absl::StrFormat("%s (%d)", description, value);
}

}  // namespace

std::string_view PardisoErrorToString(const MKL_INT error) {
  switch (error) {
    case PARDISO_NO_ERROR:
      return "success";
    case kPardisoInputInconsistent:
      return "input inconsistent";
    case kPardisoNotEnoughMemory:
      return "not enough memory";
    case kPardisoReorderingProblem:
      return "reordering problem";
    case kPardisoNumericalFactorizationProblem:
      return "zero pivot or numerical factorization problem";
    case kPardisoInternalProblem:
      return "internal problem";
    case kPardisoNonsymmetricReorderingFailed:
      return "reordering failed for a nonsymmetric matrix";
    case kPardisoSingularDiagonalMatrix:
      return "singular diagonal matrix";
    case kPardisoIntegerOverflow:
      return "32-bit integer overflow";
    case kPardisoNotEnoughOutOfCoreMemory:
      return "not enough memory for out-of-core execution";
    case kPardisoOutOfCoreFileOpenError:
      return "error opening out-of-core files";
    case kPardisoOutOfCoreFileAccessError:
      return "read or write error with out-of-core files";
    default:
      return "unknown";
  }
}

bool CheckPardisoStatus(const MKL_INT actual,
                        const char* operation,
                        std::string* message) {
  if (actual == PARDISO_NO_ERROR) {
    return true;
  }

  *message = absl::StrFormat("PARDISO %s returned %s (%d), expected %s (%d).",
                             operation,
                             PardisoErrorToString(actual),
                             actual,
                             PardisoErrorToString(PARDISO_NO_ERROR),
                             PARDISO_NO_ERROR);
  return false;
}

LinearSolverTerminationType PardisoErrorToTerminationType(const MKL_INT error) {
  switch (error) {
    case PARDISO_NO_ERROR:
      return LinearSolverTerminationType::SUCCESS;
    case kPardisoNumericalFactorizationProblem:
    case kPardisoSingularDiagonalMatrix:
      return LinearSolverTerminationType::FAILURE;
    default:
      return LinearSolverTerminationType::FATAL_ERROR;
  }
}

const std::vector<PardisoIparmEntry>& GetPardisoIparmEntries() {
  static const std::vector<PardisoIparmEntry> entries = {
      {kPardisoUseDefaultValuesParameter, "Defaults for unspecified entries"},
      {kPardisoOrderingParameter, "Fill-in reducing ordering"},
      {kPardisoPreconditionedCgParameter,
       "Iterative factorization or solution"},
      {kPardisoUserPermutationParameter, "Supplied permutation"},
      {kPardisoRefinementStepsPerformedParameter,
       "Iterative refinement steps performed"},
      {kPardisoIterativeRefinementParameter,
       "Maximum iterative refinement steps"},
      {kPardisoPerturbedPivotsParameter, "Pivots modified for stability"},
      {kPardisoPeakSymbolicMemoryParameter, "Peak symbolic memory in KiB"},
      {kPardisoPermanentSymbolicMemoryParameter,
       "Permanent symbolic memory in KiB"},
      {kPardisoFactorMemoryParameter, "Factor memory in KiB"},
      {kPardisoFactorNonzerosParameter, "Factor nonzeros"},
      {kPardisoFactorOperationsParameter,
       "Factorization operations in millions"},
      {kPardisoCgDiagnosticsParameter, "CG or CGS diagnostics"},
      {kPardisoParallelFactorizationParameter, "Parallel factorization"},
      {kPardisoParallelSolveParameter, "Parallel solve"},
      {kPardisoMatrixCheckerParameter, "Sparse matrix validation"},
      {kPardisoIndexBaseParameter, "CSR row and column index base"}};
  return entries;
}

std::vector<int> GetChangedPardisoIparmIndices(
    const std::array<MKL_INT, kPardisoParameterCount>& previous,
    const std::array<MKL_INT, kPardisoParameterCount>& current) {
  const std::vector<PardisoIparmEntry>& entries = GetPardisoIparmEntries();
  std::vector<int> changed_indices;
  changed_indices.reserve(entries.size());
  for (const PardisoIparmEntry& entry : entries) {
    if (previous[entry.index] != current[entry.index]) {
      changed_indices.push_back(entry.index);
    }
  }
  return changed_indices;
}

std::string_view PardisoPhaseToString(const MKL_INT phase) {
  switch (phase) {
    case kPardisoAnalysis:
      return "symbolic analysis";
    case kPardisoNumericalFactorization:
      return "numerical factorization";
    case kPardisoSolve:
      return "solve";
    case kPardisoRelease:
      return "release";
    default:
      return "unknown phase";
  }
}

std::string FormatPardisoIparmTable(
    const std::string_view title,
    const std::vector<int>& indices,
    const std::array<MKL_INT, kPardisoParameterCount>& current,
    const std::array<MKL_INT, kPardisoParameterCount>* previous) {
  const std::vector<PardisoIparmEntry>& entries = GetPardisoIparmEntries();
  std::string table = absl::StrFormat("\n%s\n", title);
  if (previous == nullptr) {
    absl::StrAppendFormat(&table,
                          "  %*s %-70s %s\n",
                          kPardisoIparmColumnWidth,
                          "iparm[]",
                          "Parameter",
                          "Value");
  } else {
    absl::StrAppendFormat(&table,
                          "  %*s %-70s %-30s %s\n",
                          kPardisoIparmColumnWidth,
                          "iparm[]",
                          "Parameter",
                          "Previous",
                          "Current");
  }
  for (const int index : indices) {
    const auto entry = std::find_if(entries.begin(),
                                    entries.end(),
                                    [index](const PardisoIparmEntry& value) {
                                      return value.index == index;
                                    });
    if (entry == entries.end()) {
      continue;
    }
    const std::string current_value =
        FormatPardisoIparmValue(index, current[index]);
    if (previous == nullptr) {
      absl::StrAppendFormat(&table,
                            "  %*d %-70s %s\n",
                            kPardisoIparmColumnWidth,
                            entry->index,
                            entry->name,
                            current_value);
    } else {
      const std::string previous_value =
          FormatPardisoIparmValue(index, (*previous)[index]);
      absl::StrAppendFormat(&table,
                            "  %*d %-70s %-30s %s\n",
                            kPardisoIparmColumnWidth,
                            entry->index,
                            entry->name,
                            previous_value,
                            current_value);
    }
  }
  return table;
}

}  // namespace ceres::internal

#endif  // CERES_NO_MKL
