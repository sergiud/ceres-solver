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

#include <string>

#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/internal/config.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

#ifndef CERES_NO_MKL
namespace {

CompressedRowSparseMatrix CreateUpperPositiveDefiniteMatrix() {
  CompressedRowSparseMatrix matrix(2, 2, 3);
  matrix.mutable_rows()[0] = 0;
  matrix.mutable_rows()[1] = 2;
  matrix.mutable_rows()[2] = 3;
  matrix.mutable_cols()[0] = 0;
  matrix.mutable_cols()[1] = 1;
  matrix.mutable_cols()[2] = 1;
  matrix.mutable_values()[0] = 4.0;
  matrix.mutable_values()[1] = 1.0;
  matrix.mutable_values()[2] = 3.0;
  matrix.set_storage_type(
      CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
  return matrix;
}

}  // namespace

TEST(MklPardiso, SolvesWithTwoLevelFactorization) {
  auto matrix = CreateUpperPositiveDefiniteMatrix();
  auto solver = MklSparseCholesky::Create(OrderingType::AMD, 2, true);
  std::string message;
  ASSERT_EQ(solver->Factorize(&matrix, &message),
            LinearSolverTerminationType::SUCCESS)
      << message;

  const double rhs[] = {1.0, 2.0};
  double solution[] = {0.0, 0.0};
  ASSERT_EQ(solver->Solve(rhs, solution, &message),
            LinearSolverTerminationType::SUCCESS)
      << message;
  EXPECT_NEAR(solution[0], 1.0 / 11.0, 1e-12);
  EXPECT_NEAR(solution[1], 7.0 / 11.0, 1e-12);
}

TEST(MklPardiso, RejectsChangedNonzeroCountAfterAnalysis) {
  auto matrix = CreateUpperPositiveDefiniteMatrix();
  auto solver = MklSparseCholesky::Create(OrderingType::AMD, 1, false);
  std::string message;
  ASSERT_EQ(solver->Factorize(&matrix, &message),
            LinearSolverTerminationType::SUCCESS)
      << message;

  CompressedRowSparseMatrix diagonal(2, 2, 2);
  diagonal.mutable_rows()[0] = 0;
  diagonal.mutable_rows()[1] = 1;
  diagonal.mutable_rows()[2] = 2;
  diagonal.mutable_cols()[0] = 0;
  diagonal.mutable_cols()[1] = 1;
  diagonal.mutable_values()[0] = 4.0;
  diagonal.mutable_values()[1] = 3.0;
  diagonal.set_storage_type(
      CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
  EXPECT_EQ(solver->Factorize(&diagonal, &message),
            LinearSolverTerminationType::FATAL_ERROR);
  EXPECT_THAT(message, ::testing::HasSubstr("structure changed"));
}

TEST(MklPardiso, ReportsIndefiniteMatrixAsNumericalFailure) {
  auto matrix = CreateUpperPositiveDefiniteMatrix();
  matrix.mutable_values()[0] = -4.0;

  auto solver = MklSparseCholesky::Create(OrderingType::AMD, 1, false);
  std::string message;
  EXPECT_EQ(solver->Factorize(&matrix, &message),
            LinearSolverTerminationType::FAILURE);
}

#endif  // CERES_NO_MKL

}  // namespace ceres::internal
