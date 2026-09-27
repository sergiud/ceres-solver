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

#include "ceres/mkl_sparse_matrix.h"

#include <string>

#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/internal/config.h"
#include "gtest/gtest.h"

#ifndef CERES_NO_MKL
#include "mkl.h"
#endif  // CERES_NO_MKL

namespace ceres::internal {

#ifndef CERES_NO_MKL
TEST(MklSparseMatrix, MklThreadScopeDoesNotExceedMklThreadLimit) {
  const int limit = mkl_get_max_threads();
  {
    const MklThreadScope scope(limit + 1);
    EXPECT_EQ(mkl_get_max_threads(), limit);
  }
  EXPECT_EQ(mkl_get_max_threads(), limit);
}

TEST(MklSparseMatrix, CreateLeavesUnsortedInputUntouched) {
  constexpr int kNumRows = 1;
  constexpr int kNumCols = 3;
  constexpr int kNumNonzeros = 2;
  constexpr int kFirstColumn = 2;
  constexpr int kSecondColumn = 0;
  CompressedRowSparseMatrix input(kNumRows, kNumCols, kNumNonzeros);
  input.mutable_rows()[0] = 0;
  input.mutable_rows()[1] = kNumNonzeros;
  input.mutable_cols()[0] = kFirstColumn;
  input.mutable_cols()[1] = kSecondColumn;
  input.mutable_values()[0] = 1.0;
  input.mutable_values()[1] = 2.0;

  MklCsrMatrix matrix;
  std::string message;
  ASSERT_TRUE(matrix.Create(input, &message)) << message;
  EXPECT_NE(matrix.get(), nullptr);
  EXPECT_EQ(input.cols()[0], kFirstColumn);
  EXPECT_EQ(input.cols()[1], kSecondColumn);
  EXPECT_EQ(input.values()[0], 1.0);
  EXPECT_EQ(input.values()[1], 2.0);
}

#endif  // CERES_NO_MKL

}  // namespace ceres::internal
