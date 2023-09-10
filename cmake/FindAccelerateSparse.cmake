# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2023 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: alexs.mac@gmail.com (Alex Stewart)
#

# FindAccelerateSparse.cmake - Find the sparse solvers in Apple's Accelerate
#                              framework, introduced in Xcode 9.0 (2017).
#                              Note that this is distinct from the Accelerate
#                              framework on its own, which existed in previous
#                              versions but without the sparse solvers.
#
# This module defines the following variables which should be referenced
# by the caller to use the library.
#
# AccelerateSparse_FOUND: TRUE iff an Accelerate framework including the sparse
#                         solvers, and all dependencies, has been found.
# AccelerateSparse_INCLUDE_DIRS: Include directories for Accelerate framework.
# AccelerateSparse_LIBRARIES: Libraries for Accelerate framework and all
#                             dependencies.
#
# The following variables are also defined by this module, but in line with
# CMake recommended FindPackage() module style should NOT be referenced directly
# by callers (use the plural variables detailed above instead).  These variables
# do however affect the behaviour of the module via FIND_[PATH/LIBRARY]() which
# are NOT re-called (i.e. search for library is not repeated) if these variables
# are set with valid values _in the CMake cache_. This means that if these
# variables are set directly in the cache, either by the user in the CMake GUI,
# or by the user passing -DVAR=VALUE directives to CMake when called (which
# explicitly defines a cache variable), then they will be used verbatim,
# bypassing the HINTS variables and other hard-coded search locations.
#
# AccelerateSparse_INCLUDE_DIR: Include directory for Accelerate framework, not
#                               including the include directory of any
#                               dependencies.
# AccelerateSparse_LIBRARY: Accelerate framework, not including the libraries of
#                           any dependencies.

find_path (AccelerateSparse_INCLUDE_DIR
  NAMES Accelerate.h
  DOC "AccelerateSparse include directory"
)

find_library (AccelerateSparse_LIBRARY
  NAMES Accelerate
  DOC "AccelerateSparse library"
)

mark_as_advanced (AccelerateSparse_INCLUDE_DIR AccelerateSparse_LIBRARY)

# Determine if the Accelerate framework detected includes the sparse solvers.
include (CheckCXXSourceCompiles)
include (CMakePushCheckState)
include (FindPackageHandleStandardArgs)

if (AccelerateSparse_INCLUDE_DIR AND AccelerateSparse_LIBRARY)
  cmake_push_check_state (RESET)
  set (CMAKE_REQUIRED_INCLUDES ${AccelerateSparse_INCLUDE_DIR})
  set (CMAKE_REQUIRED_LIBRARIES ${AccelerateSparse_LIBRARY})

  check_cxx_source_compiles (
    "#include <Accelerate.h>
     int main() {
       SparseMatrix_Double A;
       SparseFactor(SparseFactorizationCholesky, A);
       return 0;
     }"
     ACCELERATE_FRAMEWORK_HAS_SPARSE_SOLVER)
  cmake_pop_check_state ()

  if (NOT ACCELERATE_FRAMEWORK_HAS_SPARSE_SOLVER)
    set (_MORE_FPHSA_ARGS FAIL_MESSAGE
      "Could not find AccelerateSparse because the Accelerate framework ${AccelerateSparse_LIBRARY} does not include the sparse solvers")
  endif (NOT ACCELERATE_FRAMEWORK_HAS_SPARSE_SOLVER)
endif (AccelerateSparse_INCLUDE_DIR AND AccelerateSparse_LIBRARY)

# Handle REQUIRED / QUIET optional arguments and version.
# Report variables that should be set by the user not the legacy convience variables.
find_package_handle_standard_args (AccelerateSparse
  REQUIRED_VARS AccelerateSparse_INCLUDE_DIR AccelerateSparse_LIBRARY
  ${_MORE_FPHSA_ARGS}
)

unset (_MORE_FPHSA_ARGS)

if (AccelerateSparse_FOUND)
  set(AccelerateSparse_INCLUDE_DIRS ${AccelerateSparse_INCLUDE_DIR})
  set(AccelerateSparse_LIBRARIES ${AccelerateSparse_LIBRARY})
endif (AccelerateSparse_FOUND)
