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

# Usage: enable_sanitizer(REQUIRED_SANITIZERS) where REQUIRED_SANITIZERS should
# contain the list of sanitizers to enable by adding the corresponding compile
# and link options to the current directory.
#
# The specified sanitizers will be checked both for compatibility with the
# current compiler and with each other as some sanitizers are mutually
# exclusive.
function(enable_sanitizer)
  # According to the Clang documentation [1] the following sanitizers are
  # mutually exclusive.
  # [1]: https://clang.llvm.org/docs/UsersManual.html#controlling-code-generation
  set(INCOMPATIBLE_SANITIZERS address thread memory)
  # Set the recommended additional common compile flags for any sanitizer to
  # get the best possible output, e.g [2] but make them visible in the cache
  # so that the user can edit them if required.
  # [2]: https://clang.llvm.org/docs/AddressSanitizer.html#usage
  set(COMMON_SANITIZER_COMPILE_OPTIONS
    "-g -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Common compile flags enabled for any sanitizer")

  # Check that the specified list of sanitizers to enable does not include
  # multiple entries from the incompatible list.
  string(JOIN "|" INCOMPATIBLE_SANITIZERS_REGEX ${INCOMPATIBLE_SANITIZERS})
  set(REQUESTED_INCOMPATIBLE_SANITIZERS ${ARGN})
  list(REMOVE_DUPLICATES REQUESTED_INCOMPATIBLE_SANITIZERS)
  list(FILTER REQUESTED_INCOMPATIBLE_SANITIZERS INCLUDE REGEX
    "^(${INCOMPATIBLE_SANITIZERS_REGEX})$")
  list(LENGTH REQUESTED_INCOMPATIBLE_SANITIZERS
    REQUESTED_INCOMPATIBLE_SANITIZER_COUNT)
  if (REQUESTED_INCOMPATIBLE_SANITIZER_COUNT GREATER 1)
    include(PrettyPrintCMakeList)
    pretty_print_cmake_list(REQUESTED_SANITIZERS ${ARGN})
    pretty_print_cmake_list(
      PRETTY_INCOMPATIBLE_SANITIZERS ${INCOMPATIBLE_SANITIZERS})
    message(FATAL_ERROR "Found incompatible sanitizers in requested set: "
      "${REQUESTED_SANITIZERS}. The following sanitizers are mutually "
      "exclusive: ${PRETTY_INCOMPATIBLE_SANITIZERS}")
  endif()

  include(CheckCXXCompilerFlag)

  set(SANITIZER_FLAGS)
  foreach(REQUESTED_SANITIZER IN LISTS ARGN)
    set(SANITIZER_FLAG -fsanitize=${REQUESTED_SANITIZER})
    # The check caches its result. Each sanitizer therefore requires its own
    # result variable.
    string(MAKE_C_IDENTIFIER "HAVE_SANITIZER_${REQUESTED_SANITIZER}"
      HAVE_SANITIZER)
    # The sanitizer flag must be passed to the linker as well.
    set(CMAKE_REQUIRED_LINK_OPTIONS ${SANITIZER_FLAG})
    check_cxx_compiler_flag(${SANITIZER_FLAG} ${HAVE_SANITIZER})
    if (NOT ${HAVE_SANITIZER})
      message(FATAL_ERROR "Specified sanitizer: ${REQUESTED_SANITIZER} is not "
        "supported by the compiler.")
    endif()
    message(DEBUG "Enabling sanitizer: ${REQUESTED_SANITIZER}")
    list(APPEND SANITIZER_FLAGS ${SANITIZER_FLAG})
  endforeach()

  if (SANITIZER_FLAGS)
    separate_arguments(COMMON_SANITIZER_COMPILE_OPTIONS NATIVE_COMMAND
      "${COMMON_SANITIZER_COMPILE_OPTIONS}")
    # As per the Clang documentation, the sanitizer flags must be added to both
    # the compiler and linker flags.
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:C,CXX>:${SANITIZER_FLAGS};${COMMON_SANITIZER_COMPILE_OPTIONS}>")
    add_link_options(${SANITIZER_FLAGS})
  endif()
endfunction()
