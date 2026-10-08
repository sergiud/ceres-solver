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
# Authors: keir@google.com (Keir Mierle)
#          alexs.mac@gmail.com (Alex Stewart)

# Set up the git hook to make Gerrit Change-Id: lines in commit messages.
function(ADD_GERRIT_COMMIT_HOOK SOURCE_DIR BINARY_DIR)
  if (NOT (EXISTS ${SOURCE_DIR} AND IS_DIRECTORY ${SOURCE_DIR}))
    message(FATAL_ERROR "Specified SOURCE_DIR: ${SOURCE_DIR} does not exist, "
      "or is not a directory, cannot add Gerrit commit hook.")
  endif()
  if (NOT (EXISTS ${BINARY_DIR} AND IS_DIRECTORY ${BINARY_DIR}))
    message(FATAL_ERROR "Specified BINARY_DIR: ${BINARY_DIR} does not exist, "
      "or is not a directory, cannot add Gerrit commit hook.")
  endif()
  # Only a Git checkout of Ceres itself receives the hook. A copy of Ceres
  # inside the source tree of another repository must not modify its hooks.
  if (NOT EXISTS ${SOURCE_DIR}/.git)
    return()
  endif()

  find_package(Git QUIET)
  if (NOT Git_FOUND)
    return()
  endif()

  # Let Git locate the hooks directory, which also accounts for submodules,
  # worktrees and a custom core.hooksPath.
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --git-path hooks
    WORKING_DIRECTORY ${SOURCE_DIR}
    RESULT_VARIABLE GIT_RESULT
    OUTPUT_VARIABLE GIT_HOOKS_DIRECTORY
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (NOT GIT_RESULT EQUAL 0)
    return()
  endif()
  get_filename_component(GIT_HOOKS_DIRECTORY "${GIT_HOOKS_DIRECTORY}" ABSOLUTE
    BASE_DIR ${SOURCE_DIR})

  if (NOT EXISTS ${GIT_HOOKS_DIRECTORY}/commit-msg)
    # Hook installation is an internal setup detail.
    message(DEBUG "Adding commit hook for Gerrit to: ${GIT_HOOKS_DIRECTORY}")
    # Download the hook only if it is not already present.
    set(COMMIT_HOOK_URL
      https://ceres-solver-review.googlesource.com/tools/hooks/commit-msg)
    file(DOWNLOAD ${COMMIT_HOOK_URL} ${BINARY_DIR}/commit-msg
      STATUS COMMIT_HOOK_DOWNLOAD_STATUS)
    list(GET COMMIT_HOOK_DOWNLOAD_STATUS 0 COMMIT_HOOK_DOWNLOAD_ERROR)

    # An incomplete hook must not be installed since its presence prevents
    # the download from being retried.
    if (COMMIT_HOOK_DOWNLOAD_ERROR EQUAL 0)
      # Make the downloaded file executable, since it is not by default.
      file(COPY ${BINARY_DIR}/commit-msg
        DESTINATION ${GIT_HOOKS_DIRECTORY}/
        FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)
    else()
      list(GET COMMIT_HOOK_DOWNLOAD_STATUS 1 COMMIT_HOOK_DOWNLOAD_MESSAGE)
      message(WARNING "Downloading the Gerrit commit hook from "
        "${COMMIT_HOOK_URL} failed with status "
        "${COMMIT_HOOK_DOWNLOAD_ERROR} (${COMMIT_HOOK_DOWNLOAD_MESSAGE}), "
        "expected 0. The hook was not installed.")
    endif()
  endif()
endfunction()
