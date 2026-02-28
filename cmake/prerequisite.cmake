# ##############################################################################
# Microsoft GSL (Guidelines Support Library)
# On HPC we cannot rely on FetchContent (no outbound https / broken git-remote-https),
# so we use a preinstalled package (e.g. Spack "cppgsl") or a vendored copy.
# ##############################################################################

include_guard(GLOBAL)

# Try to find headers like: <prefix>/include/gsl/gsl
find_path(MS_GSL_INCLUDE_DIR
  NAMES gsl/gsl
  HINTS
    ENV CPATH
    ENV CMAKE_PREFIX_PATH
  PATH_SUFFIXES include
)

if(MS_GSL_INCLUDE_DIR)
  message(STATUS "Found Microsoft GSL headers in: ${MS_GSL_INCLUDE_DIR}")
  add_library(GSL INTERFACE)
  target_include_directories(GSL INTERFACE "${MS_GSL_INCLUDE_DIR}")
else()
  message(FATAL_ERROR
    "Microsoft GSL headers not found (expected gsl/gsl).\n"
    "Fix options:\n"
    "  1) spack install cppgsl && spack load cppgsl\n"
    "  2) vendor microsoft/GSL into third_party/GSL and ensure third_party/GSL/include/gsl/gsl exists.\n")
endif()
