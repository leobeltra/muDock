# ##############################################################################
# ### Handle third party dependenies
# ##############################################################################

# Third-party libraries that must be manually installed
find_package(Boost CONFIG REQUIRED COMPONENTS program_options graph fiber)
find_package(OpenBabel3 REQUIRED)
find_package(MPI REQUIRED C)


# ##############################################################################
# Microsoft GSL (Guidelines Support Library)
# On HPC we cannot rely on FetchContent (no outbound https / broken git-remote-https),
# so we use a preinstalled package (e.g. Spack "cppgsl") or a vendored copy.
# ##############################################################################

# TODO remove GSL depndencies
include(FetchContent)
FetchContent_Declare(
  GSL
  GIT_REPOSITORY "https://github.com/microsoft/GSL"
  GIT_TAG "v4.0.0"
  GIT_SHALLOW ON)
FetchContent_MakeAvailable(GSL)


