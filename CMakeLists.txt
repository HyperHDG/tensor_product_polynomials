cmake_minimum_required(VERSION 3.16)

# Detect if we are running bundled into a larger project
if(NOT DEFINED PROJECT_NAME)
  set(SUBPROJECT OFF)
else()
  set(SUBPROJECT ON)
endif()

project(TensorProductPolynomials LANGUAGES CXX)

# Define an option to build documentation
option(BUILD_DOC "Whether to build a Doxygen documentation" ON)

# Create a target that transports all the necessary information to downstream projects
add_library(tpp INTERFACE)

# Add the correct include directories
target_include_directories(
  tpp
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Add the necessary C++ standard requirement
target_compile_features(
  tpp
  INTERFACE
  cxx_std_20
)

# Add tests if requested through the BUILD_TESTING option (default: ON)
include(CTest)
if(BUILD_TESTING AND NOT SUBPROJECT)
  enable_testing()
  add_subdirectory(tests)
endif()

# Build documentation if requested through the BUILD_DOC option (default: ON)
if(BUILD_DOC AND NOT SUBPROJECT)
  add_subdirectory(doxygen)
endif()

# Add a namespaced alias target for the export
add_library(HyperHDG::tpp ALIAS tpp)

# Install targets, headers and configuration
include(GNUInstallDirs)
install(
  TARGETS tpp
  EXPORT tpp-config
)

install(
  EXPORT tpp-config
  NAMESPACE HyperHDG::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tpp
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
