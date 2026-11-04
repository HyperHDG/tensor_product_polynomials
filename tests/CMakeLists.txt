# Add a the test executable and CTest test
add_executable(one_dimensional one_dimensional.cxx)
target_link_libraries(one_dimensional PUBLIC tpp)
add_test(COMMAND one_dimensional)

# Gauss--Legendre rules in double, long double and std::float128_t (where available)
add_executable(quadrature quadrature.cxx)
target_link_libraries(quadrature PUBLIC tpp)
target_compile_features(quadrature PUBLIC cxx_std_23)
add_test(NAME quadrature COMMAND quadrature)
