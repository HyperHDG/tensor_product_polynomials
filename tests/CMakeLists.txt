# Add a the test executable and CTest test
add_executable(one_dimensional one_dimensional.cxx)
target_link_libraries(one_dimensional PUBLIC tpp)
add_test(NAME one_dimensional COMMAND one_dimensional)

# Gauss--Legendre rules, one executable per floating-point type (built in parallel):
# double, long double and std::float128_t where the compiler has it (e.g. GCC >= 13 in C++23).
include(CheckCXXSourceCompiles)
set(CMAKE_CXX_STANDARD 23)
check_cxx_source_compiles("#include <stdfloat>
int main() { std::float128_t x = 1; return x > 0 ? 0 : 1; }" TPP_HAVE_FLOAT128)
unset(CMAKE_CXX_STANDARD)
set(quadrature_types "double=double" "long_double=long double")
if(TPP_HAVE_FLOAT128)
  list(APPEND quadrature_types "float128=std::float128_t")
endif()
foreach(entry IN LISTS quadrature_types)
  string(REGEX MATCH "^([^=]*)=(.*)$" match "${entry}")
  set(name quadrature_${CMAKE_MATCH_1})
  add_executable(${name} quadrature.cxx)
  target_link_libraries(${name} PUBLIC tpp)
  target_compile_features(${name} PUBLIC cxx_std_23)
  target_compile_definitions(${name} PRIVATE "REAL=${CMAKE_MATCH_2}")
  add_test(NAME ${name} COMMAND ${name})
endforeach()
