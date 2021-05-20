# Add a the test executable and CTest test
add_executable(one_dimensional one_dimensional.cxx)
target_link_libraries(one_dimensional PUBLIC tpp)
add_test(COMMAND one_dimensional)
