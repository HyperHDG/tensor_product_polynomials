PROJECT     	= TensorProductPolynomials
.PHONY:       	doxygen format

doxygen:
	cd doxygen; doxygen Doxyfile

format:
	clang-format -i include/tpp/*.hxx include/tpp/*/*.hxx
