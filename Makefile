PROJECT     	= TensorProductPolynomials
.PHONY:       	clean doxygen format

clean:
	rm -rf doxygen/html doxygen/latex

doxygen:
	cd doxygen; doxygen Doxyfile

format:
	clang-format -i include/tpp/*.hxx include/tpp/*/*.hxx
