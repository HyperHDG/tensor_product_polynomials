PROJECT     	= TensorProductPolynomials
.PHONY:       	clean doxygen format test

clean:
	rm -rf build doxygen/html doxygen/latex

doxygen:
	cd doxygen; doxygen Doxyfile

format:
	clang-format -i include/tpp/*.hxx include/tpp/*/*.hxx tests/*.cxx

test:
	$(MAKE) clean
	mkdir -p build
	cd build; cmake ..
	cd build; make
	cd build; make test
