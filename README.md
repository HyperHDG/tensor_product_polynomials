# Welcome to Tensor Product Polynomials (TPP)!

It contains a small, C++ based library implementing the evaluation and numerical integration of
tensor-product type polynomials on tensor-product type volumes and (arbitrary dimensional) surfaces.

This library is heavily used by our main (yet unpublished) library "HyperHDG", which solves partial
differential equations on HyperGraphs using hybrid discontinuous Galerkin methods. Thus, TPP
provides all functions that deal with polynomials and are needed for HyperHDG.

Some of the design goals are to:
- exploit the tensor product structures of polynomials and domains to efficiently evaluate, derive,
and integrate.
- enable as much compile time computation as possible, and therefore reduce runtime overhead,
- provide a header-only library without configuration step which can be used easily in other
projects.


# How to use TPP?

## Compiler requirements

TPP is a fully header based library that uses the C++20 standard. It has been tested to work for all
compilers that are allowed for HyperHDG (cf. HyperHDG compilers). Thus, the first step to use the
library consists in obtaining one of these compilers.


## Install TPP

Since TPP is a header-only library, you can simply download it from GitHub and put the contents of
the `include` directory somewhere where the compiler finds it. There is no installation necessary.


## Obtain the documentation

The documentation can either be obtained via the command `make doxygen` after cloning and entering
the repository, or by having a look at the `namespace TPP` within the [online documentation](
https://hyperhdg.github.io/auto_pages/doxygen/namespaceTPP.html) of HyperHDG. Please note that the
online documentation might be outdated if HyperHDG does not use the latest version of TPP.
