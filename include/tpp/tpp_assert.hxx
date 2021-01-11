/*!*************************************************************************************************
 * \file    tpp_assert.hxx
 * \brief   This file provides the function \c tpp_assert.
 *
 * This is a wrapper file to provide a function that allows to use assertions that are similar to
 * those provided by cassert. That is, we define a macro \c tpp_assert that implements assert.
 * If a user wants to use assertions, it is recommended to use \c tpp_assert(\c Expr, \c Msg). The
 * use of the function \c __TPP_Assert is \b not recommended.
 *
 * Whether this functionality is active or not can be deduced via setting \c NDEBUG, when the code
 * is compiled. Using this functionality makes your program significantly slower. However, usage is
 * highly recommended for testing.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/

#pragma once  // Ensure that file is included only once in a single compilation.

#ifndef NDEBUG
#include <iostream>
#endif

namespace TPP
{
#ifndef NDEBUG

/*!*************************************************************************************************
 * \brief   The assertion to be used within tpp --- deactivate using -DNDEBUG compile flag.
 *
 * \param   Expr  C++ Expression that can be evaluated to \c true or \c false.
 * \param   Msg   Message that is to be displayed if \c Expr is evaluated to \c false.
 **************************************************************************************************/
#define tpp_assert(Expr, Msg) __TPP_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   This function is not (never) to be used.
 *
 * This function is \b not to be used in regular code. It only / solely is defined to allow the use
 * of function \c tpp_assert( \c Expr, \c Msg) which is implemented as a macro in file HyAssert.hxx.
 **************************************************************************************************/
constexpr void __TPP_Assert(const char* exp_str,
                            bool exp,
                            const char* file,
                            int line,
                            const char* msg)
{
  if (!exp)
  {
    std::cerr << "Assert failed:  " << msg << std::endl
              << "Expected:       " << exp_str << std::endl
              << "Source:         " << file << ", line " << line << std::endl;
    abort();
  }
}

#else  // alternative branch of ifndef NDEBUG
#define tpp_assert(Expr, Msg) ;
#endif  // end of ifndef NDEBUG

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------

}  // end of namespace TPP
