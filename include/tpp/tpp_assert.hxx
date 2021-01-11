/*!*************************************************************************************************
 * \file    tpp_assert.hxx
 * \brief   This file provides the function \c tpp_assert.
 *
 * This is a wrapper file to provide a function that allows to use assertions that are similar to
 * those provided by cassert. That is, we define a macro \c hy_assert that implements assert.
 * If a user wants to use assertions, it is recommended to use \c hy_assert(\c Expr, \c Msg). The
 * use of the function \c __TPP_Assert is \b not recommended.
 *
 * Function \c tpp_assert takes two arguments. The first argument is evaluated to a \c boolean and
 * if this returns \c true, nothing is done. If the argument is \c false, the running program is
 * terminated and the second argument is displayed as part of an error message. Here, the second
 * argument is handled as a \c stringstream (without initial \c <<). Thus, for two integers a and b,
 * a function call might look like: tpp_assert( a == b , "Integers have not been the same, since a
 * turned out to be " << a << " and b was " << b << "." );
 *
 * Whether this functionality is active or not can be deduced via setting \c NDEBUG, when the code
 * is compiled. Using this functionality makes your program significantly slower. However, usage is
 * highly recommended for testing.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/

#pragma once  // Ensure that file is included only once in a single compilation.

namespace TPP
{
#ifndef NDEBUG

#include <iostream>

/*!*************************************************************************************************
 * \brief   The assertion to be used within HyperHDG --- deactivate using -DNDEBUG compile flag.
 *
 * \param   Expr  C++ Expression that can be evaluated to \c true or \c false.
 * \param   Msg   Message that is to be displayed if \c Expr is evaluated to \c false.
 **************************************************************************************************/
#define tpp_assert(Expr, Msg) __TPP_Assert(#Expr, Expr, __FILE__, __LINE__, __hy_assertion_text)

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   This function is not (never) to be used.
 *
 * This function is \b not to be used in regular code. It only / solely is defined to allow the use
 * of function \c hy_assert( \c Expr, \c Msg) which is implemented as a macro in file HyAssert.hxx.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
inline void __TPP_Assert(const char* exp_str, bool exp, const char* file, int line, const char* msg)
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
