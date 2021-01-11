#pragma once  // Ensure that file is included only once in a single compilation.

#include <tpp/tpp_assert.hxx>

#include <limits>

namespace TPP
{
/*!*************************************************************************************************
 * \brief   Calculate the non-negative square-root of a non-negative number at compile time.
 *
 * \tparam  float_t The floating point type with respect to which the square root is evaluated.
 * \param   square  The number whose square-root is evaluated.
 * \retval  root    The square root of the given number.
 **************************************************************************************************/
template <typename float_t>
static constexpr float_t heron_root(const float_t square)
{
  tpp_assert(square >= 0., "Each square of a number must be non-negative!");

  if (square == 0.)
    return 0.;

  constexpr float_t corr_fac = 1. - 10. * std::numeric_limits<float_t>::epsilon();
  float_t upper_root = 0.5 * (square + 1.);
  float_t lower_root = square / upper_root;
  while (lower_root < upper_root * corr_fac)
  {
    upper_root = 0.5 * (lower_root + upper_root);
    lower_root = square / upper_root;
  }

  return 0.5 * (lower_root + upper_root);
}

}  // end of namespace TPP
