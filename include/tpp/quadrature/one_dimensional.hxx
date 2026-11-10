#pragma once  // Ensure that file is included only once in a single compilation.

#include <tpp/compile_time_tricks.hxx>

#include <array>
#include <limits>

namespace TPP
{
namespace Quadrature
{
/*!*************************************************************************************************
 * \brief   Gauss--Legendre quadrature points and weights in one spatial dimension.
 *
 * \tparam  quad_deg  Desired degree of accuracy.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/
template <unsigned int quad_deg>
class GaussLegendre
{
 public:
  /*!***********************************************************************************************
   * \brief   Calculate the amount of quadrature points at compile time.
   *
   * Naive implementation to calculate the amount of needed quadrature points if a rule of accuracy
   * \c quad_deg is desired in one dimension.
   *
   * \retval  n_quad_points       Amount of needed quadrature points.
   ************************************************************************************************/
  static constexpr unsigned int n_points()
  {
    unsigned int amount = 1;
    for (; 2 * amount - 1 < quad_deg; ++amount)
      ;
    return amount;
  }

 private:
  /*!***********************************************************************************************
   * \brief   Transform array of quadrature points in interval \f$[-1,1]\f$ to \f$[0,1]\f$.
   ************************************************************************************************/
  template <typename return_t>
  static constexpr std::array<return_t, n_points()> transform_points(
    std::array<return_t, n_points()> points)
  {
    for (unsigned int index = 0; index < n_points(); ++index)
      points[index] = 0.5 * (points[index] + 1.);
    return points;
  }
  /*!***********************************************************************************************
   * \brief   Transform array of quadrature weights for interval \f$[-1,1]\f$ to \f$[0,1]\f$.
   ************************************************************************************************/
  template <typename return_t>
  static constexpr std::array<return_t, n_points()> transform_weights(
    std::array<return_t, n_points()> weights)
  {
    for (unsigned int index = 0; index < n_points(); ++index)
      weights[index] *= 0.5;
    return weights;
  }

 private:
  /*!***********************************************************************************************
   * \brief   Absolute value, usable in constant expressions for any floating point type.
   ************************************************************************************************/
  template <typename return_t>
  static constexpr return_t absolute(const return_t x)
  {
    return x < (return_t)0. ? -x : x;
  }
  /*!***********************************************************************************************
   * \brief   Cosine by its Taylor series, usable in constant expressions; only for the initial
   *          guesses of the points, so its accuracy is irrelevant.
   ************************************************************************************************/
  template <typename return_t>
  static constexpr return_t cosine(const return_t x)
  {
    return_t term = 1., sum = 1.;
    for (unsigned int k = 1; k < 100 && absolute(term) > std::numeric_limits<return_t>::epsilon();
         ++k)
    {
      term *= -x * x / (return_t)((2 * k - 1) * (2 * k));
      sum += term;
    }
    return sum;
  }
  /*!***********************************************************************************************
   * \brief   Legendre polynomial P_n on [-1,1] and its derivative at x, n = n_points().
   ************************************************************************************************/
  template <typename return_t>
  static constexpr std::array<return_t, 2> legendre(const return_t x)
  {
    return_t p_previous = 1., p = x;
    for (unsigned int k = 2; k <= n_points(); ++k)
    {
      const return_t p_next =
        ((return_t)(2 * k - 1) * x * p - (return_t)(k - 1) * p_previous) / (return_t)k;
      p_previous = p;
      p = p_next;
    }
    if constexpr (n_points() == 1)
      p_previous = 1.;
    return {p, (return_t)n_points() * (x * p - p_previous) / (x * x - (return_t)1.)};
  }
  /*!***********************************************************************************************
   * \brief   Points (ascending, in [-1,1]) and weights of the Gauss--Legendre rule, in return_t.
   *
   * The points are the roots of P_n, found by Newton's method from the classical initial guesses
   * -cos(pi (i + 3/4) / (n + 1/2)) to the precision of return_t, and the weights are
   * 2 / ((1 - x^2) P_n'(x)^2). All in return_t, so that the rule is as accurate as the type, and
   * usable in constant expressions.
   ************************************************************************************************/
  template <typename return_t>
  static constexpr std::array<std::array<return_t, n_points()>, 2> rule()
  {
    constexpr return_t pi = (return_t)3.14159265358979323846264338327950288419716939937510L;
    std::array<std::array<return_t, n_points()>, 2> rule;
    for (unsigned int i = 0; i < n_points(); ++i)
    {
      return_t x =
        -cosine<return_t>(pi * ((return_t)i + (return_t)0.75) / ((return_t)n_points() + 0.5));
      for (unsigned int iteration = 0; iteration < 100; ++iteration)
      {
        const std::array<return_t, 2> p = legendre<return_t>(x);
        const return_t step = p[0] / p[1];
        x -= step;
        if (absolute(step) <= std::numeric_limits<return_t>::epsilon())
          break;
      }
      const return_t derivative = legendre<return_t>(x)[1];
      rule[0][i] = x;
      rule[1][i] = (return_t)2. / (((return_t)1. - x * x) * derivative * derivative);
    }
    return rule;
  }

 public:
  /*!***********************************************************************************************
   * \brief   Gauss--Legendre quadrature points on one-dimensional unit interval.
   *
   * Returns the quadrature points of the quadrature rule with accuracy order \c quad_deg on the
   * one-dimensional unit interval \f$[0,1]\f$, in ascending order, computed in \c return_t to its
   * precision.
   *
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  quad_points         \c std::array containing the quadrature points.
   ************************************************************************************************/
  template <typename return_t = double>
  static constexpr std::array<return_t, n_points()> points()
  {
    return transform_points<return_t>(rule<return_t>()[0]);
  }
  /*!***********************************************************************************************
   * \brief   Gauss--Legendre quadrature weights on one-dimensional unit interval.
   *
   * Returns the quadrature weights of the quadrature rule with accuracy order \c quad_deg on the
   * one-dimensional unit interval \f$[0,1]\f$, in the order of the points, computed in
   * \c return_t to its precision.
   *
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  quad_weights        \c std::array containing the quadrature weights.
   ************************************************************************************************/
  template <typename return_t = double>
  static constexpr std::array<return_t, n_points()> weights()
  {
    return transform_weights<return_t>(rule<return_t>()[1]);
  }
};  // end of struct GaussLegendre

}  // end of namespace Quadrature

}  // end of namespace TPP
