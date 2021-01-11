#pragma once  // Ensure that file is included only once in a single compilation.

#include <tpp/compile_time_tricks.hxx>

#include <array>

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

 public:
  /*!***********************************************************************************************
   * \brief   Gauss--Legendre quadrature points on one-dimensional unit interval.
   *
   * Returns the quadrature points of the quadrature rule with accuracy order \c quad_deg on the
   * one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  quad_points         \c std::array containing the quadrature points.
   ************************************************************************************************/
  template <typename return_t = double>
  static constexpr std::array<return_t, n_points()> points()
  {
    static_assert(n_points() < 10, "Amount of points needs to be smaller than 10!");

    if constexpr (n_points() == 1)
      return transform_points(std::array<return_t, n_points()>({static_cast<return_t>(0.)}));
    if constexpr (n_points() == 2)
      return transform_points(std::array<return_t, n_points()>(
        {-heron_root<return_t>(1. / 3.), heron_root<return_t>(1. / 3.)}));
    if constexpr (n_points() == 3)
      return transform_points(
        std::array<return_t, n_points()>({-heron_root<return_t>(3. / 5.), static_cast<return_t>(0.),
                                          heron_root<return_t>(3. / 5.)}));
    if constexpr (n_points() == 4)
      return transform_points(std::array<return_t, n_points()>(
        {-heron_root<return_t>(3. / 7. + 2. / 7. * heron_root<return_t>(6. / 5.)),
         -heron_root<return_t>(3. / 7. - 2. / 7. * heron_root<return_t>(6. / 5.)),
         heron_root<return_t>(3. / 7. - 2. / 7. * heron_root<return_t>(6. / 5.)),
         heron_root<return_t>(3. / 7. + 2. / 7. * heron_root<return_t>(6. / 5.))}));
    if constexpr (n_points() == 5)
      return transform_points(std::array<return_t, n_points()>(
        {static_cast<return_t>(-heron_root<return_t>(5. + 2. * heron_root<return_t>(10. / 7.)) /
                               3.),
         static_cast<return_t>(-heron_root<return_t>(5. - 2. * heron_root<return_t>(10. / 7.)) /
                               3.),
         static_cast<return_t>(0.),
         static_cast<return_t>(heron_root<return_t>(5. - 2. * heron_root<return_t>(10. / 7.)) / 3.),
         static_cast<return_t>(heron_root<return_t>(5. + 2. * heron_root<return_t>(10. / 7.)) /
                               3.)}));
    if constexpr (n_points() == 6)
      return transform_points(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.6612093864662645), static_cast<return_t>(-0.6612093864662645),
         static_cast<return_t>(-0.2386191860831969), static_cast<return_t>(0.2386191860831969),
         static_cast<return_t>(-0.9324695142031521), static_cast<return_t>(0.9324695142031521)}));
    if constexpr (n_points() == 7)
      return transform_points(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.0000000000000000), static_cast<return_t>(0.4058451513773972),
         static_cast<return_t>(-0.4058451513773972), static_cast<return_t>(-0.7415311855993945),
         static_cast<return_t>(0.7415311855993945), static_cast<return_t>(-0.9491079123427585),
         static_cast<return_t>(0.9491079123427585)}));
    if constexpr (n_points() == 8)
      return transform_points(std::array<return_t, n_points()>(
        {static_cast<return_t>(-0.1834346424956498), static_cast<return_t>(0.1834346424956498),
         static_cast<return_t>(-0.5255324099163290), static_cast<return_t>(0.5255324099163290),
         static_cast<return_t>(-0.7966664774136267), static_cast<return_t>(0.7966664774136267),
         static_cast<return_t>(-0.9602898564975363), static_cast<return_t>(0.9602898564975363)}));
    if constexpr (n_points() == 9)
      return transform_points(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.0000000000000000), static_cast<return_t>(-0.8360311073266358),
         static_cast<return_t>(0.8360311073266358), static_cast<return_t>(-0.9681602395076261),
         static_cast<return_t>(0.9681602395076261), static_cast<return_t>(-0.3242534234038089),
         static_cast<return_t>(0.3123470770400029), static_cast<return_t>(0.2606106964029354),
         static_cast<return_t>(0.2606106964029354)}));
  }
  /*!***********************************************************************************************
   * \brief   Gauss--Legendre quadrature weights on one-dimensional unit interval.
   *
   * Returns the quadrature weights of the quadrature rule with accuracy order \c quad_deg on the
   * one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  quad_weights        \c std::array containing the quadrature weights.
   ************************************************************************************************/
  template <typename return_t = double>
  static constexpr std::array<return_t, n_points()> weights()
  {
    static_assert(n_points() < 10, "Amount of points needs to be smaller than 10!");

    if constexpr (n_points() == 1)
      return transform_weights(std::array<return_t, n_points()>({static_cast<return_t>(2.)}));
    if constexpr (n_points() == 2)
      return transform_weights(
        std::array<return_t, n_points()>({static_cast<return_t>(1.), static_cast<return_t>(1.)}));
    if constexpr (n_points() == 3)
      return transform_weights(std::array<return_t, n_points()>({static_cast<return_t>(5. / 9.),
                                                                 static_cast<return_t>(8. / 9.),
                                                                 static_cast<return_t>(5. / 9.)}));
    if constexpr (n_points() == 4)
      return transform_weights(std::array<return_t, n_points()>(
        {static_cast<return_t>(1. / 36. * (18. - heron_root<return_t>(30.))),
         static_cast<return_t>(1. / 36. * (18. + heron_root<return_t>(30.))),
         static_cast<return_t>(1. / 36. * (18. + heron_root<return_t>(30.))),
         static_cast<return_t>(1. / 36. * (18. - heron_root<return_t>(30.)))}));
    if constexpr (n_points() == 5)
      return transform_weights(std::array<return_t, n_points()>(
        {static_cast<return_t>(1. / 900. * (322. - 13. * heron_root<return_t>(70.))),
         static_cast<return_t>(1. / 900. * (322. + 13. * heron_root<return_t>(70.))),
         static_cast<return_t>(1. / 900. * (322. + 190.)),
         static_cast<return_t>(1. / 900. * (322. + 13. * heron_root<return_t>(70.))),
         static_cast<return_t>(1. / 900. * (322. - 13. * heron_root<return_t>(70.)))}));
    if constexpr (n_points() == 6)
      return transform_weights(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.3607615730481386), static_cast<return_t>(0.3607615730481386),
         static_cast<return_t>(0.4679139345726910), static_cast<return_t>(0.4679139345726910),
         static_cast<return_t>(0.1713244923791704), static_cast<return_t>(0.1713244923791700)}));
    if constexpr (n_points() == 7)
      return transform_weights(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.4179591836734694), static_cast<return_t>(0.3818300505051189),
         static_cast<return_t>(0.3818300505051189), static_cast<return_t>(0.2797053914892766),
         static_cast<return_t>(0.2797053914892766), static_cast<return_t>(0.1294849661688697),
         static_cast<return_t>(0.1294849661688697)}));
    if constexpr (n_points() == 8)
      return transform_weights(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.3626837833783620), static_cast<return_t>(0.3626837833783620),
         static_cast<return_t>(0.3137066458778873), static_cast<return_t>(0.3137066458778873),
         static_cast<return_t>(0.2223810344533745), static_cast<return_t>(0.2223810344533745),
         static_cast<return_t>(0.1012285362903763), static_cast<return_t>(0.1012285362903763)}));
    if constexpr (n_points() == 9)
      return transform_weights(std::array<return_t, n_points()>(
        {static_cast<return_t>(0.3302393550012598), static_cast<return_t>(0.1806481606948574),
         static_cast<return_t>(0.1806481606948574), static_cast<return_t>(0.0812743883615744),
         static_cast<return_t>(0.0812743883615744), static_cast<return_t>(0.3123470770400029),
         static_cast<return_t>(0.3123470770400029), static_cast<return_t>(0.2606106964029354),
         static_cast<return_t>(0.2606106964029354)}));
  }
};  // end of struct GaussLegendre

}  // end of namespace Quadrature

}  // end of namespace TPP
