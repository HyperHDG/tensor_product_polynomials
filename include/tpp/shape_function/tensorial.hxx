#pragma once  // Ensure that file is included only once in a single compilation.

#include <tpp/hypercube.hxx>
#include <tpp/shape_function/one_dimensional.hxx>
#include <tpp/tpp_assert.hxx>

#include <array>

namespace TPP
{
namespace ShapeType
{
/*!*************************************************************************************************
 * \brief   Struct that handles the evaluation of tensorial shape functions.
 *
 * \tparam  shape_t   The typename of the one dimensional shape functions.
 * \tparam  dimT      Local dimension of the shape function's domain or the point.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/
template <typename shape_t, unsigned int dimT>
struct Tensorial
{
  /*!***********************************************************************************************
   * \brief   Make type of one-dimensional shape functions accessable to everyone.
   ************************************************************************************************/
  typedef shape_t shape_fun_1d;
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_fun() { return Hypercube<dimT>::pow(shape_fun_1d::n_fun()); }
  /*!***********************************************************************************************
   * \brief   Dimension of abscissas of the polynomials.
   ************************************************************************************************/
  static constexpr unsigned int dim() { return dimT; }
  /*!***********************************************************************************************
   * \brief   Maximum degree of all polynomials.
   ************************************************************************************************/
  static constexpr unsigned int degree() { return shape_t::degree(); }
  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c point.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  point_t       Type of abscissa value.
   * \param   index         Index of evaluated shape function.
   * \param   point         Abscissa / function argument.
   * \retval  fct_value     Evaluated value of shape function.
   ************************************************************************************************/
  template <typename return_t, typename point_t>
  static constexpr return_t fct_val(__attribute__((unused)) const unsigned int index,
                                    __attribute__((unused)) const point_t& point)
  {
    static_assert(point_t::size() == dim(), "Point needs to have correct dimension.");

    return_t value = 1.;
    if constexpr (dimT == 0)
      return (return_t)1.;
    else
    {
      std::array<unsigned int, std::max(dimT, 1U)> index_dim =
        Hypercube<dimT>::index_decompose(index, shape_fun_1d::n_fun());
      for (unsigned int dim = 0; dim < dimT; ++dim)
        value *= shape_fun_1d::template fct_val<return_t>(index_dim[dim], point[dim]);
    }
    return value;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  point_t       Type of abscissa value.
   * \param   index         Index of evaluated shape function.
   * \param   point         Abscissa / function argument.
   * \param   der_dim       Dimension with respect to which derivative is evaluated.
   * \retval  fct_value     Evaluated value of shape function's derivative.
   ************************************************************************************************/
  template <typename return_t, typename point_t>
  static constexpr return_t der_val(const unsigned int index,
                                    const point_t& point,
                                    const unsigned int der_dim)
  {
    static_assert(point_t::size() == dim(), "Point needs to have correct dimension.");
    tpp_assert(der_dim < dimT, "The derivative needs to be with respect to a valid dimension.");

    return_t value = 1.;
    std::array<unsigned int, std::max(dimT, 1U)> index_dim =
      Hypercube<dimT>::index_decompose(index, shape_fun_1d::n_fun());
    for (unsigned int dim = 0; dim < dimT; ++dim)
      if (der_dim == dim)
        value *= shape_fun_1d::template der_val<return_t>(index_dim[dim], point[dim]);
      else
        value *= shape_fun_1d::template fct_val<return_t>(index_dim[dim], point[dim]);

    return value;
  }
};  // end of struct Tensorial

}  // end of namespace ShapeType

}  // end of namespace TPP
