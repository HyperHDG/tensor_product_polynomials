#pragma once  // Ensure that file is included only once in a single compilation.

#include <cmath>

#include <tpp/compile_time_tricks.hxx>
#include <tpp/tpp_assert.hxx>

namespace TPP
{
namespace ShapeType
{
/*!*************************************************************************************************
 * \brief   Struct that handles the evaluation of one-dimensional Legendre polynomials.
 *
 * We use Clenshaw's algorithm (cf. https://en.wikipedia.org/wiki/Clenshaw_algorithm; date: Jan. 09,
 * 2021) to evaluate the three-term recusion formula defining Legendre polynomials as defined on
 * https://en.wikipedia.org/wiki/Legendre_polynomials; date Jan. 09, 2021).
 *
 * The shape functions in this struct, however, are defined with respect to the unit interval [0,1]
 * and normalized to be L^2 orthonormal (not just orthogonal).
 *
 * \tparam  poly_deg  The maximum degree of the polynomial.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/
template <unsigned int poly_deg>
struct Legendre
{
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_fun() { return degree() + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of abscissas of the polynomials.
   ************************************************************************************************/
  static constexpr unsigned int dim() { return 1U; }
  /*!***********************************************************************************************
   * \brief   Maximum degree of all polynomials.
   ************************************************************************************************/
  static constexpr unsigned int degree() { return poly_deg; }

  /*!***********************************************************************************************
   * \brief   Evaluate variable alpha of the Clenshaw algorithm.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t clenshaw_alpha(const unsigned int k, const input_t& x)
  {
    return (((return_t)(2 * k + 1)) / ((return_t)(k + 1))) * ((return_t)x);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate variable beta of the Clenshaw algorithm.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t clenshaw_beta(const unsigned int k, const input_t&)
  {
    return -((return_t)k) / ((return_t)(k + 1));
  }
  /*!***********************************************************************************************
   * \brief   Evaluate variable b of the Clenshaw algorithm.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t clenshaw_b(const unsigned int k, const unsigned int n, const input_t& x)
  {
    if (k > n)
      return (return_t)0.;
    else if (k == n)
      return (return_t)1.;
    else
      return clenshaw_alpha<return_t>(k, x) * clenshaw_b<return_t>(k + 1, n, x) +
             clenshaw_beta<return_t>(k + 1, x) * clenshaw_b<return_t>(k + 2, n, x);
  }

  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t fct_val(const unsigned int index,
                                    const input_t& x_val,
                                    const bool normalized = true)
  {
    // Transform x value from reference interval [0,1] -> [-1,1].
    const return_t x = 2. * (return_t)x_val - 1.;

    // Evaluate the value of the Legendre polynomial at x.
    return_t legendre_val;
    switch (index)
    {
      case 0:
        legendre_val = 1.;
        break;
      case 1:
        legendre_val = x;
        break;
      default:
        legendre_val =
          x * clenshaw_b<return_t>(1, index, x) - 0.5 * clenshaw_b<return_t>(2, index, x);
    }

    // Return L^2 normalized value.
    if (normalized)
      legendre_val *= TPP::heron_root((return_t)(2. * index + 1.));

    return legendre_val;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate value of the derivative of orthonormal shape function.
   *
   * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape
   * function on the reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function's derivative.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t der_val(const unsigned int index,
                                    const input_t& x_val,
                                    const bool normalized = true)
  {
    if (index == 0)
      return (return_t)0.;

    // Transform x value from reference interval [0,1] -> [-1,1].
    const return_t x = 2. * (return_t)x_val - 1.;

    // Non-recusrive version of:
    // return_t legendre_val = ((return_t)index) * fct_val<return_t>(index - 1, x_val, false) +
    //                         x * der_val<return_t>(index - 1, x_val, false, false);
    return_t legendre_val = 0., x_pot = 1.;
    for (unsigned int p = index; p > 0; --p)
    {
      legendre_val += ((return_t)p) * x_pot * fct_val<return_t>(p - 1, x_val, false);
      x_pot *= x;
    }

    // Return L^2 normalized value.
    if (normalized)
      legendre_val *= TPP::heron_root((return_t)(2. * index + 1.));

    return 2. * legendre_val;
  }
};  // end of struct Legendre

/*!*************************************************************************************************
 * \brief   Struct that handles the evaluation of one-dimensional Lobatto polynomials.
 *
 * The \f$n\f$-th Lobatto polynomial \f$L_n\f$ is defined to be the definite integral of the
 * \f$n-1\f$-st Legendre polynomial \f$P_n$, i.e.,  \f$L_n(x) = \int_0^x P_{n-1}(t) \; dt\f$. Thus,
 * the Lobatto polynomials turn out to be \f$H^1\f$ orthonormal. They can be calculated according to
 * \f$L_n(x) = \frac{1}{n} (x P_{n-1}(x) - P_{n-2}(x))\f$, which can be deduced from the derivation
 * rules of the Legendre polynomials. The zeroth polynomial is constant.
 *
 * \tparam  poly_deg  The maximum degree of the polynomial.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/
template <unsigned int poly_deg>
struct Lobatto
{
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_fun() { return degree() + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of abscissas of the polynomials.
   ************************************************************************************************/
  static constexpr unsigned int dim() { return 1U; }
  /*!***********************************************************************************************
   * \brief   Maximum degree of all polynomials.
   ************************************************************************************************/
  static constexpr unsigned int degree() { return poly_deg; }

  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t fct_val(const unsigned int index,
                                    const input_t& x_val,
                                    const bool normalized = true)
  {
    if (index == 0)
      return 1.;
    else if (index == 1)
      return x_val;

    // Transform x value from reference interval [0,1] -> [-1,1].
    const return_t x = 2. * (return_t)x_val - 1.;

    // Evaluate the value of the Lobatto polynomial at x.
    return_t lobatto_val =
      (x * Legendre<poly_deg>::template fct_val<return_t>(index - 1, x_val, false) -
       Legendre<poly_deg>::template fct_val<return_t>(index - 2, x_val, false)) /
      (return_t)(index);

    // Return L^2 normalized value.
    if (normalized)
      lobatto_val *= TPP::heron_root((return_t)(2. * index - 1.));

    return 0.5 * lobatto_val;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate value of the derivative of orthonormal shape function.
   *
   * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape
   * function on the reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function's derivative.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t der_val(const unsigned int index,
                                    const input_t& x_val,
                                    const bool normalized = true)
  {
    if (index == 0)
      return 0.;

    return Legendre<poly_deg>::template fct_val<return_t>(index - 1, x_val, normalized);
  }
};  // end of struct Lobatto

/*!*************************************************************************************************
 * \brief   Struct that handles the evaluation of one-dimensional trigonometric functions.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2024.
 **************************************************************************************************/
template <unsigned int fourier_deg>
struct Foruier
{
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_fun() { return degree() + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of abscissas of the polynomials.
   ************************************************************************************************/
  static constexpr unsigned int dim() { return 1U; }
  /*!***********************************************************************************************
   * \brief   Maximum degree of all polynomials.
   ************************************************************************************************/
  static constexpr unsigned int degree() { return fourier_deg; }

  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t fct_val(const unsigned int index, const input_t& x_val)
  {
    switch (index)
    {
      case 0:
        return 1.;
      case index % 2 == 1:
        return cos(.5 * M_PI * x_val * (double)(index+1));
      default:
        return sin(.5 * M_PI * x_val * (double)(index));
    }
  }
  /*!***********************************************************************************************
   * \brief   Evaluate value of the derivative of orthonormal shape function.
   *
   * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape
   * function on the reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function's derivative.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t der_val(const unsigned int index, const input_t& x_val)
  {
    switch (index)
    {
      case 0:
        return 0.;
      case index % 2 == 1:
        return -0.5 * M_PI * (double)(index+1) * sin(.5 * M_PI * (double)(index+1) * x_val);
      default:
        return 0.5 * M_PI * (double)(index + 1) * cos(.5 * M_PI * (double)(i + 1) * x_val);
    }
  }
};  // end of struct Fourier

}  // end of namespace ShapeType

}  // end of namespace TPP
