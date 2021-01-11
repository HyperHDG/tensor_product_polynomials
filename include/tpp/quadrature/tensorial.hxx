#pragma once  // Ensure that file is included only once in a single compilation.

#include <tpp/quadrature/one_dimensional.hxx>
#include <tpp/tpp_assert.hxx>

#include <array>
#include <cmath>
#include <numeric>

namespace TPP
{
namespace Quadrature
{
/*!*************************************************************************************************
 * \brief   General integrator class on tensorial hypergraphs.
 *
 * \tparam  quadrature_t      The one-dimensional quadrature rule applied.
 * \param   shape_t           Type of shape functions.
 * \tparam  return_t          Floating type specification. Default is double.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/
template <typename quadrature_t, typename shape_t, typename return_t = double>
class Tensorial
{
 public:
  /*!***********************************************************************************************
   * \brief   Make type of shape functions publicly available.
   ************************************************************************************************/
  typedef shape_t shape_fun_t;
  /*!***********************************************************************************************
   * \brief   The dimension of the Lebesque measure with respect to which we integrate.
   ************************************************************************************************/
  static constexpr unsigned int dim() { return shape_t::dim(); }
  /*!***********************************************************************************************
   * \brief   Calculate the amount of quadrature points.
   ************************************************************************************************/
  static constexpr unsigned int n_quad_points()
  {
    return Hypercube<dim()>::pow(quadrature_t::n_quad_points());
  }

  /*!***********************************************************************************************
   * \brief   Quadrature points on one-dimensional unit interval.
   ************************************************************************************************/
  static constexpr std::array<return_t, quadrature_t::n_points()> quad_points =
    quadrature_t::template points<return_t>();
  /*!***********************************************************************************************
   * \brief   Quadrature weights on one-dimensional unit interval.
   ************************************************************************************************/
  static constexpr std::array<return_t, quadrature_t::n_points()> quad_weights =
    quadrature_t::template weights<return_t>();

 private:
  /*!***********************************************************************************************
   * \brief   Amount of shape functions with respect to a single spatial dimension.
   ************************************************************************************************/
  static constexpr unsigned int n_fun_1D = shape_t::shape_fun_t::shape_fun_1d::n_fun();
  /*!***********************************************************************************************
   * \brief   Shape functions evaluated at quadrature points of unit interval.
   ************************************************************************************************/
  static constexpr std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D>
  shape_fcts_at_quad_points()
  {
    std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < quadrature_t::n_points(); ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template fct_val<return_t>(
          fct, quadrature_t::points()[pt]);

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Derivatives of shape functions evaluated at quadrature points of unit interval.
   ************************************************************************************************/
  std::array<std::array<return_t, quadrature_t::n_points()>,
             n_fun_1D> static constexpr shape_ders_at_quad_points()
  {
    std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < quadrature_t::n_points(); ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template der_val<return_t>(
          fct, quadrature_t::points()[pt]);

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Shape functions evaluated at boundaries of unit interval.
   ************************************************************************************************/
  std::array<std::array<return_t, 2>, n_fun_1D> static constexpr shape_fcts_at_bdrs()
  {
    std::array<std::array<return_t, 2>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < 2; ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template fct_val<return_t>(fct, pt);

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Derivatives of shape functions evaluated at boundaries of unit interval.
   ************************************************************************************************/
  std::array<std::array<return_t, 2>, n_fun_1D> static constexpr shape_ders_at_bdrs()
  {
    std::array<std::array<return_t, 2>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < 2; ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template der_val<return_t>(fct, pt);

    return result;
  }

 public:
  /*!***********************************************************************************************
   * \brief   Shape functions evaluated at quadrature points of unit interval.
   ************************************************************************************************/
  static constexpr std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D>
    shape_fcts_at_quad = shape_fcts_at_quad_points();
  /*!***********************************************************************************************
   * \brief   Derivatives of shape functions evaluated at quadrature points of unit interval.
   ************************************************************************************************/
  static constexpr std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D>
    shape_ders_at_quad = shape_ders_at_quad_points();
  /*!***********************************************************************************************
   * \brief   Shape functions evaluated at boundaries of unit interval.
   ************************************************************************************************/
  static constexpr std::array<std::array<return_t, 2>, n_fun_1D> shape_fcts_at_bdr =
    shape_fcts_at_bdrs();
  /*!***********************************************************************************************
   * \brief   Derivatives of shape functions evaluated at boundaries of unit interval.
   ************************************************************************************************/
  static constexpr std::array<std::array<return_t, 2>, n_fun_1D> shape_ders_at_bdr =
    shape_ders_at_bdrs();
  /*!***********************************************************************************************
   * \brief   Integrate product of one-dimensional shape functions.
   *
   * \param   i             Local index of local one-dimensional shape function.
   * \param   j             Local index of local one-dimensional shape function.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_1D_phiphi(const unsigned int i, const unsigned int j)
  {
    tpp_assert(i < shape_fcts_at_quad.size() && j < shape_fcts_at_quad.size(),
               "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights.size(); ++q)
      result += quad_weights[q] * shape_fcts_at_quad[i][q] * shape_fcts_at_quad[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of one-dimensional shape function and one derivative.
   *
   * \param   i             Local index of local one-dimensional shape function.
   * \param   j             Local index of local one-dimensional shape function (with derivative).
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_1D_phiDphi(const unsigned int i, const unsigned int j)
  {
    tpp_assert(i < shape_fcts_at_quad.size() && j < shape_fcts_at_quad.size(),
               "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights.size(); ++q)
      result += quad_weights[q] * shape_fcts_at_quad[i][q] * shape_ders_at_quad[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of one-dimensional shape function and one derivative.
   *
   * \param   i             Local index of local one-dimensional shape function (with derivative).
   * \param   j             Local index of local one-dimensional shape function.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_1D_Dphiphi(const unsigned int i, const unsigned int j)
  {
    tpp_assert(i < shape_fcts_at_quad.size() && j < shape_fcts_at_quad.size(),
               "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights.size(); ++q)
      result += quad_weights[q] * shape_ders_at_quad[i][q] * shape_fcts_at_quad[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of two one-dimensional shape functions' derivatives.
   *
   * \param   i             Local index of local one-dimensional shape function (with derivative).
   * \param   j             Local index of local one-dimensional shape function (with derivative).
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_1D_DphiDphi(const unsigned int i, const unsigned int j)
  {
    tpp_assert(i < shape_fcts_at_quad.size() && j < shape_fcts_at_quad.size(),
               "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights.size(); ++q)
      result += quad_weights[q] * shape_ders_at_quad[i][q] * shape_ders_at_quad[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over dimT-dimensional unit volume.
   *
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_vol_phiphi(const unsigned int i, const unsigned int j)
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function amd derivative over dimT-dimensional unit volume.
   *
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function (with derivative).
   * \param   dim_der       Dimension of the derivative.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_vol_phiDphi(const unsigned int i,
                                                  const unsigned int j,
                                                  const unsigned int dim_der)
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (dim_der == dim_fct)
        integral *= integrate_1D_phiDphi(dec_i[dim_fct], dec_j[dim_fct]);
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function and derivative over dimT-dimensional unit volume.
   *
   * \param   i             Local index of local shape function (with derivative).
   * \param   j             Local index of local shape function.
   * \param   dim_der       Dimension of the derivative.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_vol_Dphiphi(const unsigned int i,
                                                  const unsigned int j,
                                                  const unsigned int dim_der)
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (dim_der == dim_fct)
        integral *= integrate_1D_Dphiphi(dec_i[dim_fct], dec_j[dim_fct]);
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over dimT-dimensional volume's boundary.
   *
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_bdr_phiphi(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int bdr)
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (bdr_dim == dim_fct)
        integral *=
          shape_fcts_at_bdr[dec_i[dim_fct]][bdr_ind] * shape_fcts_at_bdr[dec_j[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function of volume times shape function of volume's face
   *          over dimT-dimensional volume's boundary.
   *
   * \param   i             Local index of local volume shape function.
   * \param   j             Local index of local boundary shape function.
   * \param   bdr           Boundary face index.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  static constexpr return_t integrate_bdr_phipsi(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int bdr)
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(dim() - 1, 1U)> dec_j =
      Hypercube<dim() - 1>::index_decompose(j, n_fun_1D);
    unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (bdr_dim == dim_fct)
        integral *= shape_fcts_at_bdr[dec_i[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct - (dim_fct > bdr_dim)]);
    return integral;
  }

  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions times some function over some geometry.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function that is also to be integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_vol_phiphifunc(const unsigned int i,
                                           const unsigned int j,
                                           geom_t& geom,
                                           const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, dim()> dec_q;
    smallVec_t quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_val *= quad_weights[dec_q[dim]] * shape_fcts_at_quad[dec_i[dim]][dec_q[dim]] *
                    shape_fcts_at_quad[dec_j[dim]][dec_q[dim]];
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), f_param) * quad_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions times some function over some geometry.
   *
   * \tparam  smallVec_t    Type of vector which is returned by the function.
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function that is also to be integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   dimension     Local dimension with respect to which the vector-function is integrated.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            point_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_vol_phiphivecfunc(const unsigned int i,
                                              const unsigned int j,
                                              const unsigned int dimension,
                                              geom_t& geom,
                                              const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, dim()> dec_q;
    smallVec_t quad_pt;
    const auto mat_q = geom.mat_q();

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_val *= quad_weights[dec_q[dim]] * shape_fcts_at_quad[dec_i[dim]][dec_q[dim]] *
                    shape_fcts_at_quad[dec_j[dim]][dec_q[dim]];
      }
      integral +=
        scalar_product(fun(geom.map_ref_to_phys(quad_pt), f_param), mat_q.get_column(dimension)) *
        quad_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over some geometry.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename geom_t>
  static return_t integrate_vol_phiphi(const unsigned int i, const unsigned int j, geom_t& geom)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of linear combinations of shape functions over some geometry.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  array_size    Size of arrays containing coefficients of linear combinations.
   * \tparam  floating_t    The floating point type for the calculation.
   * \param   is            Coefficients of local shape functions.
   * \param   js            Coefficients of local shape functions.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of lineat combinations of shape functions.
   ************************************************************************************************/
  template <typename geom_t, std::size_t array_size, typename floating_t>
  static return_t integrate_vol_phiphi(const std::array<floating_t, array_size>& is,
                                       const std::array<floating_t, array_size>& js,
                                       geom_t& geom)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val, is_val, js_val, val_helper;

    std::array<unsigned int, dim()> dec_k, dec_q;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      is_val = 0.;
      js_val = 0.;

      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
        quad_val *= quad_weights[dec_q[dim]];

      for (unsigned int k = 0; k < array_size; ++k)
      {
        dec_k = Hypercube<dim()>::index_decompose(k, n_fun_1D);
        val_helper = 1.;
        for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
          val_helper *= shape_fcts_at_quad[dec_k[dim]][dec_q[dim]];
        is_val += is[k] * val_helper;
        js_val += js[k] * val_helper;
      }
      integral += quad_val * is_val * js_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function times some function over some geometry.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function whose product with shape function is integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_vol_phifunc(const unsigned int i,
                                        geom_t& geom,
                                        const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_q;
    smallVec_t quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_val *= quad_weights[dec_q[dim]] * shape_fcts_at_quad[dec_i[dim]][dec_q[dim]];
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), f_param) * quad_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Average integral of product of shape function times some function over some geometry.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function whose product with shape function is integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_volUni_phifunc(const unsigned int i,
                                           geom_t& geom,
                                           const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_q;
    smallVec_t quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_val *= quad_weights[dec_q[dim]] * shape_fcts_at_quad[dec_i[dim]][dec_q[dim]];
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), f_param) * quad_val;
    }
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape function times other shape function over some geometry.
   *
   * \note    poly_deg_i and poly_deg_j must be set to the maximum polynomial degree (which is also
   *          their default value) if the basis is not hierarchical.
   *
   * \tparam  smallVec_t    Type of the return value.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  poly_deg_i    Polynomial degree of shape functions associated to i.
   * \tparam  poly_deg_j    Polynomial degree of shape functions associated to j.
   * \param   i             Local index of local shape function with gradient.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename smallVec_t,
            typename geom_t,
            unsigned int poly_deg_i = shape_t::degree(),
            unsigned int poly_deg_j = shape_t::degree()>
  static smallVec_t integrate_vol_nablaphiphi(const unsigned int i,
                                              const unsigned int j,
                                              geom_t& geom)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    static_assert(poly_deg_i <= shape_t::degree() && poly_deg_j <= shape_t::degree(),
                  "The maximum polynomial degrees must be larger than or equal to the given ones.");
    smallVec_t integral(1.);
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, poly_deg_i + 1);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, poly_deg_j + 1);
    for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
        if (dim == dim_fct)
          integral[dim] *= integrate_1D_Dphiphi(dec_i[dim_fct], dec_j[dim_fct]);
        else
          integral[dim] *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return geom.area() * integral / transposed(geom.mat_r());
  }
  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape functions times other function over some geometry.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   f_param       Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_vol_nablaphinablaphifunc(const unsigned int i,
                                                     const unsigned int j,
                                                     geom_t& geom,
                                                     return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_weight;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, dim()> dec_q;
    smallVec_t quad_pt, nabla_phi_i, nabla_phi_j;
    const auto rrT = mat_times_transposed_mat(geom.mat_r(), geom.mat_r());

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      nabla_phi_i = 1.;
      nabla_phi_j = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_weight *= quad_weights[dec_q[dim]];
        for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct)
          {
            nabla_phi_i[dim_fct] *= shape_ders_at_quad[dec_i[dim]][dec_q[dim]];
            nabla_phi_j[dim_fct] *= shape_ders_at_quad[dec_j[dim]][dec_q[dim]];
          }
          else
          {
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad[dec_i[dim]][dec_q[dim]];
            nabla_phi_j[dim_fct] *= shape_fcts_at_quad[dec_j[dim]][dec_q[dim]];
          }
        }
      }
      integral += quad_weight * fun(geom.map_ref_to_phys(quad_pt), f_param) *
                  scalar_product(nabla_phi_i, nabla_phi_j / rrT);
    }
    return geom.area() * integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate derivative of shape function times other function over some geometry.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   dim_der       Dimension with respect to which derivative is calculated.
   * \param   geom          Geometrical information.
   * \param   f_param       Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_vol_derphifunc(const unsigned int i,
                                           const unsigned int dim_der,
                                           geom_t& geom,
                                           return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_weight;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_q;
    smallVec_t quad_pt, nabla_phi_i;
    const auto rT = transposed(geom.mat_r());

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      nabla_phi_i = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_weight *= quad_weights[dec_q[dim]];
        for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct)
            nabla_phi_i[dim_fct] *= shape_ders_at_quad[dec_i[dim]][dec_q[dim]];
          else
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad[dec_i[dim]][dec_q[dim]];
        }
      }
      integral +=
        quad_weight * fun(geom.map_ref_to_phys(quad_pt), f_param) * (nabla_phi_i / rT)[dim_der];
    }
    return geom.area() * integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape function times shape function times other function times
   *          normal over some geometry's boundary.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function with gradient.
   * \param   j             Local index of local shape function.
   * \param   bdr           Index of the boundatry face to integrate over.
   * \param   geom          Geometrical information.
   * \param   f_param       Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_bdr_nablaphiphinufunc(const unsigned int i,
                                                  const unsigned int j,
                                                  const unsigned int bdr,
                                                  geom_t& geom,
                                                  return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_weight, phi_j;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, std::max(1U, dim() - 1)> dec_q;
    smallVec_t quad_pt, nabla_phi_i, normal;
    const auto rT = transposed(geom.mat_r());
    const unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<dim() - 1>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      phi_j = 1.;
      nabla_phi_i = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        if (dim == bdr_dim)
        {
          normal[dim] = 2. * bdr_ind - 1.;
          quad_pt[dim] = (return_t)bdr_ind;
          phi_j *= shape_fcts_at_bdr[dec_j[dim]][bdr_ind];
        }
        else
        {
          normal[dim] = 0.;
          quad_pt[dim] = quad_points[dec_q[dim - (dim > bdr_dim)]];
          phi_j *= shape_fcts_at_quad[dec_j[dim]][dec_q[dim - (dim > bdr_dim)]];
          quad_weight *= quad_weights[dec_q[dim - (dim > bdr_dim)]];
        }
        for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct && dim == bdr_dim)
            nabla_phi_i[dim_fct] *= shape_ders_at_bdr[dec_i[dim]][bdr_ind];
          else if (dim == dim_fct)
            nabla_phi_i[dim_fct] *= shape_ders_at_quad[dec_i[dim]][dec_q[dim - (dim > bdr_dim)]];
          else if (dim == bdr_dim)
            nabla_phi_i[dim_fct] *= shape_fcts_at_bdr[dec_i[dim]][bdr_ind];
          else
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad[dec_i[dim]][dec_q[dim - (dim > bdr)]];
        }
      }
      integral += quad_weight * fun(geom.map_ref_to_phys(quad_pt), f_param) * phi_j *
                  scalar_product(normal, nabla_phi_i / rT);
    }
    return geom.face_area(bdr) * integral;
  }

  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape function times shape function times other function times
   *          normal over some geometry's boundary.
   *
   * \tparam  point_t       Type of point which is the first argument of the function.
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function with gradient.
   * \param   j             Local index of local shape function.
   * \param   bdr           Index of the boundatry face to integrate over.
   * \param   geom          Geometrical information.
   * \param   f_param       Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_bdr_nablaphipsinufunc(const unsigned int i,
                                                  const unsigned int j,
                                                  const unsigned int bdr,
                                                  geom_t& geom,
                                                  return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_weight, phi_j;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(1U, dim() - 1)> dec_j =
      Hypercube<dim() - 1>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, std::max(1U, dim() - 1)> dec_q;
    smallVec_t quad_pt, nabla_phi_i, normal;
    const auto rT = transposed(geom.mat_r());
    const unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<dim() - 1>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      phi_j = 1.;
      nabla_phi_i = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        if (dim == bdr_dim)
        {
          normal[dim] = 2. * bdr_ind - 1.;
          quad_pt[dim] = (return_t)bdr_ind;
        }
        else
        {
          normal[dim] = 0.;
          quad_pt[dim] = quad_points[dec_q[dim - (dim > bdr_dim)]];
          phi_j *= shape_fcts_at_quad[dec_j[dim - (dim > bdr_dim)]][dec_q[dim - (dim > bdr_dim)]];
          quad_weight *= quad_weights[dec_q[dim - (dim > bdr_dim)]];
        }
        for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct && dim == bdr_dim)
            nabla_phi_i[dim_fct] *= shape_ders_at_bdr[dec_i[dim]][bdr_ind];
          else if (dim == dim_fct)
            nabla_phi_i[dim_fct] *= shape_ders_at_quad[dec_i[dim]][dec_q[dim - (dim > bdr_dim)]];
          else if (dim == bdr_dim)
            nabla_phi_i[dim_fct] *= shape_fcts_at_bdr[dec_i[dim]][bdr_ind];
          else
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad[dec_i[dim]][dec_q[dim - (dim > bdr_dim)]];
        }
      }
      integral += quad_weight * fun(geom.map_ref_to_phys(quad_pt), f_param) * phi_j *
                  scalar_product(normal, nabla_phi_i / rT);
    }
    return geom.face_area(bdr) * integral;
  }

  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over boundary face.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename geom_t>
  static return_t integrate_bdr_phiphi(const unsigned int i,
                                       const unsigned int j,
                                       const unsigned int bdr,
                                       geom_t& geom)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    unsigned int dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
      if (dim == dim_fct)
        integral *=
          shape_fcts_at_bdr[dec_i[dim_fct]][bdr_ind] * shape_fcts_at_bdr[dec_j[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral * geom.face_area(bdr);
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions of volumen and skeletal over boundary face.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \param   i             Local index of local volumne shape function.
   * \param   j             Local index of local skeletal shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename geom_t>
  static return_t integrate_bdr_phipsi(const unsigned int i,
                                       const unsigned int j,
                                       const unsigned int bdr,
                                       geom_t& geom)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(1U, dim() - 1)> dec_j =
      Hypercube<dim() - 1>::index_decompose(j, n_fun_1D);
    unsigned int dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < geom_t::hyEdge_dim(); ++dim_fct)
      if (dim == dim_fct)
        integral *= shape_fcts_at_bdr[dec_i[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct - (dim_fct > dim)]);
    return integral * geom.face_area(bdr);
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions times some function over boundary face.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function that is multiplied by shape function.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_bdr_phifunc(const unsigned int i,
                                        const unsigned int bdr,
                                        geom_t& geom,
                                        const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(1U, dim() - 1)> dec_q;
    smallVec_t quad_pt;
    unsigned int dim_bdr = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<dim() - 1>::index_decompose(q, n_fun_1D);
      quad_val = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        if (dim == dim_bdr)
        {
          quad_pt[dim] = bdr_ind;
          quad_val *= shape_fcts_at_bdr[dec_i[dim]][bdr_ind];
        }
        else
        {
          quad_pt[dim] = quad_points[dec_q[dim - (dim > dim_bdr)]];
          quad_val *= quad_weights[dec_q[dim - (dim > dim_bdr)]] *
                      shape_fcts_at_quad[dec_i[dim]][dec_q[dim - (dim > dim_bdr)]];
        }
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), f_param) * quad_val;
    }
    return integral * geom.face_area(bdr);
  }
  /*!***********************************************************************************************
   * \brief   Average integral of product of skeletal shape functions times some function.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function that is multiplied by shape function.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   i             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t = point_t>
  static return_t integrate_bdrUni_psifunc(const unsigned int i,
                                           const unsigned int bdr,
                                           geom_t& geom,
                                           const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_val;
    std::array<unsigned int, std::max(1U, dim() - 1)> dec_q,
      dec_i = Hypercube<dim() - 1>::index_decompose(i, n_fun_1D);
    smallVec_t quad_pt;
    unsigned int dim_bdr = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<dim() - 1>::index_decompose(q, n_fun_1D);
      quad_val = 1.;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        if (dim == dim_bdr)
          quad_pt[dim] = bdr_ind;
        else
        {
          quad_pt[dim] = quad_points[dec_q[dim - (dim > dim_bdr)]];
          quad_val *=
            quad_weights[dec_q[dim - (dim > dim_bdr)]] *
            shape_fcts_at_quad[dec_i[dim - (dim > dim_bdr)]][dec_q[dim - (dim > dim_bdr)]];
        }
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), f_param) * quad_val;
    }
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Squared L2 distance of some function and an discrete function on volume.
   *
   * \tparam  geom_t        Geometry which is the integration domain.
   * \tparam  fun           Function whose distance is measured.
   * \tparam  smallVec_t    Type of local point with respect to hyperedge. Defaults to point_t.
   * \param   coeffs        Coefficients of discrete function.
   * \param   geom          Geometrical information.
   * \param   f_param       Function parameter (e.g. time) with respect to which it is evaluated.
   * \retval  integral      Squared distance of functions.
   ************************************************************************************************/
  template <typename point_t,
            typename geom_t,
            return_t fun(const point_t&, const return_t),
            typename smallVec_t,
            std::size_t n_coeff>
  static return_t integrate_vol_diffsquare_discana(const std::array<return_t, n_coeff> coeffs,
                                                   geom_t& geom,
                                                   const return_t f_param = 0.)
  {
    static_assert(geom_t::hyEdge_dim() == dim(), "Dimension of hyperedge must fit to quadrature!");
    return_t integral = 0., quad_weight;
    std::array<unsigned int, dim()> dec_i, dec_q;
    std::array<return_t, n_coeff> quad_val;
    smallVec_t quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights.size(), geom_t::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<dim()>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      quad_val = coeffs;
      for (unsigned int dim = 0; dim < geom_t::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points[dec_q[dim]];
        quad_weight *= quad_weights[dec_q[dim]];
        for (unsigned int i = 0; i < n_coeff; ++i)
        {
          dec_i = Hypercube<geom_t::hyEdge_dim()>::index_decompose(i, n_fun_1D);
          quad_val[i] *= shape_fcts_at_quad[dec_i[dim]][dec_q[dim]];
        }
      }
      integral += quad_weight * std::pow(fun(geom.map_ref_to_phys(quad_pt), f_param) -
                                           std::accumulate(quad_val.begin(), quad_val.end(), 0.),
                                         2);
    }
    return integral * geom.area();
  }
};  // end of class Integrator

}  // end of namespace Quadrature

}  // end of namespace TPP
