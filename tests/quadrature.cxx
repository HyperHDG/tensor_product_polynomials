#include <tpp/quadrature/one_dimensional.hxx>

#include <cmath>
#include <cstdio>
#include <limits>
#include <utility>

#if __has_include(<stdfloat>)
#include <stdfloat>
#endif

// Gauss--Legendre rules with 1 to 12 points in double, long double and (where the compiler has it)
// std::float128_t: points ascending in (0, 1), positive weights summing to one, and every monomial
// x^k with k <= 2 n - 1 integrated exactly, up to a few epsilon of the type.

int failures = 0;

template <typename T>
T absolute(const T x)
{
  return x < (T)0. ? -x : x;
}

template <unsigned int n_points, typename T>
void check(const char* type)
{
  using Rule = TPP::Quadrature::GaussLegendre<2 * n_points - 1>;
  static_assert(Rule::n_points() == n_points);
  constexpr auto points = Rule::template points<T>();
  constexpr auto weights = Rule::template weights<T>();
  const T tolerance = (T)16. * std::numeric_limits<T>::epsilon();
  bool ok = true;
  T sum = 0.;
  for (unsigned int i = 0; i < n_points; ++i)
  {
    ok = ok && points[i] > (T)0. && points[i] < (T)1. && weights[i] > (T)0.;
    ok = ok && (i == 0 || points[i - 1] < points[i]);
    sum += weights[i];
  }
  T worst = absolute(sum - (T)1.);
  for (unsigned int k = 0; k <= 2 * n_points - 1; ++k)
  {
    T integral = 0.;
    for (unsigned int i = 0; i < n_points; ++i)
    {
      T power = 1.;
      for (unsigned int j = 0; j < k; ++j)
        power *= points[i];
      integral += weights[i] * power;
    }
    worst = std::max(worst, absolute(integral - (T)1. / (T)(k + 1)));
  }
  ok = ok && worst <= tolerance;
  if (!ok)
    std::printf("FAIL  %2u points, %s: deviation %.3e (tolerance %.3e)\n", n_points, type,
                (double)worst, (double)tolerance);
  failures += !ok;
}

template <typename T>
void check_all(const char* type)
{
  [&]<unsigned int... n>(std::integer_sequence<unsigned int, n...>)
  { (check<n + 1, T>(type), ...); }(std::make_integer_sequence<unsigned int, 12>());
}

int main()
{
  check_all<double>("double");
  check_all<long double>("long double");
#if defined(__STDCPP_FLOAT128_T__)
  check_all<std::float128_t>("float128_t");
  std::printf("float128_t checked\n");
#endif
  std::printf("%s\n", failures ? "FAILED" : "ok");
  return failures;
}
