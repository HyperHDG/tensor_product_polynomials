#include <tpp/shape_function/one_dimensional.hxx>
#include <tpp/tpp_assert.hxx>

#include <array>
#include <cmath>

using namespace std;
using namespace TPP;
using namespace TPP::ShapeType;

template <typename float_t>
int legendre()
{
  auto legendre_val = [](const float_t x, const unsigned int index) -> float_t
  {
    switch (index)
    {
      case 0:
        return 1.;
      case 1:
        return -sqrt(3) * (1. - 2. * x);
      case 2:
        return sqrt(5) * ((6. * x - 6.) * x + 1.);
      case 3:
        return sqrt(7) * (((20. * x - 30.) * x + 12.) * x - 1.);
      case 4:
        return sqrt(9) * ((((70. * x - 140.) * x + 90.) * x - 20.) * x + 1.);
      case 5:
        return sqrt(11) * (((((252 * x - 630.) * x + 560.) * x - 210.) * x + 30.) * x - 1.);
    }
    return 0.;
  };

  auto legendre_der = [](const float_t x, const unsigned int index) -> float_t
  {
    switch (index)
    {
      case 0:
        return 0.;
      case 1:
        return sqrt(12);
      case 2:
        return sqrt(5) * (12. * x - 6.);
      case 3:
        return sqrt(7) * ((60. * x - 60.) * x + 12.);
      case 4:
        return sqrt(9) * (((280. * x - 420.) * x + 180.) * x - 20.);
      case 5:
        return sqrt(11) * ((((1260. * x - 2520.) * x + 1680.) * x - 420.) * x + 30.);
    }
    return 0.;
  };

  std::array<float_t, 5> x_val = {0., 0.25, 0.5, 0.75, 1.};
  float_t distance, tolerance = 1e-13;
  bool successful = true;

  for (unsigned int deg = 0; deg < 6; ++deg)
    for (unsigned int val = 0; val < x_val.size(); ++val)
    {
      distance = abs(legendre_val(x_val[val], deg) -
                     Legendre<5>::template fct_val<float_t>(deg, x_val[val]));
      tpp_assert(distance < tolerance,
                 ("Legendre evaluation failed for polynomial index " + to_string(deg)).c_str());
      if (distance >= tolerance)
        successful = false;
    }

  for (unsigned int deg = 0; deg < 6; ++deg)
    for (unsigned int val = 0; val < x_val.size(); ++val)
    {
      distance = abs(legendre_der(x_val[val], deg) -
                     Legendre<5>::template der_val<float_t>(deg, x_val[val]));
      tpp_assert(distance < tolerance,
                 ("Legendre evaluation failed for polynomial index " + to_string(deg)).c_str());
      if (distance >= tolerance)
        successful = false;
    }

  return successful - 1;
}

template <typename float_t>
int lobatto()
{
  auto lobatto_val = [](const float_t x, const unsigned int index) -> float_t
  {
    switch (index)
    {
      case 0:
        return 1.;
      case 1:
        return x;
      case 2:
        return -sqrt(3) * (x - pow(x, 2));
      case 3:
        return sqrt(5) * (2. * pow(x, 3) - 3. * pow(x, 2) + x);
      case 4:
        return sqrt(7) * (5. * pow(x, 4) - 10. * pow(x, 3) + 6. * pow(x, 2) - x);
      case 5:
        return sqrt(9) * (14. * pow(x, 5) - 35. * pow(x, 4) + 30. * pow(x, 3) - 10 * pow(x, 2) + x);
      case 6:
        return sqrt(11) * (42. * pow(x, 6) - 126 * pow(x, 5) + 140. * pow(x, 4) - 70. * pow(x, 3) +
                           15. * pow(x, 2) - x);
    }
    return 0.;
  };

  auto lobatto_der = [](const float_t x, const unsigned int index) -> float_t
  {
    switch (index)
    {
      case 0:
        return 0.;
      case 1:
        return 1.;
      case 2:
        return -sqrt(3) * (1. - 2. * x);
      case 3:
        return sqrt(5) * ((6. * x - 6.) * x + 1.);
      case 4:
        return sqrt(7) * (((20. * x - 30.) * x + 12.) * x - 1.);
      case 5:
        return sqrt(9) * ((((70. * x - 140.) * x + 90.) * x - 20.) * x + 1.);
      case 6:
        return sqrt(11) * (((((252 * x - 630.) * x + 560.) * x - 210.) * x + 30.) * x - 1.);
    }
    return 0.;
  };

  std::array<float_t, 5> x_val = {0., 0.25, 0.5, 0.75, 1.};
  float_t distance, tolerance = 1e-14;
  bool successful = true;

  for (unsigned int deg = 0; deg < 7; ++deg)
    for (unsigned int val = 0; val < x_val.size(); ++val)
    {
      distance =
        abs(lobatto_val(x_val[val], deg) - Lobatto<6>::template fct_val<float_t>(deg, x_val[val]));
      tpp_assert(distance < tolerance,
                 ("Legendre evaluation failed for polynomial index " + to_string(deg)).c_str());
      if (distance >= tolerance)
        successful = false;
    }

  for (unsigned int deg = 0; deg < 7; ++deg)
    for (unsigned int val = 0; val < x_val.size(); ++val)
    {
      distance =
        abs(lobatto_der(x_val[val], deg) - Lobatto<6>::template der_val<float_t>(deg, x_val[val]));
      tpp_assert(distance < tolerance,
                 ("Legendre evaluation failed for polynomial index " + to_string(deg)).c_str());
      if (distance >= tolerance)
        successful = false;
    }

  return successful - 1;
}

int main()
{
  int result = 0;

  result += legendre<double>() + lobatto<double>();

  return result;
}
