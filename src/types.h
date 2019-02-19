#ifndef _TYPES_
#define _TYPES_

#include <complex>
#include <vector>

typedef long double LD;

typedef std::complex<LD> CD;

typedef std::vector<CD> poly_type;

typedef std::vector<double> VD;

namespace 
{
  const CD rho = 6.0;
  const long double eps = 1e-10;
}

#endif
