#ifndef _POLY_
#define _POLY_

#include <vector>
#include <complex>
#include <iostream>

#include "types.h"

class Poly
{
 public:
  poly_type dat;

  Poly(){}

 Poly(poly_type::value_type a0):dat(1){dat[0] = a0;}
 Poly(poly_type::value_type a0, poly_type::value_type a1):dat(2){dat[0] = a0;dat[1] = a1;}
 Poly(poly_type::value_type a0, poly_type::value_type a1, poly_type::value_type a2):dat(3)
    {dat[0] = a0; dat[1] = a1; dat[2] = a2;} 

  Poly operator/(const Poly&)const;
  Poly operator*(const Poly&)const;
  Poly operator+(const Poly&)const;
  Poly operator-(const Poly&)const;
  poly_type::value_type operator()(poly_type::value_type)const;
  
  template<typename T>
    void fst(T x, T *pz, T *pf, T *ps)const;
  int deg()const{return dat.size() - 1;}
};

std::ostream& operator<<(std::ostream& os, const Poly& r);
poly_type all_roots(const Poly& r);
#endif
