#ifndef _MAT3POLY_
#define _MAT3POLY_

#include <iostream>

#include "poly.h"
#include "types.h"

class Mat3poly
{
 public:
  std::vector<Poly> dat;
 Mat3poly():dat(9){}
  std::vector<Poly>::iterator operator[](int ind){return dat.begin() + ind * 3;}
  std::vector<Poly>::const_iterator operator[](int ind)const{return dat.begin() + ind * 3;}
};

#endif
