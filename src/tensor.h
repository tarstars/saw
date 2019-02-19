#ifndef _TENSOR_
#define _TENSOR_

#include "types.h"

#include <iostream>

class Coeff;
class Mat3;
class Vec3;
class Poly;

class Tensor
{
 public:
  CD dat[3][3][3][3];
  Tensor();
  Tensor(const Coeff&);
  Tensor rotate(const Mat3&);
  Tensor rotate_bidir(const Mat3& st, const Mat3& inv);
  void test(std::ostream&);
  Mat3 crist(const Vec3&)const;
  Poly poly_by_n3(double n1, double n2, double g)const;
  void mupad_style_output(std::ostream&)const;
};

std::ostream& operator<<(std::ostream& os, const Tensor&);

#endif
