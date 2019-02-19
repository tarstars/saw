#ifndef _VECCOMP_
#define _VECCOMP_

#include <iostream>

#include "types.h"

class VecComp3
{
public:
  poly_type::value_type x, y, z;
  explicit VecComp3(poly_type::value_type ax = 0, poly_type::value_type ay = 0, poly_type::value_type az = 0): x(ax), y(ay), z(az) {}

  CD abs()const;

  VecComp3 norm()const 
  {
    poly_type::value_type r=abs(); 
    return VecComp3(x / r, y / r, z / r);
  }

  poly_type::value_type operator()(int t){return (t == 0) ? x : ((t==1)?y:z);}
};

CD operator*(const VecComp3& l, const VecComp3& r);
VecComp3 operator&(const VecComp3 &r1, const VecComp3 &r2);
poly_type::value_type operator*(const VecComp3 &r1, const VecComp3 &r2);
VecComp3 operator*(const VecComp3 &r1, poly_type::value_type r2);
VecComp3 operator*(poly_type::value_type r1, const VecComp3 &r2);
VecComp3 operator+(const VecComp3 &r1, const VecComp3 &r2);
std::ostream& operator<<(std::ostream& os, const VecComp3 &r);

#endif
