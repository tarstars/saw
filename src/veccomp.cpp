#include <cmath>

using namespace std;

#include "veccomp.h"

CD
VecComp3::abs()const 
{
  return sqrt(x * conj(x) + y * conj(y) + z * conj(z)); 
}

VecComp3 operator&(const VecComp3 &r1, const VecComp3 &r2) 
{ 
  return VecComp3(r1.y *(r2.z) - r1.z *(r2.y), 
		  r1.z *(r2.x) - r1.x *(r2.z), 
		  r1.x *(r2.y) - r1.y *(r2.x));
}

poly_type::value_type operator*(const VecComp3 &r1, const VecComp3 &r2)
{
  return r1.x * r2.x + r1.y * r2.y + r1.z * r2.z;
}

VecComp3 operator*(const VecComp3 &r1, poly_type::value_type r2)
{
    return VecComp3(r1.x * r2, r1.y * r2, r1.z * r2);
}

VecComp3 operator*(poly_type::value_type r1, const VecComp3 &r2) {return r2*r1;}

VecComp3 operator+(const VecComp3 &r1, const VecComp3 &r2)
{
    return VecComp3(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z);
}

ostream& operator<<(ostream& os, const VecComp3 &r)
{
    return os<<r.x<<" "<<r.y<<" "<<r.z;
}
