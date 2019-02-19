#include <cmath>
#include <iostream>

using namespace std;

#include "poly3.h"
#include "types.h"
#include "vec3.h"

CD 
CPoly3::find_root() const 
{
    CD x  = 1;
    CD dx = 0;
    int meter = 0;
    int max_meter = 10;

    do{
        ++meter;
        dx = val(x) / der_val(x);
        x -= dx;
    }while(abs(dx) > eps && meter < max_meter);

    return x;
}

void
CPoly3::find_two_root(CD *px1, CD *px2, CD x) const
{
    CD a = a3;
    CD b = a2 + x * a3;
    CD c = a1 + x * (a2 + x * a3);

    CD D = b * b - CD(4) * a * c;
    if(real(D) < 0) 
      {
	throw(string("Imaginary roots"));
	return;
      }

    *px1 = ( -b + sqrt(D))/(CD(2) * a);
    *px2 = ( -b - sqrt(D))/(CD(2) * a);
}

Vec3 polariz(Vec3 r1, Vec3 r2, CD x)
{
  r1.dat[0] -= real(x);
  r2.dat[1] -= real(x);

  return (r1 & r2).norm();
}

std::ostream& operator<<(std::ostream& os, const CPoly3& r)
{
  return os << r.a0 << " " << r.a1 << " " << r.a2 << " " << r.a3;
}
