#include <cmath>
#include <iostream>

#include "mat3.h"
#include "vec3.h"

using namespace std;

namespace farn{

  Vec3 operator&(const Vec3 &r1, const Vec3 &r2) 
  { 
    double r1x = r1.dat[0];
    double r1y = r1.dat[1];
    double r1z = r1.dat[2];

    double r2x = r2.dat[0];
    double r2y = r2.dat[1];
    double r2z = r2.dat[2];

    return Vec3(r1y * r2z - r1z * r2y, r1z * r2x - r1x * r2z, r1x * r2y - r1y * r2x);
  }

  double operator*(const Vec3 &r1, const Vec3 &r2)
  {
    double s = 0; 

    for(int p = 0; p < 3; ++p)
      s += r1.dat[p] * r2.dat[p];

    return s;
  }

  Vec3 operator*(const Vec3 &r1, double r)
  {
    Vec3 ret;

    for(int p = 0; p < int(r1.dat.size()); ++p)
      ret[p] = r1.dat[p] * r;

    return ret;
  }



  Vec3 operator+(const Vec3 &r1, const Vec3 &r2)
  {
    Vec3 ret;

    for(int p = 0; p < 3; ++p)
      ret.dat[p] = r1.dat[p] + r2.dat[p];

    return ret;
  }

  Vec3
  Vec3::theta_phi2vec(double theta, double phi)
  {
    return Vec3(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta) );
  }

  ostream& operator<<(ostream& os, const Vec3 &r)
  {
    return os << r.dat[0] << " " << r.dat[1] << " " << r.dat[2];
  }

  double& 
  Vec3::operator[](int t)
  {
    return dat[t];
  }
}
