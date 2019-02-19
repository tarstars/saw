#include <complex>
#include <iostream>

#include "mat3.h"
#include "vec3c.h"
#include "vec3.h"

using namespace farn;
using namespace std;

Vec3c::Vec3c(const Vec3& r) : dat(3) {
  dat[0] = r.x();
  dat[1] = r.y();
  dat[2] = r.z();
}

ostream& 
farn::operator<<(std::ostream& os, const Vec3c &r) {
  return os << r.x() << " " << r.y() << " " << r.z();
}

double
Vec3c::getSuperModule()const {
  return sqrt(real(x() * conj(x()) + y() * conj(y()) + z() * conj(z())));
}

Vec3c
Vec3c::operator&(const Vec3c& r) const {
  complex<double> r1x = dat[0];
  complex<double> r1y = dat[1];
  complex<double> r1z = dat[2];
  
  complex<double> r2x = r.dat[0];
  complex<double> r2y = r.dat[1];
  complex<double> r2z = r.dat[2];

  return Vec3c(r1y * r2z - r1z * r2y, r1z * r2x - r1x * r2z, r1x * r2y - r1y * r2x);
}

void
Vec3c::normalize() {
  double m = abs();

  dat[0] /= m;
  dat[1] /= m;
  dat[2] /= m;
}

double 
Vec3c::abs() const {
  return sqrt(real(dat[0] * conj(dat[0]) + dat[1] * conj(dat[1]) + dat[2] * conj(dat[2])));
}

Vec3c&
Vec3c::operator+=(const Vec3c& r) {
  dat[0] += r.dat[0];
  dat[1] += r.dat[1];
  dat[2] += r.dat[2];

  return *this;
}

Vec3c 
farn::operator*(const Vec3c& vec, const Mat3& mat) {
  Vec3c ret;

  for(int q=0; q<3; q++) {
    ret[q] = 0;
    for(int p=0; p<3; p++) {
      ret[q] += vec[p] * mat[p][q];
    }
  }

  return ret;
}
