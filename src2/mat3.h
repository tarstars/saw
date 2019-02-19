#ifndef _MAT3_
#define _MAT3_

#include <iostream>
#include "types.h"

namespace farn{

  class Vec3;
  class Vec3c;

  class Mat3
  {
  public:
    CD dat[3][3];
    Mat3();
    Mat3(double, int);
    Mat3(const Vec3c&, const Vec3c&, const Vec3c&);
    CD* operator[](int ind){return &dat[ind][0];}
    const CD* operator[](int ind) const{return &dat[ind][0];}
    Vec3c row(int) const;
    Vec3c calcPol() const;
    CD det() const;
    bool isZeroDet() const;
  };
  
  Vec3c operator*(const Mat3 &l, const Vec3 &r);
  Vec3c operator*(const Mat3 &l, const Vec3c &r);
  Mat3 operator*(const Mat3 &l, const Mat3 &r);

  std::ostream& operator<<(std::ostream&, const Mat3&);
}
#endif
