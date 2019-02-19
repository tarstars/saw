#ifndef _MAT3_
#define _MAT3_

#include <iostream>
#include "types.h"

class Vec3;

class Mat3
{
 public:
  CD dat[3][3];
  Mat3();
  Mat3(double, int);
  CD* operator[](int ind){return &dat[ind][0];}
  const CD* operator[](int ind)const{return &dat[ind][0];}
  Vec3 row(int)const;
  CD det()const;
};

Vec3 operator*(const Mat3 &l, const Vec3 &r);
Mat3 operator*(const Mat3 &l, const Mat3 &r);

std::ostream& operator<<(std::ostream&, const Mat3&);
#endif
