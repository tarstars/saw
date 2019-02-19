#ifndef _VEC3_
#define _VEC3_

#include <iostream>
#include <vector>
#include "types.h"

class Vec3
{
public:
  std::vector<LD> dat;
  
  explicit Vec3(LD ax = 0, LD ay = 0, LD az = 0): 
  dat(3)
  {
    dat[0] = ax; 
    dat[1] = ay; 
    dat[2] = az;
  }

  LD abs()const 
  {
    LD s = 0; 
    for(int t = 0; t < 3; ++t)
      s += dat[t] * dat[t];
    return sqrt(s); }

  Vec3 norm()const 
  {LD r=abs(); return Vec3(dat[0] / r, dat[1] / r, dat[2] / r);}

  LD& operator[](int t);
  const LD& operator[](int t)const{return dat[t];}

  static Vec3 theta_phi2vec(double theta, double phi);
};

Vec3 operator&(const Vec3 &r1, const Vec3 &r2);
LD operator*(const Vec3 &r1, const Vec3 &r2);
Vec3 operator*(const Vec3 &r1, CD r2);
Vec3 operator*(CD r1, const Vec3 &r2);
Vec3 operator+(const Vec3 &r1, const Vec3 &r2);

std::ostream& operator<<(std::ostream& os, const Vec3 &r);


#endif
