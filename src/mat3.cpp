#include <cmath>

using namespace std;

#include "mat3.h"
#include "vec3.h"

Mat3::Mat3()
{
  for(int p = 0; p < 3; ++p)
    for(int q = 0; q < 3; ++q)
      dat[p][q] = 0;
}

Mat3::Mat3(double phi, int ind)
{
  dat[(ind + 0) % 3][(ind + 0) % 3] = 1; 
  dat[(ind + 0) % 3][(ind + 1) % 3] = 0; 
  dat[(ind + 0) % 3][(ind + 2) % 3] = 0;

  dat[(ind + 1) % 3][(ind + 0) % 3] = 0; 
  dat[(ind + 1) % 3][(ind + 1) % 3] = cos(phi); 
  dat[(ind + 1) % 3][(ind + 2) % 3] = -sin(phi);

  dat[(ind + 2) % 3][(ind + 0) % 3] = 0; 
  dat[(ind + 2) % 3][(ind + 1) % 3] = sin(phi); 
  dat[(ind + 2) % 3][(ind + 2) % 3] = cos(phi);
}

Vec3
Mat3::row(int ind)const
{
  return Vec3(real(dat[ind][0]), real(dat[ind][1]), real(dat[ind][2]));
}

Vec3
operator*(const Mat3 &l, const Vec3 &r)
{

  Vec3 ret;
  int p,q;

  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      ret[p] += real(l[p][q]) * r[q];

  return ret;
}

ostream& operator<<(ostream& os, const Mat3 &r)
{
  int p, q;
  for(p = 0; p < 3; ++p)
    {
      for(q = 0; q < 3; ++q)
	os << r.dat[p][q] << " ";
      os << endl;
    }
  return os;
}

CD
Mat3::det()const
{
  int prog[6][7] = 
    {
      {0,0,1,1,2,2,1},
      {0,1,1,2,2,0,1},
      {0,2,1,0,2,1,1},
      {2,0,1,1,0,2,-1},
      {2,1,1,2,0,0,-1},
      {2,2,1,0,0,1,-1}};

  CD ret;
  for(int t = 0; t < 6; ++t)
    {
      CD summary(prog[t][6]);
      for(int p = 0; p < 3; ++p)
	{
	  summary = summary * dat[prog[t][2 * p]][prog[t][2 * p + 1]];
	}
      ret = ret + summary;
    }
  return ret;
}


Mat3 operator*(const Mat3 &lm, const Mat3 &rm)
{
  Mat3 ret;

  int p;
  int q;
  int r;

  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      for(r = 0; r < 3; ++r)
	ret[p][q] += lm[p][r] * rm[r][q];

  return ret;
}
