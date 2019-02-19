#include <cmath>

#include "mat3.h"
#include "vec3.h"
#include "vec3c.h"
#include "util.h"

using namespace std;
namespace farn{

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

  Mat3::Mat3(const Vec3c &r1, const Vec3c &r2, const Vec3c &r3){
    Mat3& tst = *this;
    tst[0][0]= r1[0]; tst[0][1]= r1[1];    tst[0][2]= r1[2];
    tst[1][0]= r2[0];  tst[1][1]= r2[1];    tst[1][2]= r2[2];
    tst[2][0]= r3[0];  tst[2][1]= r3[1];    tst[2][2]= r3[2];
  }

  Vec3c
  Mat3::row(int ind) const {
    return Vec3c(dat[ind][0], dat[ind][1], dat[ind][2]);
  }

  Vec3c
  operator*(const Mat3 &l, const Vec3 &r) {
    Vec3c ret;
    int p,q;

    for(p = 0; p < 3; ++p) {
      for(q = 0; q < 3; ++q) {
	ret[p] += l[p][q] * r[q];
      }
    }

    return ret;
  }

  Vec3c
  operator*(const Mat3 &l, const Vec3c &r)
  {

    Vec3c ret;
    int p,q;

    for(p = 0; p < 3; ++p) {
      for(q = 0; q < 3; ++q) {
	ret[p] += l[p][q] * r[q];
      }
    }
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
    for(int t = 0; t < 6; ++t) {
      CD summary(prog[t][6]);
      for(int p = 0; p < 3; ++p) {
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

    for(p = 0; p < 3; ++p) {
      for(q = 0; q < 3; ++q) {
	for(r = 0; r < 3; ++r) {
	  ret[p][q] += lm[p][r] * rm[r][q];
	}
      }
    }

    return ret;
  }

  Vec3c
  Mat3::calcPol()const{
    #ifdef DEBUG
    cerr << "calcPol" << endl;
    #endif

    Vec3c a1(dat[0][0], dat[0][1], dat[0][2]);
    Vec3c a2(dat[1][0], dat[1][1], dat[1][2]);
    Vec3c a3(dat[2][0], dat[2][1], dat[2][2]);

    a1.normalize();
    a2.normalize();
    a3.normalize();

    #ifdef DEBUG
    cerr << "\t a1 = " << a1 << endl;
    cerr << "\t a2 = " << a2 << endl;
    cerr << "\t a3 = " << a3 << endl << endl;
    #endif

    Vec3c r1 = (a1 & a2);
    Vec3c r2 = (a2 & a3);
    Vec3c r3 = (a3 & a1);

    #ifdef DEBUG
    cerr << "\t r1 = " << r1 << endl;
    cerr << "\t r2 = " << r2 << endl;
    cerr << "\t r3 = " << r3 << endl << endl;
    #endif

    double v1 = r1.abs();
    double v2 = r2.abs();
    double v3 = r3.abs();

    #ifdef DEBUG
    cerr << "\tv1 = " << v1 << endl;
    cerr << "\tv2 = " << v2 << endl;
    cerr << "\tv3 = " << v3 << endl << endl;
    #endif

    if (v1 >= max(v2, v3)){
      #ifdef DEBUG
      cerr << "\t a1 a2 selected" << endl;
      #endif
      return r1.norm();
    }

    if (v2 >= max(v1, v3)){
      #ifdef DEBUG
      cerr << "\t a2 a3 selected" << endl;
      #endif
      return r2.norm();
    }

    #ifdef DEBUG
    cerr << "\t a3 a1 selected" << endl;
    #endif
    return r3.norm();
  }

  bool
  Mat3::isZeroDet() const {
    Mat3 cp;

    double maxVal = 0;
    for(int p = 0; p < 3; ++p) {
      for(int q = 0; q < 3; ++q) {
	maxVal = max(maxVal, abs(dat[p][q]));
      }
    }

    for(int p = 0; p < 3; ++p) {
      for(int q = 0; q < 3; ++q) {
	cp.dat[p][q] = dat[p][q] / maxVal;
      }
    }

    double val = abs(cp.det());

    return val < 5e-2;
  }
}

