#ifndef _VEC3_
#define _VEC3_

#include <iostream>
#include <vector>
#include <complex>

namespace farn{

  class Mat3;

  class Vec3
  {
  public:
    std::vector<double> dat;
  
    explicit Vec3(double ax = 0, double ay = 0, double az = 0): 
    dat(3)
      {
	dat[0] = ax; 
	dat[1] = ay; 
	dat[2] = az;
      }

    double abs()const 
    {
      double s = 0; 
      for(int t = 0; t < 3; ++t)
	s += dat[t] * dat[t];
      return sqrt(s); }

    Vec3 norm()const 
    {double r=abs(); return Vec3(dat[0] / r, dat[1] / r, dat[2] / r);}

    void normalize(){
      double r = abs();
      for(int t = 0; t < int(dat.size()); ++t)
	dat[t] /= r;
    }

    double& x(){return dat[0];}
    double& y(){return dat[1];}
    double& z(){return dat[2];}

    const double& x()const{return dat[0];}
    const double& y()const{return dat[1];}
    const double& z()const{return dat[2];}

    double& operator[](int t);
    const double& operator[](int t)const{return dat[t];}

    static Vec3 theta_phi2vec(double theta, double phi);
  };

  Vec3 operator&(const Vec3 &r1, const Vec3 &r2);
  double operator*(const Vec3 &r1, const Vec3 &r2);
  Vec3 operator*(const Vec3 &r1, double);
  Vec3 operator+(const Vec3 &r1, const Vec3 &r2);

  std::ostream& operator<<(std::ostream& os, const Vec3 &r);
}

#endif
