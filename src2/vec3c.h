#ifndef _VEC3C_
#define _VEC3C_

#include <iostream>
#include <vector>
#include <complex>

namespace farn{

  class Mat3;
  class Vec3;

  class Vec3c
  {
  public:
    std::vector<std::complex<double> > dat;
  
    explicit Vec3c(std::complex<double> ax = 0, std::complex<double> ay = 0, std::complex<double> az = 0): 
    dat(3)
      {
	dat[0] = ax; 
	dat[1] = ay; 
	dat[2] = az;
      }

    explicit Vec3c(const Vec3&);

    double getSuperModule()const;

    std::complex<double>& x(){return dat[0];}
    std::complex<double>& y(){return dat[1];}
    std::complex<double>& z(){return dat[2];}

    const std::complex<double>& x()const{return dat[0];}
    const std::complex<double>& y()const{return dat[1];}
    const std::complex<double>& z()const{return dat[2];}

    Vec3c operator*(double r) const {return Vec3c(r * dat[0], r * dat[1], r * dat[2]);}
    Vec3c operator*(const std::complex<double>& r) const {return Vec3c(r * dat[0], r * dat[1], r * dat[2]);}
    Vec3c operator+(const Vec3c& r) const {return Vec3c(dat[0] + r.dat[0], dat[1] + r.dat[1], dat[2] + r.dat[2]);}
    Vec3c& operator+=(const Vec3c& r);
    Vec3c operator-(const Vec3c& r) const {
      return Vec3c(dat[0] - r.dat[0], dat[1] - r.dat[1], dat[2] - r.dat[2]);
    }

    double abs() const;
    Vec3c norm() const {return (*this) * (1.0 / abs());}
    void normalize();

    std::complex<double>& operator[](int ind) {return dat[ind];}
    const std::complex<double>& operator[](int ind) const {return dat[ind];}

    Vec3c operator&(const Vec3c& r)const;
  };

  std::ostream& operator<<(std::ostream& os, const Vec3c &r);
  Vec3c operator*(const Vec3c &, const Mat3&);
  
}

#endif
