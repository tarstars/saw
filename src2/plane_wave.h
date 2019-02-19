#ifndef PLANE_WAVE
#define PLANE_WAVE

#include "vec3.h"
#include "vec3c.h"
#include "mat3.h"
#include "util.h"

namespace farn{

  class Storage;

  class PlaneWave{
    Vec3c slow, q, k;
    Mat3 S, T;

  public:
    PlaneWave(){}
    PlaneWave(double s1, 
	      double s2, 
	      const std::complex<double>& s3, 
	      const PolyMatrix& poly_mat, 
	      const MaterialTensor& tens, 
	      double omega);
    
    Vec3c getLastColumn() const;

    static int getDimensions();
    void incrementStorage(Storage& stor, 
			  const std::complex<double> & ampl, 
			  int p, 
			  int q) const;

    std::complex<double> getKz()const;

    Vec3c getPointing() const;
    bool  isHomo() const;
    Vec3c getQ() const;
    void  printDxPosition(ostream&) const;

    friend class CriteriaSheet;
    friend std::ostream& operator << (std::ostream& os, const PlaneWave& r);
  };

  std::ostream& operator << (std::ostream& os, const PlaneWave& r);
}

#endif
