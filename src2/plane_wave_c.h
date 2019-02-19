#ifndef PLANE_WAVEC
#define PLANE_WAVEC

#include "vec3c.h"
#include "mat3.h"
#include "util.h"

namespace farn{

  class Storage;

  class PlaneWaveC{
    farn::Vec3c q;
    Mat3 S, T;
  public:
    PlaneWaveC(){}

    void setTestValue(int p, int q);

    void takeFromStorage(const Storage&, int, int);
    double getV()const;

    friend std::ostream& operator << (std::ostream& os, const PlaneWaveC& r);
  };

  std::ostream& operator << (std::ostream& os, const PlaneWaveC& r);
}

#endif
