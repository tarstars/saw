#ifndef SPACIAL_MATRIX
#define SPACIAL_MATRIX

#include <vector>

#include "plane_wave_c.h"

namespace farn{

  class MatrixFFTW;

  class SpacialMatrix{
    int h, w;
    std::vector<PlaneWaveC> dat;

  public:

  SpacialMatrix():h(0), w(0){}
  SpacialMatrix(int hh, int ww):h(hh), w(ww), dat(hh * ww){}

    PlaneWaveC& operator()(int p, int q){return dat[p * w + q];}
    const PlaneWaveC& operator()(int p, int q)const{return  dat[p * w + q];}
    
    int height()const{return h;}
    int width()const{return w;}


    void printv(std::ostream&)const;
    void fillSliceWithV(int t, Storage&);
    void fillMatrixWithV(MatrixFFTW&);

    friend std::ostream& operator<<(std::ostream&, const SpacialMatrix&);
  };

  std::ostream& operator<<(std::ostream&, const SpacialMatrix&);
}

#endif
