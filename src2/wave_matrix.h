#ifndef WAVE_MATRIX
#define WAVE_MATRIX

#include <iostream>
#include <vector>

#include "util.h"
#include "composite_wave.h"

namespace farn{
  
  //  class Eliminate;
  class MatrixFFTW;
  class Storage;

  class WaveMatrix{
    int n;
    std::vector<CompositeWave> dat;
  public:
    WaveMatrix(int);
    CompositeWave& operator()(int p, int q){return dat[p * n + q];}
    const CompositeWave& operator()(int p, int q)const{return dat[p * n + q];}
    void loadFftwMatrix(const MatrixFFTW&);
    Storage getStorage() const;
    void logState(std::ostream&) const;
    void makeZShift(double);
    int dimension()const{return n;}
    void filter(const Criteria&);
    void setOnes(void);
    MatrixFFTW realRootsNumberMap() const;
    MatrixFFTW transmission() const;
    MatrixFFTW eliminateMap() const;
    void dxExport(std::ostream&) const;
    std::vector<TransmissPoint> gist() const;

    friend std::ostream& operator<<(std::ostream&, const WaveMatrix&);
  };

  std::ostream& operator<<(std::ostream&, const WaveMatrix&);
}

#endif
