#include <iostream>

#include "matrix_fftw.h"
#include "spacial_matrix.h"
#include "storage.h"

using namespace farn;
using namespace std;

void
SpacialMatrix::printv(ostream& os) const {
  for(int p = 0; p < height(); ++p) {
    for(int q = 0; q < width(); ++q) {
      os << (*this)(p,q).getV() << " ";
    }
    os << endl;
  }
}

ostream& 
farn::operator<<(std::ostream& os, const SpacialMatrix& r){
  os << "h = " << r.h << " w = " << r.w << endl;
  
  for(int p = 0; p < r.height(); ++p){
    for(int q = 0; q < r.width(); ++q){
      os << "(" << p << "," << q << "): " << endl;
      os << r(p,q) << endl << endl;
    }
  }
  return os;
}

void
SpacialMatrix::fillSliceWithV(int t, Storage& storage) {
  for(int p = 0; p < height(); ++p) {
    for(int q = 0; q < width(); ++q) {
      storage(t, p, q) = (*this)(p, q).getV();
    }
  }
}

void
SpacialMatrix::fillMatrixWithV(MatrixFFTW& mat) {
  for(int p = 0; p < height(); ++p) {
    for(int q = 0; q < width(); ++q) {
      mat(p, q) = (*this)(p, q).getV();
    }
  }
}

