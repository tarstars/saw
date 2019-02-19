#include "matrix_fftw.h"
#include "storage.h"
#include "wave_matrix.h"

using namespace farn;

WaveMatrix::WaveMatrix(int nn) : n(nn), dat(nn*nn) {

}

void
WaveMatrix::loadFftwMatrix(const MatrixFFTW& r) {
  if(r.width() != n || r.height() != n)
    throw(string("WaveMatrix::loadFftwMatrix: dimensions mismatch"));

  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      (*this)(p,q).setAmplitude(r(p,q));
    }
  }
}

std::ostream& 
farn::operator<<(std::ostream& os, const WaveMatrix &r) {
  os << "Wave Matrix: " << endl;
  for(int p = 0; p < r.n; ++p) {
    for(int q = 0; q < r.n; ++q) {
      os << "element(" << p << ", " << q << "): " << endl;
      os << r(p, q) << endl;
    }
    os << endl;
  }
  return os;
}

Storage
WaveMatrix::getStorage()const{
  Storage ret(PlaneWave::getDimensions(), n, n);

  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      (*this)(p,q).incrementStorage(ret, p, q);
    }
  }

  return ret;
}

void
WaveMatrix::makeZShift(double dz){
  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      (*this)(p,q).makeZShift(dz);
    }
  }
}

void
WaveMatrix::logState(ostream& os)const{
  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      os << "(p,q) = (" << p << "," << q << ") ";
      (*this)(p, q).logState(os);
      os << endl;
    }
  } 
}

void
WaveMatrix::filter(const Criteria& crit) {
  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      (*this)(p,q).filter(crit);
    }
  }
}

void
WaveMatrix::setOnes(void) {
  for(int p = 0; p < n; ++p){
    for(int q = 0; q < n; ++q){
      (*this)(p,q).setAmplitudeOne();
    }
  }
}

MatrixFFTW 
WaveMatrix::realRootsNumberMap()const {
  MatrixFFTW ret(n, n);

  for(int p = 0; p < n; ++p){
    for(int q = 0; q < n; ++q){
      ret(p,q) = (*this)(p,q).getRealRootsNumber();
    }
  }

  return ret;
}

MatrixFFTW
WaveMatrix::transmission() const {
  MatrixFFTW ret(n, n);

  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      ret(p, q) = (*this)(p, q).getTransmission();
    }
  }

  return ret;
}

vector<TransmissPoint>
WaveMatrix::gist() const {
  vector<TransmissPoint> ret;

  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      ret.push_back((*this)(p, q).getTransmissPoint());
    }
  }

  return ret;
}

MatrixFFTW
WaveMatrix::eliminateMap() const {
  MatrixFFTW ret(n, n);

  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      ret(p, q) = (*this)(p, q).getEliminateCode();
    }
  }

  return ret;
}

int fi(int p, int n) {
  return (p + n/2 + (n%2))%n;
}

void
WaveMatrix::dxExport(ostream& os) const {
  os << "object \"sheets pos1\" class array type float rank 1 shape 3 items " << n * n << " data follows" << endl;
  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      //os << p << " " << q << " " << 0 << endl;
      dat[p * n + q].printDxPosition(os, 0);
    }
  }
  os << endl;

  os << "object \"sheets pos2\" class array type float rank 1 shape 3 items " << n * n <<  " data follows" << endl;
  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      //os << p << " " << q << " " << 1 << endl;
      dat[p * n + q].printDxPosition(os, 1);
    }
  }
  os << endl;

  os << "object \"sheets pos3\" class array type float rank 1 shape 3 items " << n * n << " data follows" << endl;
  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      //os << p << " " << q << " " << 2 << endl;
      dat[p * n + q].printDxPosition(os, 2);
    }
  }
  os << endl;

  os << "object \"sheets connect\" class array type int rank 1 shape 3 items " << 2 * (n - 1) * (n - 1) << " data follows" << endl;
  for(int p = 0; p < n - 1; ++p) {
    for(int q = 0; q < n - 1; ++q) {
      os << (fi(p,n) * n + fi(q, n)) << " " << (fi(p,n) * n + fi(q + 1, n)) << " " << (fi((p + 1), n) * n + fi(q + 1, n)) << endl;
      os << (fi(p,n) * n + fi(q, n)) << " " << (fi((p + 1), n) * n + fi(q + 1,n)) << " " << (fi((p + 1),n) * n + fi(q,n)) << endl;
    }
  }
  os << "attribute \"element type\" string \"triangles\"" << endl;
  os << endl;

  for(int ind = 0; ind < 3; ++ind) {
    os << "object \"sheets data" << (ind + 1) << "\" class array type float rank 0 items " << n * n << " data follows" << endl;
    for(int p = 0; p < n; ++p) {
      for(int q = 0; q < n; ++q) {
	os << dat[p * n + q].displacementMagnitude(ind) << " ";
      }
    }
    os << endl << endl;
  }

  os << "object \"sheets field 1\" class field" << endl;
  os << "component \"positions\" \"sheets pos1\"" << endl;
  os << "component \"connections\" \"sheets connect\"" << endl;
  os << "component \"data\" \"sheets data1\"" << endl;
  os << endl << " " << endl;

  os << "object \"sheets field 2\" class field" << endl;
  os << "component \"positions\" \"sheets pos2\"" << endl;
  os << "component \"connections\" \"sheets connect\"" << endl;
  os << "component \"data\" \"sheets data2\"" << endl;
  os << endl << " " << endl;

  os << "object \"sheets field 3\" class field" << endl;
  os << "component \"positions\" \"sheets pos3\"" << endl;
  os << "component \"connections\" \"sheets connect\"" << endl;
  os << "component \"data\" \"sheets data3\"" << endl;
  os << endl << " " << endl;

  os << "object \"sheet group\" class group" << endl;
  os << "member \"sf 1\" \"sheets field 1\"" << endl;
  os << "member \"sf 2\" \"sheets field 2\"" << endl;
  os << "member \"sf 3\" \"sheets field 3\"" << endl;
  os << " " << endl;
}
