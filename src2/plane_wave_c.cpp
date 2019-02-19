#include <iostream>

#include "plane_wave_c.h"
#include "storage.h"

using namespace farn;
using namespace std;

ostream& 
farn::operator << (std::ostream& os, const PlaneWaveC& r) {
  os << "q = " << r.q << endl;
  os << "S = " << r.S << endl;
  os << "T = " << r.T << endl;
  return os;
}

void 
PlaneWaveC::takeFromStorage(const Storage& stor, int p, int q) {
  this -> q.x() = stor( 0, p, q);
  this -> q.y() = stor( 1, p, q);
  this -> q.z() = stor( 2, p, q);

  S[0][0] = stor( 3, p, q);
  S[1][1] = stor( 4, p, q);
  S[2][2] = stor( 5, p, q);
  S[1][2] = stor( 6, p, q);
  S[0][2] = stor( 7, p, q);
  S[0][1] = stor( 8, p, q);

  T[0][0] = stor( 9, p, q);
  T[1][1] = stor(10, p, q);
  T[2][2] = stor(11, p, q);
  T[1][2] = stor(12, p, q);
  T[0][2] = stor(13, p, q);
  T[0][1] = stor(14, p, q);
}

double
PlaneWaveC::getV() const {
  return q.getSuperModule();
}

void
PlaneWaveC::setTestValue(int p, int q) {
  this->q.x() = 10 * p + q;
}
