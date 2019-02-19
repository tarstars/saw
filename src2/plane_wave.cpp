#include "plane_wave.h"
#include "storage.h"

using namespace farn;
using namespace std;

PlaneWave::PlaneWave(double s1, 
		     double s2, 
		     const complex<double>& s3, 
		     const PolyMatrix& poly_mat, 
		     const MaterialTensor& tens, 
		     double omega):
  slow(s1, s2, s3), 
  k(slow * omega) 
{
  ComplexMatrix curMat = evaluatePolyMatrix(poly_mat, s3);

  q = calcPol(curMat);
  S = strainFromKQ(k, q);
  T = stressFromCS(tens, S);

  #ifdef DEBUG
  cerr << "PlaneWave" << endl;
  cerr << "s1, s2, s3 = " << s1 << ", " << s2 << ", " << s3 << endl;
  cerr << "q = " << q << endl;
  cerr << "curMat * q = " << curMat * q << endl << endl;
  #endif
}

Vec3c
PlaneWave::getLastColumn() const {
  return Vec3c(T[0][2], T[1][2], T[2][2]);
}

complex<double>
PlaneWave::getKz()const {
  return k[2];
}

std::ostream& farn::operator << (ostream& os, const PlaneWave& r){
  os << "\tplane wave:" << endl;
  os << "\t\tslow = " << r.slow << endl;
  os << "\t\tq = " << r.q << endl;
  os << "\t\tk = " << r.k << endl;
  os << "\t\tS = " << r.S << endl;
  os << "\t\tT = " << r.T << endl;

  return os;
}

int
PlaneWave::getDimensions(){
  return 15;
}

void 
PlaneWave::incrementStorage(Storage& stor, const std::complex<double> & ampl, int p, int q) const {
  stor( 0, p, q) += ampl * this -> q.x();
  stor( 1, p, q) += ampl * this -> q.y();
  stor( 2, p, q) += ampl * this -> q.z();

  stor( 3, p, q) += (S[0][0]) * ampl;
  stor( 4, p, q) += (S[1][1]) * ampl;
  stor( 5, p, q) += (S[2][2]) * ampl;
  stor( 6, p, q) += (S[1][2]) * ampl;
  stor( 7, p, q) += (S[0][2]) * ampl;
  stor( 8, p, q) += (S[0][1]) * ampl;

  stor( 9, p, q) += (T[0][0]) * ampl;
  stor(10, p, q) += (T[1][1]) * ampl;
  stor(11, p, q) += (T[2][2]) * ampl;
  stor(12, p, q) += (T[1][2]) * ampl;
  stor(13, p, q) += (T[0][2]) * ampl;
  stor(14, p, q) += (T[0][1]) * ampl;
}

Vec3c
PlaneWave::getPointing() const {
  #ifdef DEBUG
  cerr << "T = " << endl << T << endl;
  cerr << "q = " << endl << q << endl;
  cerr << "Pointing = " << endl << T * q << endl << endl;
  #endif

  return T * q; 
}

bool
PlaneWave::isHomo() const {
  //cerr << slow.z() << " ";
  //bool flag = ;
  //cerr << (flag ? "homo" : "nonhomo") << endl;
  return abs(imag(slow.z())) < 1e-6;
}

Vec3c
PlaneWave::getQ() const {
  return q;
}

void
PlaneWave::printDxPosition(ostream& os) const {
  double m = 1;
  double x = m*real(slow[0]);
  double y = m*real(slow[1]);
  double val = -real(slow[2]);

  if (abs(imag(slow[2])) > 1e-5) {
    val = 0;
  }

  os << x  << " " << y << " " << val << endl; 
}
