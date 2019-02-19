#include "criteria.h"
#include "composite_wave.h"
#include "plane_wave.h"
#include "util.h"

#include <algorithm>
#include <iostream>

using namespace farn;
using namespace std;

bool simple_complex_compare(const complex<double> & left, const complex<double> & right) {
  return real(left) < real(right);
}

CompositeWave::CompositeWave(double xs1, double xs2, 
			     const MaterialTensor& tens, 
			     double rho, double omega, 
			     const Vec3c& unitForceReal, 
			     CalculationType ct): s1(xs1), s2(xs2), eliminateAsSaw(false) {

#ifdef DEBUG
  cerr << "composite wave s1 = " << s1 << " s2 = " << s2 << endl;
  cerr << "material tensor = " << tens << endl;
#endif

  PolyMatrix poly_mat = makePolyMatrix(tens, rho, s1, s2);
#ifdef DEBUG
  cerr << "poly matrix done " << endl;
#endif

  Poly pol = det(poly_mat);
#ifdef DEBUG
  cerr << "determinant done " << endl;
#endif

  Poly::RootVec roots = pol.all_roots();
#ifdef DEBUG
  cerr << "roots done " << endl;
  cerr << "all roots: ";
  for(int t = 0; t < (int)roots.size(); ++t) {
    cerr << roots[t] << " ";
  }
  cerr << endl;
#endif
  
  Vec3c unitForce(unitForceReal);
#ifdef DEBUG
  cerr << "unit force done " << endl;
#endif
  
  double eps = 1e-15;

  vector<complex<double> > s3s;

#ifdef DEBUG
  cerr << "ready to fill s3c " << endl;
  int real_roots_number = 0;
#endif

  for(int t = 0; t < int(roots.size()); ++t) {
    complex<double> val = roots[t];
#ifdef DEBUG
    if (abs(imag(val)) < eps) {
      ++real_roots_number;
    }
#endif
    
    if ((abs(imag(val)) < eps && rightPointing(s1, s2, real(val), poly_mat, tens, omega)) || imag(val) < -eps) {
      s3s.push_back(val);
    }
  }

#ifdef DEBUG
  cerr << "selected roots number =  " << real_roots_number << endl;
#endif
  
  sort(s3s.begin(), s3s.end(), simple_complex_compare);

  if (s3s.size() != 3){
    stringstream msg;

    msg << "CompositeWave: number of selected roots is not equal 3. Roots are: " << endl;
    for(int t = 0; t < (int)roots.size(); ++t)
      msg << roots[t] << endl;
    msg << endl;

    msg << "selected roots are:" << endl;
    for(int t = 0; t < (int)s3s.size(); ++t)
      msg << s3s[t] << endl;

    throw(msg.str());
  }
  
#ifdef DEBUG
  cerr << "\ts3 roots are " << s3s << endl;
#endif

  for(int t = 0; t < 3; ++t) {
    wav[t] = PlaneWave(s1, s2, s3s[t], poly_mat, tens, omega);
  }

  Vec3c vecs[3];
  for(int t = 0; t < 3; ++t) {
    if (ct == FixedForce) {
      vecs[t] = wav[t].getLastColumn();
    }
    if (ct == FixedDisplacement) {
      vecs[t] = wav[t].getQ();
    }
  }

#ifdef DEBUG
  cerr << "vecs:" << endl;
  cerr << "1 -> " << vecs[0] << endl;
  cerr << "2 -> " << vecs[1] << endl;
  cerr << "3 -> " << vecs[2] << endl;
#endif

  Mat3 cm  = columns2matrix(vecs[0], vecs[1], vecs[2]);
  Mat3 cm1 = columns2matrix(unitForce, vecs[1], vecs[2]);
  Mat3 cm2 = columns2matrix(vecs[0], unitForce, vecs[2]);
  Mat3 cm3 = columns2matrix(vecs[0], vecs[1], unitForce);
  
  complex<double> denum = cm.det();

  /*if (cm.isZeroDet()) {
    eliminateAsSaw = true;
    Ak[0] = Ak[1] = Ak[2] = 0;
    } else {*/
    Ak[0] = cm1.det() / denum;
    Ak[1] = cm2.det() / denum;
    Ak[2] = cm3.det() / denum;
    //}

#ifdef DEBUG
  cerr << "summ = " << (vecs[0] * Ak[0] + vecs[1] * Ak[1] + vecs[2] * Ak[2]) << endl;
  cerr << endl << endl;
#endif

}


std::ostream& farn::operator<<(std::ostream& os, const CompositeWave& r) {
  os << "CompositeWave:" << endl;
  for(int t = 0; t < 3; ++t) {
    os << "\tweight of this wave is " << r.Ak[t] << endl;
    os << r.wav[t] << endl;
  }

  return os;
}


void
CompositeWave::incrementStorage(Storage& stor, int p, int q) const {
  for(int t = 0; t < 3; ++t) {
    complex<double> val = A * Ak[t];

#ifdef DEBUG
    if (val != val) {
      throw(string("found nan in makeZShift"));
    }
#endif
 
    wav[t].incrementStorage(stor, val, p, q);
  }
}

void
CompositeWave::makeZShift(double dz) {
  for(int t = 0; t < 3; ++t) {
    complex<double> val = exp(complex<double>(0,1) * dz * wav[t].getKz());
    
#ifdef DEBUG
    if (val != val) {
      throw(string("found nan in makeZShift"));
    }
#endif
 
    Ak[t] *= val;
  }
}

void
CompositeWave::logState(ostream& os) const {
  os << "B_1 B_2 B_3 = " << Ak[0] << " " << Ak[1] << " " << Ak[2] << " ";
}

void
CompositeWave::printSumOfForces(ostream& os) const {
  Vec3c v0 = wav[0].getLastColumn();
  Vec3c v1 = wav[1].getLastColumn();
  Vec3c v2 = wav[2].getLastColumn();

  os << "vec0 = " << v0 << endl;
  os << "vec1 = " << v1 << endl;
  os << "vec2 = " << v2 << endl;
}

void
CompositeWave::filter(const Criteria& crit) {
  for(int p = 0; p < 3; p++) {
    if (!crit.predicat(wav[p], *this)) {
      Ak[p] = 0;
    }
  }
}

int
CompositeWave::getRealRootsNumber() const {
  int ret = 0;

  for(int t = 0; t < 3; t++) {
    if (wav[t].isHomo()) {
      ++ret;
    }
  }

  return ret;
}

double
CompositeWave::getTransmission() const {
  Vec3c ret;

  for(int t = 0; t < 3; ++t) {
    ret += wav[t].getQ() * Ak[t];
  }

  return ret.getSuperModule();
}

TransmissPoint 
CompositeWave::getTransmissPoint() const {
  TransmissPoint ret;

  ret.t = getTransmission();
  ret.s1 = s1;
  ret.s2 = s2;

  return ret;
}

int
CompositeWave::getEliminateCode() const {
  if (eliminateAsSaw) {
    return 1;
  }
  return 0;
}

void
CompositeWave::setAmplitude(const complex<double>& r) {
  A = r;
}

void
CompositeWave::setAmplitudeOne() {
  A = 1;
}

void
CompositeWave::printDxPosition(ostream& os, int ind) const {
  wav[ind].printDxPosition(os);
}

double
CompositeWave::displacementMagnitude(int ind) const {
  return abs(A * Ak[ind]);
}
