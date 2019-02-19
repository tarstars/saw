#include "matrix_fftw.h"
#include "mat3.h"
#include "povray_export.h"
#include "storage.h"
#include "vec3.h"

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace farn;

double nullifier(double x){
  if( abs(x) < 1e-110)
    return 0;
  return x;
}

ostream& 
farn::operator<<(ostream& os, const Storage& r){
  for(int t = 0; t < r.h; ++t){
    os << t << ")" << endl;
    for(int p = 0; p < r.d; ++p){
      for(int q = 0; q < r.w; ++q){
	complex<double> v = r(t,p,q);
	os << "(" << nullifier(real(v)) << "," << nullifier(imag(v))  << ") ";
      }
      os << endl;
    }
    os << endl;
  }
  return os;
}


MatrixFFTW 
Storage::sliceHW(int ind)const{
  MatrixFFTW ret(h, w);
  
  for(int p = 0; p < h; ++p) {
    for(int q = 0; q < w; ++q) {
      ret(p, q) = (*this)(p, ind, q);
    }
  }

  return ret;
}

MatrixFFTW 
Storage::sliceHD(int ind)const{
  MatrixFFTW ret(h, d);

  for(int p = 0; p < h; ++p){
    for(int q = 0; q < d; ++q){
      ret(p, q) = (*this)(p, q, ind);
    }
  }

  return ret;
}


MatrixFFTW 
Storage::sliceDW(int ind)const{
  MatrixFFTW ret(d, w);

  for(int p = 0; p < d; ++p){
    for(int q = 0; q < w; ++q){
      ret(p, q) = (*this)(ind, p, q);
    }
  }

  return ret;
}


void
Storage::printShortModules(ostream& os) const {
  os.precision(3);

  for(int p = 0; p < height(); ++p) {
    os << p << ")" << endl;
    for(int q = 0; q < depth(); ++q) {
      os << "\t";
      for(int r = 0; r < width(); ++r) {
	os << abs((*this)(p, q, r)) << " ";
      }
      os << endl;
    }
  }
}

void
Storage::printJson(ostream& os) const {
  os.precision(10);

  os << "[" << endl;
  for(int p = 0; p < height(); ++p) {
    os << " [";
    for(int q = 0; q < depth(); ++q) {
      os << "  [";
      for(int r = 0; r < width(); ++r) {
	complex<double> val = (*this)(p, q, r);
	os << "(" << real(val) << "+" << imag(val) << "j), ";
      }
      os << "]," << endl;
    }
    os << " ]," << endl;
  }
  os << "]" << endl;
}

void
Storage::findNan(ostream& os) const {
  for(int p = 0; p < height(); ++p) {
    for(int q = 0; q < depth(); ++q) {
      for(int r = 0; r < width(); ++r) {
	double re = real((*this)(p, q, r));
	double im = imag((*this)(p, q, r));

	if (re!=re || im != im) {
	  os << "nan at position (p, q, r) = " << p << ", " << q << ", " << r << endl;
	}
      }
    }
  }
}

void
Storage::saveForDx(const string& flnm_base) {
  string flnm_dx = flnm_base + ".dx";
  string flnm_data = flnm_base+ ".dat";

  ofstream destDx(flnm_dx);
  destDx << "object 1 class gridpositions counts " << height() << " " << depth() << " " << width() << endl;
  destDx << "origin 0 0 0" << endl;
  destDx << "delta 0 0 1" << endl;
  destDx << "delta 0 1 0" << endl;
  destDx << "delta 1 0 0" << endl;
  destDx << "attribute \"dep\" string \"positions\"" << endl;
  destDx << "#" << endl;
  destDx << "object 2 class gridconnections counts " << height() << " " << depth() << " " << width() << endl;
  destDx << "attribute \"element type\" string \"cubes\"" << endl;
  destDx << "attribute \"dep\" string \"connections\"" << endl;
  destDx << "attribute \"ref\" string \"positions\"" << endl;
  destDx << "#" << endl;
  destDx << "object 3 class array type double rank 0 items " << height() * depth() * width() << " lsb binary data file " << flnm_data << ",0" << endl;
  destDx << "#" << endl;
  destDx << "object \"my test\" class field" << endl;
  destDx << "component \"positions\" value 1" << endl;
  destDx << "component \"connections\" value 2" << endl;
  destDx << "component \"data\" value 3" << endl;
  destDx << "#" << endl;
  destDx << "end " << endl;

  ofstream destData(flnm_data);
  for(int p = 0; p < height(); ++p) {
    for(int q = 0; q < depth(); ++q) {
      for(int r = 0; r < width(); ++r) {
	double val = abs((*this)(p, q, r));
	destData.write((char*)&val, sizeof(val));
      }
    }
  }
}

void
Storage::save(const string& flnm) {
  cerr << "save storage" << endl;
  cerr << "\tdimensions = " << h << " " << d << " " << w << endl;
  ofstream dest(flnm.c_str());
  if (!dest) {
    cerr << "can not save storage in text format to path '" << flnm << "'" << endl;
    return;
  }

  dest.precision(10);

  dest << h << " " << d << " " << w << endl;
  for(int r = 0; r < h; r++) {
    for(int p = 0; p < d; p++) {
      for(int q = 0; q < w; q++) {
	dest << (*this)(r, p, q) << " ";
      }
    }
  }
}

void
Storage::load(const string& flnm) {
  cerr << "load storage" << endl;
  ifstream source(flnm.c_str());

  if (source >> h >> d >> w) {}
  else {
    cerr << "can not open file " << flnm << endl;
  }

  cerr << "\tdimensions = " << h << " " << d << " " << w << endl;
  
  (*this) = Storage(h, d, w);

  for(int r = 0; r < h; r++) {
    cerr << r << " ";
    for(int p = 0; p < d; p++) {
      for(int q = 0; q < w; q++) {
	source >> (*this)(r, p, q);
      }
    }
  }
}

namespace farn {
  double brezenDeviat(const Vec3&, const Vec3&);
}

MatrixFFTW
Storage::accumulateAlongDirection(const Vec3& dir) const {
  cerr << "accumulate along direction" << endl;
  MatrixFFTW ret(d, w);

  for(int destp = 0; destp < d; destp++) {
    for(int destq = 0; destq < w; destq++) {

      complex<double> accumVal = 0;
      double amount = 0;

      int r = 0, p = 0, q = 0;

      while(r < h) {
	bool firstTime = true;
	double deviat = 0;
	int rr = 0, pp = 0, qq = 0; 

	for(int dp = -1; dp < 2; dp++) {
	  for(int dq = -1; dq < 2; dq++) {
	    for(int dr = 0; dr < 2; dr++) {
	      if (!(dp == 0 && dq == 0 && dr == 0)) {
		int nr = r + dr, np = p + dp, nq = q + dq;
		double val = brezenDeviat(Vec3(nr, np, nq), dir);
		if (firstTime) {
		  firstTime = false;
		  deviat = val;
		  rr = nr;
		  pp = np; 
		  qq = nq;
		} else {
		  if (val < deviat) {
		    deviat = val;
		    pp = np;
		    qq = nq;
		    rr = nr;
		  }
		}
	      }
	    }
	  }
	}
    
	int ranged_p = (destp + p) % d;
	int ranged_q = (destq + q) % w;

	ranged_p = (ranged_p + d)%d;
	ranged_q = (ranged_q + w)%w;

	accumVal += (*this)(r, ranged_p, ranged_q);
	amount++;

	r = rr;
	p = pp; 
	q = qq;
      }

      ret(destp, destq) = accumVal / amount;
    }
  }

  return ret;
}

bool
Storage::operator==(const Storage& rt) const {
  if (rt.h != h) {
    return false;
  }

  if (rt.d != d) {
    return false;
  }

  if (rt.w != w) {
    return false;
  }

  for(int r = 0; r < 3; r++) {
    for(int p = 0; p < 3; p++) {
      for(int q = 0; q < 3; q++) {
	if (abs((*this)(r, p, q) - rt(r, p, q)) > 1e-8) {
	  cerr.precision(10);
	  cerr << "\t" << r << " " << p << " " << q << " " << (*this)(r, p, q) << " " << rt(r, p, q) << endl;
	  return false;
	}
      }
    }
  }

  return true;
}

double Storage::calcRealThreshold(double thr) const {
  double minVal = 0;
  double maxVal = 0; 

  bool firstTime = true;
  for(int r = 0; r < h; r++) {
    for(int p = 0; p < d; p++) {
      for(int q = 0; q < w; q++) {
	double val = abs((*this)(r, p, q));
	if (firstTime) {
	  minVal = val;
	  maxVal = val;
	  firstTime = false;
	} else {
	  maxVal = max(val, maxVal);
	  minVal = min(val, minVal);
	}
      }
    }
  }

  return minVal + (maxVal - minVal) * thr;
}



