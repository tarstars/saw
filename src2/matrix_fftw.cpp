#include <cstring>

#include "matrix_fftw.h"

using namespace std;
using namespace farn;

#define CHECK_BOUNDS

MatrixFFTW::MatrixFFTW(){}

MatrixFFTW::MatrixFFTW(int hh, int ww):h(hh), 
				       w(ww), 
				       pDat( (fftw_complex*) fftw_malloc(hh * ww * sizeof(fftw_complex)))
{}

MatrixFFTW::~MatrixFFTW(){
  fftw_free(pDat);
}

MatrixFFTW::MatrixFFTW(const MatrixFFTW& r){
  //fftw_free(pDat);

  h = r.h;
  w = r.w;
  pDat = (fftw_complex*) fftw_malloc(h * w * sizeof(fftw_complex));

  memcpy(pDat, r.pDat, h * w * sizeof(fftw_complex));
}

void
MatrixFFTW::copyFrom(const MatrixFFTW& r){
  if(w!=r.w || h != r.h){
    throw(string("dimensions is not match in MatrixFFTW::copy"));
  }

  for(int p = 0; p < h; ++p)
    for(int q = 0; q < w; ++q)
      (*this)(p,q) = r(p,q);
}

std::complex<double>& 
MatrixFFTW::operator()(int p, int q){
  #ifdef CHECK_BOUNDS
  if(p < 0) 
    throw(string("matrixFftw::operator() p is negative"));

  if(p > h) 
    throw(string("matrixFftw::operator() p is greater then h"));

  if(q < 0) 
    throw(string("matrixFftw::operator() q is negative"));

  if(q > w) 
    throw(string("matrixFftw::operator() q is greater than w"));
  #endif

  return (std::complex<double>&)pDat[p * w + q];
}

const std::complex<double>& 
MatrixFFTW::operator()(int p, int q)const{
  #ifdef CHECK_BOUNDS
  if(p < 0) 
    throw(string("matrixFftw::operator() p is negative"));

  if(p > h) 
    throw(string("matrixFftw::operator() p is greater then h"));

  if(q < 0) 
    throw(string("matrixFftw::operator() q is negative"));

  if(q > w) 
    throw(string("matrixFftw::operator() q is greater than w"));
  #endif
  
  return (std::complex<double>&)pDat[p * w + q];
}

MatrixFFTW
MatrixFFTW::elementMultiply(const MatrixFFTW& r)const{
  if(h != r.h || w != r.w){
    throw(string("MatrixFFTW::elementMultiple: matrix dimensions is not match"));
  }

  MatrixFFTW ret(h, w);

  for(int p = 0; p < h; ++p)
    for(int q = 0; q < w; ++q)
      ret(p,q) = (*this)(p,q) * r(p,q);

  return ret;
}

MatrixFFTW
MatrixFFTW::cleanClone()const{
  return MatrixFFTW(h,w);
}

std::ostream& 
farn::operator<<(std::ostream& os, const MatrixFFTW& r){
  int w = r.width(); 
  int h = r.height();

  for(int p = 0; p < h; ++p){
    for(int q = 0; q < w; ++q){
      os << real(r(p,q)) << " ";
    }
    os << endl;
  }

  return os;
}

void
MatrixFFTW::octaveRepr(ostream& os)const{
  os << "[";
  for(int p = 0; p < h; ++p){
    for(int q = 0; q < w; ++q){
      complex<double> val = (*this)(p,q);
      os << real(val) << " + " << imag(val) << "i ";
    }
    os << ";";
  }
  os << "]";
}

void
MatrixFFTW::printShortModules(ostream& os) const {
  os.precision(3);
  for(int p = 0; p < h; ++p) {

  }
}

double
MatrixFFTW::summ() const {
  double ret = 0;

  for(int p = 0; p < h; ++p) {
    for(int q = 0; q < w; ++q) {
      ret += abs((*this)(p, q));
    }
  }

  return ret;
}
