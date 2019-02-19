#include <QImage>

#include <complex>
#include <cstring>
#include <iostream>

#include "fftw3.h"
#include "plan_fftw.h"
#include "matrix_fftw.h"
#include "util.h"

using namespace std;
using namespace farn;


void work(){
  MatrixFFTW a=loadFromPicture("a.png", 0, 1);
  MatrixFFTW b=loadFromPicture("b.png", 0, 1);

  MatrixFFTW fftwSour(a.height(), a.width());
  MatrixFFTW fftwDest(a.height(), a.width());

  PlanFFTW forw(fftwSour, fftwDest, FFTW_FORWARD, FFTW_ESTIMATE);
  PlanFFTW backw(fftwSour, fftwDest, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftwSour.copyFrom(a);
  forw.execute();
  MatrixFFTW af(fftwDest);

  fftwSour.copyFrom(b);
  forw.execute();
  MatrixFFTW bf(fftwDest);

  MatrixFFTW cf(af.elementMultiply(bf));
  fftwSour.copyFrom(cf);
  backw.execute();

  saveAsPicture(fftwDest, "c.png");
}

int main(int argc, char* argv[]){
  try{
    work();
  }catch(string msg){
    cout << "error: " << msg << endl;
  }
  return 0;
}
