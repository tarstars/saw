#include "matrix_fftw.h"
#include "plan_fftw.h"
#include <fftw3.h>

using namespace farn;

PlanFFTW::PlanFFTW(MatrixFFTW& in, MatrixFFTW& out, int dir, int flag){
  plan = fftw_plan_dft_2d(in.h, in.w, in.pDat, out.pDat, dir, flag);
}

PlanFFTW::~PlanFFTW(){
  fftw_destroy_plan(plan);
}

void
PlanFFTW::execute()const{
  fftw_execute(plan);
}

