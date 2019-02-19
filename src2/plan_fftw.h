#ifndef PLAN_FFTW
#define PLAN_FFTW

#include <fftw3.h>

namespace farn{
  class MatrixFFTW;

  class PlanFFTW{
    fftw_plan plan;
    PlanFFTW(const PlanFFTW&);
    PlanFFTW& operator=(const PlanFFTW&);
  public:
    PlanFFTW(MatrixFFTW& in, MatrixFFTW& out, int dir, int flag);
    ~PlanFFTW();
    void execute()const;
  };

}

#endif
