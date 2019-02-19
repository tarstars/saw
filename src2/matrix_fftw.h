#ifndef MATRIXFFTW
#define MATRIXFFTW

#include <complex>
#include <fftw3.h>


namespace farn{

  class MatrixFFTW{
    int h, w;
    fftw_complex *pDat;

    MatrixFFTW& operator=(const MatrixFFTW&);

  public:
    MatrixFFTW();
    MatrixFFTW(int hh, int ww);

    ~MatrixFFTW();

    MatrixFFTW(const MatrixFFTW& r);

    std::complex<double>& operator()(int p, int q);
    const std::complex<double>& operator()(int p, int q)const;
    double summ() const;

    void octaveRepr(std::ostream&)const;

    void copyFrom(const MatrixFFTW&);
    MatrixFFTW cleanClone()const;

    int width()const{return w;}
    int height()const{return h;}
    MatrixFFTW elementMultiply(const MatrixFFTW&)const;

    void printShortModules(std::ostream&)const;

    friend class PlanFFTW;
  };

  std::ostream& operator<<(std::ostream&, const MatrixFFTW& r);
}
#endif
