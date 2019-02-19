#ifndef STORAGE
#define STORAGE

#include <vector>
#include <complex>

namespace farn{

  class Mat3;
  class MatrixFFTW;
  class Vec3;

  class Storage{
    int h, d, w;
    std::vector<std::complex<double> > dat;
  public:

  Storage():h(0), d(0), w(0) {}
  Storage(int hh, int dd, int ww):h(hh), d(dd), w(ww), dat(hh * dd * ww){}

    int height() const{return h;}
    int depth()  const{return d;}
    int width()  const{return w;}
  
    std::complex<double>& operator()(int t, int p, int q) {return dat.at(t * d * w + p * w + q);}
    const std::complex<double>& operator()(int t, int p, int q)const {return dat.at(t * d * w + p * w + q);}

    bool operator==(const Storage&) const;

    MatrixFFTW sliceHW(int ind) const;
    MatrixFFTW sliceHD(int ind) const;
    MatrixFFTW sliceDW(int ind) const;
    MatrixFFTW accumulateAlongDirection(const Vec3&) const;

    void findNan(std::ostream&) const;

    double calcRealThreshold(double) const;
    void printShortModules(std::ostream&) const;
    void printJson(std::ostream&) const;
    void saveForDx(const std::string&);

    void save(const std::string&);
    void load(const std::string&);

    friend std::ostream& operator<<(std::ostream& os, const Storage& r);
  };
  std::ostream& operator<<(std::ostream& os, const Storage& r);
}



#endif
