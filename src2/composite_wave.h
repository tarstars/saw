#ifndef COMPOSITE_WAVE
#define COMPOSITE_WAVE

#include <iostream>

#include "util.h"
#include "plane_wave.h"

namespace farn{
  class Storage;
  class Criteria;

  class CompositeWave{
    std::complex<double> A;
    std::complex<double> Ak[3];
    PlaneWave wav[3];
    double s1, s2;
    bool eliminateAsSaw;

  public:
    CompositeWave(){}

    CompositeWave(double s1, double s2, 
		  const MaterialTensor&tens, 
		  double rho, double omega, 
		  const Vec3c& unitForce, 
		  CalculationType ct);

    void incrementStorage(Storage& stor, int p, int q) const;
    void logState(std::ostream&) const;
    void makeZShift(double dz);
    void printSumOfForces(std::ostream&) const;
    void filter(const Criteria&);
    void setAmplitude(const std::complex<double>&);
    int getRealRootsNumber() const;
    double getTransmission() const;
    TransmissPoint getTransmissPoint() const;
    int getEliminateCode() const;
    void setAmplitudeOne();
    void printDxPosition(std::ostream&, int ind) const;
    double displacementMagnitude(int) const;

    friend class CriteriaSheet;
    friend std::ostream& operator<<(std::ostream& os, const CompositeWave& r);
  };

  std::ostream& operator<<(std::ostream& os, const CompositeWave& r);
}

#endif
