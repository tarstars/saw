#ifndef PIEZOTENSOR
#define PIEZOTENSOR

#include <iostream>

namespace farn{

  class PiezoTensor{
    double dat[3][3][3];
  public:

    PiezoTensor();

    double& operator()(int, int, int);
    const double& operator()(int, int, int)const;
  };

  std::ostream& operator<<(std::ostream& os, const PiezoTensor&);
}


#endif
