#include "piezo_tensor.h"

using namespace std;

namespace farn{

  PiezoTensor::PiezoTensor(){
    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q)
	for(int r = 0; r < 3; ++r)
	  dat[p][q][r] = 0;
  }
  
  double& PiezoTensor::operator()(int p, int q, int r){return dat[p][q][r];}
  const double& PiezoTensor::operator()(int p, int q, int r)const{return dat[p][q][r];}

  std::ostream& operator<<(std::ostream& os, const PiezoTensor& ten){
    for(int p = 0; p < 3; ++p){
      for(int q = 0; q < 3; ++q){
	for(int r = 0; r < 3; ++r)
	  os << ten(p,q,r) << " ";
	os << endl;
      }
      os << endl;
    }
    return os;
  }

}
