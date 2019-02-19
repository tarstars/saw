#include "one_dir_result.h"

#include "vec3.h"
#include "vec3c.h"

#include <complex>

using namespace std;
using namespace farn;

OneDirResult::OneDirResult(const Vec3& dir, const MaterialTensor& mt, double rho) : calcDir(dir) {
  #ifdef DEBUG
  cerr << "\tOneDirResult: dir = " << dir << endl;
  #endif

  PolyMatrix cm = makeChristoffelPolyMatrix(mt, dir[0], dir[1], dir[2]);
  Poly pm = det(cm);
  Poly::RootVec ar = pm.all_roots();
  vector<double> gammas(ar.size());
  transform(ar.begin(), ar.end(), gammas.begin(), [](complex<double>& x) -> double {return abs(x);});
  sort(gammas.begin(), gammas.end());
  
  for(auto it = gammas.begin(); it != gammas.end(); ++it) {
    ComplexMatrix zdm = evaluatePolyMatrix(cm, *it);
    Vec3c q = calcPol(zdm);

    dat.push_back(VelocityElement(sqrt(*it / rho),
				  dir, 
				  Vec3(real(q[0]), real(q[1]), real(q[2]))));
  }
}

ostream&
farn::operator<<(ostream& os, const OneDirResult& r) {
  os << "\tone dir result:" << endl;
  os << "\t\t" << "dir = " << r.calcDir << endl;
  for(VelocityElement ve : r.dat) {
    os << ve << endl;
  }

  return os;
}

void
OneDirResult::exportPosition(ostream& os, int ind) const {
  dat[ind].exportPosition(os);
}

void
OneDirResult::exportPolarization(ostream& os, int ind) const {
  dat[ind].exportPolarization(os);
}
