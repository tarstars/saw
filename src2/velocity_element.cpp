#include "velocity_element.h"

#include "vec3.h"

#include <iostream>

using namespace farn;
using namespace std;

VelocityElement::VelocityElement(double vv, const Vec3& xdir, const Vec3& qq) : vel(vv), slow(1/vv), q(qq), dir(xdir) {
}

bool
VelocityElement::operator<(const VelocityElement& r) const {
  return slow < r.slow;
}

ostream& farn::operator<<(std::ostream& os, VelocityElement const & r) {
  return os << "\t\t\tvel = " << r.vel << " slow = " << r.slow << " q = " << r.q;
}

void
VelocityElement::exportPosition(ostream& os) const {
 os << dir[0] * slow << " " << dir[1] * slow << " " << dir[2] * slow << endl;
}

void
VelocityElement::exportPolarization(ostream& os) const {
  os << real(q[0]) << " " << real(q[1]) << " " << real(q[2]) << endl;
}
