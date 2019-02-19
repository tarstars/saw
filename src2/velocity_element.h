#pragma once

#include "vec3.h"

#include <iostream>

namespace farn {
  //class Vec3;

  class VelocityElement {
    double vel, slow;
    Vec3 q, dir;

  public:
    VelocityElement(double v, const Vec3& dir,  const Vec3& q);
    bool operator<(const VelocityElement&r)const;
    friend std::ostream& operator<<(std::ostream&, VelocityElement const & r);
    void exportPosition(std::ostream&) const;
    void exportPolarization(std::ostream&) const;
  };
}
