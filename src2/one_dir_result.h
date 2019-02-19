#pragma once

#include "util.h"
#include "velocity_element.h"

#include <iostream>
#include <vector>

namespace farn {
  class Vec3;

  class OneDirResult {
    Vec3 calcDir;
    std::vector<VelocityElement> dat;
  public:
    OneDirResult(){}
    OneDirResult(const Vec3&, const MaterialTensor&, double);
    void exportPosition(std::ostream&, int ind) const;
    void exportPolarization(std::ostream&, int ind) const;
    friend std::ostream& operator<<(std::ostream&, const OneDirResult&);
  }; 
}
