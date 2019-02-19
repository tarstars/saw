#pragma once

#include "one_dir_result.h"

#include "util.h"

#include <iostream>
#include <vector>
/* scheme of array indexing
  0 - north polar cap
  1..n                -  p = 1, q = 0..(n-1) - first row below polar cap
  ....
  n*(n-2)+1..n*(n-1)  -  p = (n-1), q = 0..(n-1)
  n*(n-1)+1
*/

namespace farn {
  class DxSurfaceCalculator {
    int n;
    std::vector<OneDirResult> dat;
  public:
    DxSurfaceCalculator(int point_number, const MaterialTensor&, double rho);
    void save(const char*, int ind);
    friend std::ostream& operator<<(std::ostream&, const DxSurfaceCalculator&);
  };

}
