#pragma once

#include "criteria.h"

namespace farn {
  class PlaneWave;

  class CriteriaAll:public Criteria {
  public:
    virtual bool NeedEliminate(const PlaneWave&) const {return true;}
  };
}
