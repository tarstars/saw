#pragma once

#include "criteria.h"

namespace farn {
  class PlaneWave;

  class CriteriaNone : public Criteria {
  public:
    virtual bool needEliminate(const PlaneWave&) const {return false;}
  };
}
