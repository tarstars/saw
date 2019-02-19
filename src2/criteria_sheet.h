#pragma once
#include "criteria.h"

namespace farn {
  class PlaneWave;
  class CompositeWave;

  class CriteriaSheet : public Criteria {
    int sheetNumber;
  public:
    CriteriaSheet(int);
    virtual bool predicat(const PlaneWave&, const CompositeWave&) const;
  };
}

