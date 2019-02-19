#ifndef _CRITERIA_NONHOMO_
#define _CRITERIA_NONHOMO_

#include "criteria.h"

namespace farn {
  class PlaneWave;
  class CompositeWave;

class CriteriaNonHomo : public Criteria {
 public:
  bool predicat(const PlaneWave&, const CompositeWave&) const;
};

}
#endif
