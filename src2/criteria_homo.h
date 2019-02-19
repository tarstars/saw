#ifndef _CRITERIA_HOMO_
#define _CRITERIA_HOMO_

#include "criteria.h"

namespace farn {
  class PlaneWave;
  class CompositeWave;

class CriteriaHomo : public Criteria {
 public:
  bool predicat(const PlaneWave&, const CompositeWave&) const;
};

}
#endif
