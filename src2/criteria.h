#ifndef _CRITERIA_
#define _CRITERIA_

namespace farn {
  class PlaneWave;
  class CompositeWave;

class Criteria {
 public:
  virtual bool predicat(const PlaneWave&, const CompositeWave&) const = 0;
  virtual ~Criteria(){}
};
}

#endif
