#include "criteria.h"
#include "criteria_homo.h"

#include "plane_wave.h"

using namespace farn;

bool
CriteriaHomo::predicat(const farn::PlaneWave& pw, const farn::CompositeWave&) const {
  return pw.isHomo();
}
