//#pragma once

#include "criteria.h"
#include "criteria_nonhomo.h"

#include "plane_wave.h"

using namespace farn;

bool
CriteriaNonHomo::predicat(const farn::PlaneWave& pw, const farn::CompositeWave&) const {
  return !pw.isHomo();
}
