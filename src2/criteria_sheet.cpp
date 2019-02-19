#include "composite_wave.h"
#include "criteria_sheet.h"

#include <vector>
#include <algorithm>

using namespace std;
using namespace farn;

CriteriaSheet::CriteriaSheet(int n) : sheetNumber(n) {}

bool
CriteriaSheet::predicat(const PlaneWave& cand, const CompositeWave& cw) const {
  vector<double> dat(3);
  dat[0] = real(cw.wav[0].slow[2]);
  dat[1] = real(cw.wav[1].slow[2]);
  dat[2] = real(cw.wav[2].slow[2]);
  
  sort(dat.begin(), dat.end());

  #ifdef DEBUG
  cerr << "predicat: ind = " << sheetNumber << " {" << dat[0] << ", " << dat[1] << ", " << dat[2] << "} " << " val = " << real(cand.slow[2]) << endl;
  #endif
  
  return abs(real(cand.slow[2]) - dat[sheetNumber]) < 1e-5;
}
