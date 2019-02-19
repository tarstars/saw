#ifndef HEIGHT_CALCULATOR
#define HEIGHT_CALCULATOR

#include <vector>
#include <utility>
#include "util.h"


class HeightCalculator{
  farn::MaterialTensor ten;

 public:

  HeightCalculator();
  double getHeight(double, double, int);
  std::vector<std::pair<double, double> > getAllHeights(double, double, double);
};

#endif
