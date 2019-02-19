#ifndef POVRAY_EXPORT
#define POVRAY_EXPORT

#include <fstream>
#include <string>

#include "util.h"

enum ConstructionMode {HEADER, NO_HEADER};

class PovrayExport{
  std::ofstream dest;

 public:

  PovrayExport(std::string flnm, ConstructionMode mode = HEADER );
  void sphere(double x, double y, double z, double r=0.01);
  void matrix(const farn::RotationMatrix&);
};

#endif
