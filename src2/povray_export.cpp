#include "util.h"

#include "povray_export.h"

using namespace std;
using namespace farn;

PovrayExport::PovrayExport(string flnm, ConstructionMode mode):dest(flnm.c_str()){
  if (mode == HEADER) {
  dest << "#version 3.7;\n\
\n\
global_settings {assumed_gamma 1.0}\n\
\n\
camera{\n\
  location <10 * cos(2 * pi * clock), 2, 10 * sin(2 * pi * clock)>\n\
  rotate <0, 0, 0>\n\
  angle 30\n\
  look_at <0, 0, 0>\n\
}\n\
\n\
\n\
light_source{\n\
  <10, 3, 3>\n\
  color <1, 1, 1>\n\
}\n\
\n\
light_source{\n\
  <10, 3, -3>\n\
  color <1, 1, 1>\n\
}\n\
\n\
cylinder{<-0.3, 0, 0> <2.5, 0, 0> 0.01 pigment{color <1,0,0>}}\n\
cylinder{<0, 0, -0.3> <0, 0, 2.5> 0.01 pigment{color <0,1,0>}}\n\
cylinder{<0, -0.3, 0> <0, 2.5, 0> 0.01 pigment{color <0,0,1>}}\n\
\n";
  }
}

void PovrayExport::sphere(double x, double y, double z, double r){
  dest << "sphere { <" << x << ", " << z << ", " << y << "> " << r << " pigment{color <0.3, 0.3, 0.3>}}\n";
}

void PovrayExport::matrix(const farn::RotationMatrix& r){
  dest << "cylinder{<0, 0, 0> <" << r[0][0] << "," << r[0][2] << "," << r[0][1] << "> 0.03 pigment{color <1,0,0>}}\n";
  dest << "cylinder{<0, 0, 0> <" << r[1][0] << "," << r[1][2] << "," << r[1][1] << "> 0.03 pigment{color <0,1,0>}}\n";
  dest << "cylinder{<0, 0, 0> <" << r[2][0] << "," << r[2][2] << "," << r[2][1] << "> 0.03 pigment{color <0,0,1>}}\n";
}
