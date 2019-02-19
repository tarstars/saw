#pragma once

#include "mat3.h"
#include "vec3.h"
#include "storage.h"

namespace farn {
  class VolumeData {
  public:
    Storage vd;
    Mat3 cs;
    Vec3 geometry;
    Vec3 origin;
  
  VolumeData(const Storage& xvd, 
	     const Mat3& xcs, 
	     const Vec3& xgeometry, 
	     const Vec3& xorigin
	     ) : vd(xvd), cs(xcs), geometry(xgeometry), origin(xorigin) {}

  };
}
