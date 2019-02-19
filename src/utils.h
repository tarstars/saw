#ifndef _UTILS_
#define _UTILS_

#include <vector>
#include <iostream>

#include "eigelement.h"

class Mat3;
class Vec3;
class Coeff;
class CPoly3;
class Poly;

Mat3 n2crist(const Coeff& m, const Vec3 &n);
CPoly3 mat2poly(const Mat3& m);
std::vector<CEigElement> vel_polar(const Vec3 &n);
void one_stage(double go, std::ostream& dest3d);
void generate3d();
poly_type::value_type det_val(long double  n1, long double n2, long double g);
Poly n1n2g_to_poly(long double n1, long double n2, long double g);


#endif
