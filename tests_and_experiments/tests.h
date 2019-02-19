#ifndef _TESTS_
#define _TESTS_

#include <vector>
#include <string>
#include "util.h"

void rotate_angle_explore();
void rotate_explore();
farn::RotationMatrix farn::genChungMatrix(double alpha);
std::vector<farn::RotationMatrix> genMatSeq();
void generate_teo2_xy_slice() ;
void test_composite_wave();
void make_map(const std::string& flnm, /*const std::string&,*/ const farn::RotationMatrix& rm);
void make_map_mult();
void test_poly();
void test_povray_export();
void test_mutant_christoffel_chung_matrix();
void test_mutant_christoffel_poly_matrix();
void test_mutant_christoffel_rotate();
void test_mutant_christoffel_shift();
void test_mutant_christoffel_shift_rotate();
void testZSliceNormal(const char* flnm);
void testZSliceMutant(const char* flnm);



#endif
