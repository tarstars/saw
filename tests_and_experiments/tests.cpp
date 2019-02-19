#include <fstream>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_errno.h>

#include <QColor>
#include <QImage>

#include "composite_wave.h"
#include "povray_export.h"
#include "tests.h"
#include "util.h"
#include "vec3.h"
#include "povray_export.h"


using namespace std;
//using namespace Eigen;
using namespace farn;


void rotate_angle_explore(){
  ofstream dest("phi_10.txt");
  double sq2 = 1/sqrt(2);

  double phi = 10. / 180 * M_PI;

  RotationMatrix vert = combineInMatrix(Vec3(  0,    0, 1),
					Vec3(sq2, -sq2, 0), 
					Vec3(sq2,  sq2, 0) );

  RotationMatrix ident = combineInMatrix(Vec3(1, 0, 0), 
					 Vec3(0, 1, 0), 
					 Vec3(0, 0, 1));

  RotationMatrix toaxis = combineInMatrix(Vec3( sq2, sq2, 0), 
					  Vec3(-sq2, sq2, 0), 
					  Vec3(   0,   0, 1));

  RotationMatrix custRot = combineInMatrix(Vec3( cos(phi), sin(phi), 0),
					   Vec3(-sin(phi), cos(phi), 0), 
					   Vec3(        0,        0, 1) );

  RotationMatrix inUse = custRot;

  cout << "rotation matrix:" << endl << inUse << endl << endl;

  MaterialTensor  tt = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  MaterialTensor t = rotateTensor(tt, inUse);

  //for(double theta = 0; theta < M_PI; theta += 0.1)
  {
    double theta = M_PI / 2;
    for(double phi = 0; phi < 2 * M_PI; phi += 0.01)
      {
	Vec3 n(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
	Poly::RootVec roots = det(makeChristoffelPolyMatrix(t, n)).all_roots();
	for(int t = 0; t < 3; ++t) 
	  {
	    double s = 1/(sqrt(real(roots[t]) / 5.96e3));
	    dest << n[0] * s << " " << n[1] * s << /*" " << n[2] * s <<*/ endl;
	  }
      }
  }
}

void rotate_explore(){
  RotationMatrix m = genChungMatrix(10. / 180 * M_PI);
  RotationMatrix mt = transpose(m);
  cout << "rotation matrix: " << endl << m << endl;
  
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  MaterialTensor tt = rotateTensor(t, m);

  Vec3 v(0, 0, 1);
  cout << "christoffel matrix" << endl << makeChristoffelMatrix(t, v) << endl;

  PolyMatrix pm =  makeChristoffelPolyMatrix(t, v);
  cout << "christoffel poly matrix "  << endl << pm << endl;

  Vec3 v_rot = mt * v;
  cout << "rotated vector = " << v_rot << endl;

  PolyMatrix pm_r = makeChristoffelPolyMatrix(t, v_rot);
  cout << "rotated poly matrix" << endl << pm_r << endl;

  PolyMatrix pm_ten_rot = makeChristoffelPolyMatrix(tt, v);
  cout << "poly matrix with rotated tensor" << endl << pm_ten_rot << endl;

  cout << "roots of rotation vector" << endl << det(pm_r).all_roots() << endl;
  cout << "roots of rotation tensor" << endl << det(pm_ten_rot).all_roots() << endl;
  
}

// RotationMatrix genChungMatrix(double alpha){

//   boost::array<RotationMatrix::index, 2> shape = {{3,3}};
//   RotationMatrix mat(shape);
  
//   double ca = cos(alpha);
//   double sa = sin(alpha);
//   double sq = sqrt(2);

//   return combineInMatrix( Vec3(sa / sq, sa / sq, ca),
// 			    Vec3(-1 / sq, 1 / sq,  0),
// 			    Vec3(-ca / sq, -ca / sq, sa));  
// }

vector<RotationMatrix> genMatSeq(){
  vector<RotationMatrix> ret;
  //cout << "genMatSeq" << endl;

  int n = 100;

  for(int t = 0; t < n; ++t){
    double theta = t * (M_PI / 2) / n;
    Vec3 e3(sin(theta), 0, cos(theta));
    Vec3 e2(0, 1, 0);
    Vec3 e1 = e2 & e3;

    ret.push_back(combineInMatrix(e1, e2, e3));
  }

  for(int t = 0; t < n; ++t){
    double phi = t * (M_PI / 4) / n;
    Vec3 e3(cos(phi), sin(phi), 0);
    Vec3 e1(0, 0, -1);
    Vec3 e2 = e3 & e1;

    ret.push_back(combineInMatrix(e1, e2, e3));
  }

  for(int t = 0; t < n; ++t){
    double theta = t * (M_PI / 2) / (n - 1);
    Vec3 e3(cos(theta) / sqrt(2), cos(theta) / sqrt(2), sin(theta));
    Vec3 e2(-1/sqrt(2), 1/sqrt(2), 0);
    Vec3 e1 = e2 & e3;

    ret.push_back(combineInMatrix(e1, e2, e3));
  }

  return ret;
}


void make_map(const string& flnm, /*const string& surface_name,*/ const RotationMatrix& rm){
  int n = 601;
  double A = 0.002;

  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  t = rotateTensor(t, rm);

  //PovrayExport pov_surface(surface_name);

  QImage img_dest(n, n, QImage::Format_ARGB32_Premultiplied);
  for(int p = 0; p < n; ++p){
    cerr << p << " ";
    for(int q = 0; q < n; ++q){
      double y = A * ( 1 - 2. * p / (n - 1));
      double x = A * (-1 + 2. * q / (n - 1));

      PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, x, y);
      Poly pol = det(poly_mat);
      Poly::RootVec roots = pol.all_roots();

      int realNumber = 0;
      for(int t = 0; t < int(roots.size()); ++t)
	if (abs(imag(roots[t]))<1e-10){
	  ++realNumber;
	  //pov_surface.sphere(m * x, m * y, m * real(roots[t]));
	}

      realNumber /= 2;
      
      QColor col;
      if (x==0 || y == 0)
	col = qRgb(255, 255, 255);
      else 
	switch(realNumber){
	case(3): col = qRgb(255, 0, 0); break;
	case(2): col = qRgb(0, 255, 0); break;
	case(1): col = qRgb(0,   0, 255);break;
	default: col = qRgb(0, 0, 0);
	}

      img_dest.setPixel(q, p, col.rgb());

    }
  }

  cerr << endl << endl;
  img_dest.save(flnm.c_str());
}

void make_map_mult(){
  gsl_set_error_handler_off();
  vector<RotationMatrix> vrm = genMatSeq();
  
  for(int t = 0; t < int(vrm.size()); ++t){
    stringstream name_a_prefix, name_b_prefix, name_pov_prefix;

    name_a_prefix << "mult/a" << setfill('0') << setw(3) << t ;
    name_b_prefix << "mult/c" << setfill('0') << setw(3) << t ;
    name_pov_prefix << "mult/surf" << setfill('0') << setw(3) << t;

    string pov_name = name_a_prefix.str() + ".pov";
    string map_name = name_b_prefix.str() + ".png";
    string surface_name = name_pov_prefix.str() + ".pov";

    PovrayExport pov(pov_name);
    pov.matrix(vrm[t]);

    //make_map(map_name, surface_name, vrm[t]);
  }
}

void test_poly() {
  Poly pol(4, 0, 0, 0, 1);
  cout << "poly is " << pol << endl;
}

void test_mutant_christoffel_poly_matrix() {
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0.00, 0.00);
  cout << "poly matrix is " << endl;
  Poly pol = det(poly_mat);
  cout << "polynome is " << pol << endl;
  cout << "roots are: " << pol.all_roots() << endl;
}

void test_mutant_christoffel_rotate() {
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  boost::array<RotationMatrix::index, 2> shape = {{3,3}};
  RotationMatrix mat(shape);

  mat[0][0] =  0.5;
  mat[0][1] = -0.5;
  mat[0][2] = 1 / sqrt(2);

  mat[1][0] = -0.5;
  mat[1][1] =  0.5;
  mat[1][2] = 1 / sqrt(2);

  mat[2][0] = - 1 / sqrt(2);
  mat[2][1] = - 1 / sqrt(2);
  mat[2][2] = 0;

  t = rotateTensor(t, mat);

  PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0, 0);
  Poly pol = det(poly_mat);
  Poly::RootVec roots = pol.all_roots();
  cout << "vels: " << endl;
  for(int t = 0; t < (int) roots.size(); ++t)
    cout << 1. / roots[t] << " ";
}

void test_mutant_christoffel_shift() {
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0.0001, 0.0001);
  Poly pol = det(poly_mat);
  Poly::RootVec roots = pol.all_roots();
  cout << "vels: " << endl;
  for(int t = 0; t < (int) roots.size(); ++t)
    cout << 1. / roots[t] << " ";


  PovrayExport pr("pr.pov");
  for(double x = -1.5; x <= 1.5; x += 0.01) {
    cerr << "x = " << x << endl;
    for(double y = -1.5; y <= 1.5; y += 0.01) {
      PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0.001 * x, 0.001 * y);
      Poly pol = det(poly_mat);
      Poly::RootVec roots = pol.all_roots();
      vector<double> dat;
      for(int t = 0; t < int(roots.size()); ++t) {
	if (abs(roots[t].imag()) < 1e-5) {
	  dat.push_back(roots[t].real());
	}
      }
      sort(dat.begin(), dat.end());
      if(dat.size()) {
	double z = 1000 * dat.back();
	pr.sphere(x,y,z);
      }
    }
  }
}

void test_mutant_christoffel_shift_rotate() {
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  //boost::array<MaterialTensor::index, 2> shape = {{3,3}};
  RotationMatrix rm = makeAxisRotation(0.4, 0);
  t = rotateTensor(t, rm);
  
  PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0.0001, 0.0001);
  Poly pol = det(poly_mat);
  Poly::RootVec roots = pol.all_roots();
  cout << "vels: " << endl;
  for(int t = 0; t < (int) roots.size(); ++t)
    cout << 1. / roots[t] << " ";

  PovrayExport pr("pr.pov");
  for(double x = -1.5; x <= 1.5; x += 0.01){
    cerr << "x = " << x << endl;
    for(double y = -1.5; y <= 1.5; y += 0.01){
      PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0.001 * x, 0.001 * y);
      Poly pol = det(poly_mat);
      Poly::RootVec roots = pol.all_roots();
      vector<double> dat;
      for(int t = 0; t < int(roots.size()); ++t)
	if (abs(roots[t].imag()) < 1e-5)
	  dat.push_back(roots[t].real());
      sort(dat.begin(), dat.end());
      if(dat.size()) {
	double z = 1000 * dat.back();
	pr.sphere(x,y,z);
      }
    }
  }
}

void testZSliceNormal(const char* flnm) {
  ofstream dest(flnm);

  Vec3 axis1(1, 1, 0);
  Vec3 axis2(0, 0, 1);

  axis1.normalize();
  axis2.normalize();

  MaterialTensor  matTen = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);

  int n = 500;
  double dphi = 2 * M_PI / n;
  for(int t = 0; t < n; ++t){

    double phi = t * dphi;
    double x = cos(phi);
    double y = sin(phi);

    Vec3 testDir = axis1 * x + axis2 * y;

    Poly::RootVec gs=calcGammas(det(makeChristoffelPolyMatrix(matTen, testDir)));
    for(int t = 0; t < int(gs.size()); ++t){
      double gamma = real(gs[t]);
      double s = sqrt(5.96e3 / gamma);
      dest << x * s << " " << y * s << endl;
    }
  }
}

void testZSliceMutant(const char* flnm) {
  ofstream dest(flnm);
  Vec3 dir(1, 1, 0);
  dir.normalize();

  MaterialTensor  matTen = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);

  int n = 1000;
  double a = -0.002;
  double b =  0.002;
  double h = (b - a) / n;

  for(int t = 0; t < n; ++t){
    double x = a + h * t;
    Vec3 curr = dir * x;
    PolyMatrix poly_mat = makePolyMatrix(matTen, 5.96e3, curr.x(), curr.y());
    Poly pol = det(poly_mat);
    Poly::RootVec roots= pol.all_roots();
    for(int p = 0; p < int(roots.size()); ++p){
      if(imag(roots[p]) < 1e-14){
	dest << x << " " << real(roots[p]) << endl;
      }
    }
  }
}

void test_mutant_christoffel_chung_matrix() {
  PovrayExport pr("pr.pov");
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  RotationMatrix rm = genChungMatrix(10.0 / 180 * M_PI);

  t = rotateTensor(t, rm);
  
  for(double x = -2.0; x <= 2.0; x += 0.01) {
    cerr << "x = " << x << endl;
    for(double y = -2.0; y <= 2.0; y += 0.01) {
      PolyMatrix poly_mat = makePolyMatrix(t, 5.96e3, 0.001 * x, 0.001 * y);
      Poly pol = det(poly_mat);
      Poly::RootVec roots = pol.all_roots();

      for(int t = 0; t < int(roots.size()); ++t) {
	if (abs(roots[t].imag()) < 1e-5) {
	  double z = 1000 * roots[t].real();
	  pr.sphere(x,y,z);
	}
      }
    }
  }
}

void test_povray_export() {
  PovrayExport prx("sphere_on_x.pov");
  PovrayExport pry("sphere_on_y.pov");
  PovrayExport prz("sphere_on_z.pov");

  prx.sphere(1, 0, 0);
  pry.sphere(0, 1, 0);
  prz.sphere(0, 0, 1);
}

void test_composite_wave() {
  MaterialTensor mat_tens = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  MaterialTensor rot_mat_tens = rotateTensor(mat_tens, genChungMatrix(0. / 180 * M_PI));

  CompositeWave cv(0, 0, rot_mat_tens, 5.96, 1e8, Vec3(1, 0, 0));
}

void generate_teo2_xy_slice() {
  MaterialTensor  t = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  int N = 50000;
  double dphi = 2 * M_PI / N;

  ofstream dest("teo2_xy.txt");

  for(double phi = 0; phi < 2 * M_PI; phi += dphi) {
    Vec3 n(cos(phi), sin(phi), 0);
    RotationMatrix chr = makeChristoffelMatrix(t, n);
    Poly pol =makeCharacterPoly(chr);
    Poly::RootVec rv = calcGammas(pol);

    for(int t = 0; t < (int) rv.size(); ++t) {
      double s = sqrt(rho / real(rv[t]));
      dest << n.x() * s << " " << n.y() * s << endl;
    }
  }
}
