/// @file util.cpp
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#include "errors.h"
#include "matrix_fftw.h"
#include "plan_fftw.h"
#include "piezo_tensor.h"
#include "povray_export.h"
#include "spacial_matrix.h"
#include "volume_data.h"
#include "storage.h"
#include "util.h"
#include "vec3.h"
#include "vec3c.h"
#include "wave_matrix.h"

#include <iostream>
#include <fstream>
#include <map>
#include <mutex>
#include <stack>
#include <thread>

#ifdef BEAM_STRUCT_WORK
#include <QImage>
#endif

#define _USE_MATH_DEFINES
#include <math.h>

//#define DEBUG

namespace farn {

  using namespace std;

  MaterialTensor makeZeroTensor() {
    boost::array<MaterialTensor::index, 4> shape = {{3, 3, 3, 3}};
    MaterialTensor r(shape);

    int i, j, k, l;
    for(i = 0; i < 3; ++i)
      for(j = 0; j < 3; ++j)
	for(k = 0; k < 3; ++k)
	  for(l = 0; l < 3; ++l)
	    r[i][j][k][l] = 0;

    return r;
  }

  MaterialTensor makeTetragonalMaterialTensor(double c11,
                                              double c12,
                                              double c13,
                                              double c33,
                                              double c44,
                                              double c66) {
    MaterialTensor r = makeZeroTensor();

    r[0][0][0][0] = r[1][1][1][1] = c11;
    r[0][0][1][1] = r[1][1][0][0] = c12;
    r[0][0][2][2] = r[1][1][2][2] = r[2][2][0][0] = r[2][2][1][1] = c13;
    r[2][2][2][2] = c33;
    r[1][2][1][2] = r[1][2][2][1] = r[2][1][1][2] = r[2][1][2][1] = 
      r[0][2][0][2] = r[0][2][2][0] = r[2][0][0][2] = r[2][0][2][0] = c44;
    r[0][1][1][0] = r[0][1][0][1] = r[1][0][0][1] = r[1][0][1][0] = c66;

    return r;
  }

  MaterialTensor makeTrigonalMaterialTensor(double c11,
					    double c12,
					    double c13,
					    double c14,
					    double c33,
					    double c44,
					    double c66) {
    MaterialTensor r = makeZeroTensor();

    double arr[6][6] = {
      {c11,  c12, c13,  c14,   0,   0},
      {c12,  c11, c13, -c14,   0,   0}, 
      {c13,  c13, c33,    0,   0,   0},
      {c14, -c14,   0,  c44,   0,   0},
      {  0,    0,   0,    0, c44, c14},
      {  0,    0,   0,    0, c14, c66}};

    int dict[3][3] = {
      {1, 6, 5},
      {6, 2, 4},
      {5, 4, 3}};

    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
	for(int k = 0; k < 3; ++k)
	  for(int l = 0; l < 3; ++l){
	    int p = dict[i][j]; 
	    int q = dict[k][l];
		
	    r[i][j][k][l] = arr[p - 1][q - 1];
	  }

    return r;
  }

  //rho = 5960
  MaterialTensor makeParatelluriteMaterialTensor() {
    double c11 = 5.6e10;
    double c12 = 5.145e10;
    double c13 = 2.2e10;
    double c33 = 10.6e10;
    double c44 = 2.65e10;
    double c66 = 6.6e10;
    return makeTetragonalMaterialTensor(c11, c12, c13, c33, c44, c66);
  }

  //rho = 6210
  MaterialTensor makeTelluriumMaterialTensor() {
    //from 1979_Royer
    double c11=3.257e10, c12=0.845e10, c13=2.57e10, c14=1.238e10, c33=7.17e10, c44=3.094e10, c66=1.206e10;
    return makeTrigonalMaterialTensor(c11, c12, c13, c14, c33, c44, c66);
  }

  PiezoTensor makeTrigonal3mMaterialTensor(double e15, double e21, double e31, double e33){
    PiezoTensor ret;
    
    double mat[3][6]={
      {0,    0,     0,    0,   e15,   -e21},
      {-e21, e21,   0,  e15,     0,      0},
      {e31,  e31, e33,    0,     0,      0}
    };

    int dict[3][3] = {
      {1, 6, 5},
      {6, 2, 4},
      {5, 4, 3}};
    
    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q)
	for(int r = 0; r < 3; ++r)
	  ret(p,q,r) = mat[p][dict[q][r] - 1];

    return ret;
  }


  MaterialTensor rotateTensor(MaterialTensor const& t, RotationMatrix const& m) {
    boost::array<MaterialTensor::index, 4> shape = {{3, 3, 3, 3}};
    MaterialTensor ret(shape);

    int i, j, k, l;
    int p, q, r, s;

    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	  for (l = 0; l < 3; ++l)
	    {
	      double value = 0;
	      for (p = 0; p < 3; ++p)
		for (q = 0; q < 3; ++q)
		  for (r = 0; r < 3; ++r)
		    for (s = 0; s < 3; ++s)
		      value += t[p][q][r][s] * 
			m[i][p] * m[j][q] * m[k][r] * m[l][s];
	      ret[i][j][k][l] = value;
	    }
    return ret;
  }

  PolyMatrix makeZeroPolyMatrix()
  {
    boost::array<PolyMatrix::index, 2> shape = {{3, 3}};
    return PolyMatrix(shape);
  }

  PolyMatrix makePolyMatrix(MaterialTensor const& t, double gamma, double n1, double n2) {
    boost::array<PolyMatrix::index, 2> shape = {{3, 3}};
    PolyMatrix ret(shape);

    int i, j, k, l;

    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	vector<complex<double> > poly(3);
	for (k = 0; k < 3; ++k) {
	  size_t ind = (k==2)?1:0;
	  double v1 = (k==0)?n1:(k==1)?n2:1;
	  for (l = 0; l < 3; ++l) {
	    ind += (l==2)?1:0;
	    double v2 = (l==0)?n1:(l==1)?n2:1;
	    poly[ind] += t[i][k][j][l] * v1 * v2;
	  }
	}
	ret[i][j] = Poly(poly) - ((i==j)?Poly(gamma):Poly());
      }
    }

    return ret;
  }

  RotationMatrix transpose(const RotationMatrix& m){
    boost::array<PolyMatrix::index, 2> shape = {{3,3}};
    RotationMatrix ret(shape);
    
    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q)
	ret[p][q] = m[q][p];

    return ret;
  }

  RotationMatrix makeChristoffelMatrix(MaterialTensor const& t, const Vec3& v){
    boost::array<PolyMatrix::index, 2> shape = {{3,3}};
    RotationMatrix ret(shape);

    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q){
	ret[p][q] = 0;
	for(int r = 0; r < 3; ++r)
	  for(int s = 0; s < 3; ++s)
	    ret[p][q] += t[p][r][q][s] * v[r] * v[s];
      }
    
    return ret;
  }

  RotationMatrix makePiezoChristoffelMatrix(MaterialTensor const& t, const PiezoTensor& pt, const Vec3& v, double eps11, double eps33){
    boost::array<PolyMatrix::index, 2> shape = {{3,3}};
    RotationMatrix ret(shape);

    double eps = eps11 * (v[0] * v[0] + v[1] * v[1]) + eps33 * v[2] * v[2];
    //cout << "eps = " << eps << endl;

    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q){
	//cout << "p = " << p << endl;
	//cout << "q = " << q << endl;

	ret[p][q] = 0;

	double gamma_p = 0; 
	double gamma_q = 0; 
	for(int r = 0; r < 3; ++r)
	  for(int s = 0; s < 3; ++s){
	    ret[p][q] += t[p][r][q][s] * v[r] * v[s];
	    gamma_p += pt(r, p, s) * v[r] * v[s];
	    gamma_q += pt(r, q, s) * v[r] * v[s];
	  }

	//cout << "\tGamma = " << ret[p][q] << endl;
	//cout << "\tgamma_p = " << gamma_p << endl;
	//cout << "\tgamma_q = " << gamma_q << endl;
	double add = gamma_p * gamma_q / eps;
	//cout << "\tadd = " << add << endl;
	ret[p][q] += add;
      }
    
    return ret;
  }


  PolyMatrix makeChristoffelPolyMatrix(MaterialTensor const & r, Vec3 const& v){
    return makeChristoffelPolyMatrix(r, v[0], v[1], v[2]);
  }

  PolyMatrix makeChristoffelPolyMatrix(MaterialTensor const & t, double n1, double n2, double n3){
    boost::array<PolyMatrix::index, 2> shape = {{3,3}};
    PolyMatrix ret(shape);
    double r[]={n1,n2,n3};

    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	double val = 0;
	for(int p=0; p < 3; ++p)
	  for(int q=0; q < 3; ++q)
	    val += t[i][p][j][q] * r[p] * r[q];
	ret[i][j] = Poly(val) - ((i==j)?Poly(0.,1.):Poly());
      }
    }

    return ret;
  }

  ComplexMatrix makeZeroComplexMatrix()
  {
    boost::array<ComplexMatrix::index, 2> shape = {{3, 3}};
    return ComplexMatrix(shape);
  }

  ComplexMatrix evaluatePolyMatrix(PolyMatrix const& m, complex<double> x) {
    boost::array<ComplexMatrix::index, 2> shape = {{3, 3}};
    ComplexMatrix ret(shape);

    ComplexMatrix::index i, j;
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	ret[i][j] = m[i][j](x);
      }
    }

    return ret;
  }


  Poly det(PolyMatrix const& m) {
    return m[0][0] * m[1][1] * m[2][2]
      + m[0][1] * m[1][2] * m[2][0]
      + m[0][2] * m[1][0] * m[2][1]
      - m[0][2] * m[1][1] * m[2][0]
      - m[0][0] * m[1][2] * m[2][1]
      - m[0][1] * m[1][0] * m[2][2];
  }


  complex<double> det(ComplexMatrix const& m) {
    return m[0][0] * m[1][1] * m[2][2]
      + m[0][1] * m[1][2] * m[2][0]
      + m[0][2] * m[1][0] * m[2][1]
      - m[0][2] * m[1][1] * m[2][0]
      - m[0][0] * m[1][2] * m[2][1]
      - m[0][1] * m[1][0] * m[2][2];
  }


  RotationMatrix multiply(RotationMatrix const& a, RotationMatrix const& b) {
    if (a.shape()[1] != b.shape()[0]) {
      throw Error("Incorrect matrix sizes for multiplication");
    }

    int h = a.shape()[0];
    int w = b.shape()[0];

    boost::array<RotationMatrix::index, 2> shape = {{h,w}};
    RotationMatrix r(shape);

    size_t i, j, k;
    for (i = 0; i < a.shape()[0]; ++i) {
      for (j = 0; j < b.shape()[1]; ++j) {
	double value = 0;
	for (k = 0; k < a.shape()[1]; ++k) {
	  value += a[i][k] * b[k][j];
	}
	r[i][j] = value;
      }
    }

    return r;
  }


  RotationMatrix makeAxisRotation(double phi, int axis) {
    boost::array<RotationMatrix::index, 2> shape = {{3, 3}};
    RotationMatrix dat(shape);
        
    dat[(axis + 0) % 3][(axis + 0) % 3] = 1; 
    dat[(axis + 0) % 3][(axis + 1) % 3] = 0; 
    dat[(axis + 0) % 3][(axis + 2) % 3] = 0;

    dat[(axis + 1) % 3][(axis + 0) % 3] = 0; 
    dat[(axis + 1) % 3][(axis + 1) % 3] = cos(phi); 
    dat[(axis + 1) % 3][(axis + 2) % 3] = -sin(phi);

    dat[(axis + 2) % 3][(axis + 0) % 3] = 0; 
    dat[(axis + 2) % 3][(axis + 1) % 3] = sin(phi); 
    dat[(axis + 2) % 3][(axis + 2) % 3] = cos(phi);

    return dat;
  }

  RotationMatrix makeXYSlice(double phi) {
    return combineInMatrix(Vec3(-sin(phi), cos(phi), 0),
			   Vec3(        0,        0, 1),
			   Vec3( cos(phi), sin(phi), 0));
  }


  ostream& operator<<(ostream& os, const RotationMatrix& mat)
  {
    int i, j;
    for(i = 0; i < 3; ++i)
      {
	for(j = 0; j < 3; ++j)
	  os << mat[i][j] << " ";
	os << endl;
      }
    return os;
  }

  ostream& operator<<(ostream& os, const PolyMatrix& m)
  {
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
	os << i << "," << j << ") " << m[i][j] << endl;
    os << endl;
    return os;
  }

  Poly::RootVec calcGammas(Poly const&p)
  {
    Poly::RootVec r = p.all_roots(), ret;

    for(Poly::RootVec::const_iterator it = r.begin(); it != r.end(); ++it)
      if(it -> imag() <= 0)
	ret.push_back(*it);

    return ret;            
  }

  double vecAbs(ComplexVec const& v) {
    double c1 = abs(v[0]);
    double c2 = abs(v[1]);
    double c3 = abs(v[2]);

    return sqrt(c1 * c1 + c2 * c2 + c3 * c3);
  }

  ComplexVec norm(ComplexVec const& v)
  {
    
    double r = vecAbs(v);

    ComplexVec ret = list_of(v[0] / r) (v[1] / r) (v[2] / r);

    return ret;
  }

  ComplexVec kernel(ComplexMatrix const& m) {
    ComplexVec ret = list_of(m[0][1] * m[1][2] - m[1][1] * m[0][2])
      (m[0][2] * m[1][0] - m[0][0] * m[1][2])
      (m[0][0] * m[1][1] - m[1][0] * m[0][1]);
    return ret;
  }

  ComplexVec multiply(ComplexMatrix const & m, ComplexVec const& v) {
    ComplexVec ret = list_of(0)(0)(0);

    int i, j;
    for(i = 0; i < 3; ++i)
      for( j = 0; j < 3; ++j)
	ret[i] += m[i][j] * v[j];

    return ret;
  }

  ComplexMatrix matrixOfRows(Poly::RootVec const& r0,
			     Poly::RootVec const& r1,
			     Poly::RootVec const& r2) {
    ComplexMatrix ret = makeZeroComplexMatrix();
    Poly::RootVec const* rows[] = {&r0, &r1, &r2};

    int i, j;
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	ret[i][j] = (*rows[i])[j];
      }
    }
    return ret;
  }

  ComplexMatrix makeAmplitudeMatrix(MaterialTensor const& c,
				    ComplexMatrix const& q, /* Matrix of kernels */
				    ComplexMatrix const& n /* Matrix of directions */) {

    ComplexMatrix ret = makeZeroComplexMatrix();
    size_t i, j, k, l;
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	complex<double> value = 0;
	for (k = 0; k < 3; ++k) {
	  for (l = 0; l < 3; ++l) {
	    value += c[i][2][k][l] * q[j][k] * n[j][l];
	  }
	}
	ret[i][j] = value;
      }
    }
    return ret;
  }

  complex<double> calcDet(MaterialTensor const &t, 
			  RotationMatrix const &m,
			  double gamma, 
			  double n1, 
			  double n2)
  {
    MaterialTensor t1 = rotateTensor(t, m);
    PolyMatrix pm = makePolyMatrix(t1, gamma, n1, n2);
    Poly polyDet = det(pm);
    Poly::RootVec n3s = calcGammas(polyDet);

    if(n3s.size() != 3 )
      throw Error("Insufficient roots number");
        
    ComplexVec qs[3];
    ComplexMatrix n = makeZeroComplexMatrix();
    n[0][0] = n[1][0] = n[2][0] = n1; 
    n[0][1] = n[1][1] = n[2][1] = n2; 
        

    for(int ind = 0; ind < 3; ++ind)
      {
	complex<double> n3 = n3s[ind];
	n[ind][2] = n3;
	ComplexMatrix cm = evaluatePolyMatrix(pm, n3);
	qs[ind] = kernel(cm);
      }

    ComplexMatrix q = matrixOfRows(qs[0], qs[1], qs[2]);
    ComplexMatrix dmat = makeAmplitudeMatrix(t1, q, n);

    return det(dmat);
  }

  RotationMatrix combineInMatrix(const Vec3& e1, const Vec3& e2, const Vec3& e3){
    boost::array<RotationMatrix::index, 2> shape = {{3,3}};
    RotationMatrix ret(shape);

    for(int t = 0; t < 3; ++t){
      ret[0][t] = e1[t];
      ret[1][t] = e2[t];
      ret[2][t] = e3[t];
    }

    return ret;
  }

  Vec3 operator*(const RotationMatrix& m, const Vec3& v){
    Vec3 ret;

    for(int p = 0; p < 3; ++p){
      ret[p] = 0;
      for(int q = 0; q < 3; ++q)
	ret[p] += m[p][q] * v[q];
    }

    return ret;
  }

  //todo: investigate case of even n
  double ind2slow(int p, int n, double slow){
    if (2 * p < n)
      return p / slow;

    return (p - n) / slow;
  }

  ostream& getLog(){
    static ofstream ret("sz.txt");
    return ret;
  }

#ifdef BEAM_STRUCT_WORK
  WaveMatrix create_wave_matrix(int n, double aperture, 
				double freq, 
				const MaterialTensor& tens, 
				double rho, 
				const Vec3c& unitForce, 
				CalculationType ct){
    WaveMatrix ret(n);
    double slow(aperture * freq);
    double omega(2 * M_PI * freq);

    for(int p = 0; p < n; ++p){
#ifdef DEBUG
      cout << p << " of " << n << endl;
#endif

      for(int q = 0; q < n; ++q){
#ifdef DEBUG
	cerr << "(" << p << ", " << q << ")" << endl;
#endif
	double s1 = ind2slow(q, n, slow);
	double s2 = ind2slow(p, n, slow);

	ret(p, q) = CompositeWave(s1, s2, tens, rho, omega, unitForce, ct);
      }
      //getLog() << endl;
    }
    
    return ret;
  }

  WaveMatrix create_framed_wave_matrix(int n, 
				       double sxmin, double sxmax,
				       double symin, double symax,
				       double , 
				       double freq, 
				       const MaterialTensor& tens, 
				       double rho, 
				       const Vec3c& unitForce, 
				       CalculationType ct) {
    WaveMatrix ret(n);
    //double slow(aperture * freq);
    double omega(2 * M_PI * freq);

    for(int p = 0; p < n; ++p){
      cout << p << " of " << n << endl;
#ifdef DEBUG
      cout << p << " of " << n << endl;
#endif

      for(int q = 0; q < n; ++q){
#ifdef DEBUG
	cerr << "(" << p << ", " << q << ")" << endl;
#endif
	double s1 = sxmin + (sxmax - sxmin) / n * q;
	double s2 = symin + (symax - symin) / n * p;

	ret(p, q) = CompositeWave(s1, s2, tens, rho, omega, unitForce, ct);
      }
      //getLog() << endl;
    }
    
    return ret;
  }
#endif

  //wrong chang matrix
  // RotationMatrix genChungMatrix(double alpha){

  //   boost::array<RotationMatrix::index, 2> shape = {{3,3}};
  //   RotationMatrix mat(shape);
  
  //   double ca = cos(alpha);
  //   double sa = sin(alpha);
  //   double sq = sqrt(2);

  //   return combineInMatrix( Vec3( -1 / sq,   1 / sq,   0),
  // 			    Vec3( sa / sq,  sa / sq,  ca),
  // 			    Vec3( ca / sq,  ca / sq, -sa));  
  // }

  RotationMatrix genChungMatrix(double alpha){
    boost::array<RotationMatrix::index, 2> shape = {{3,3}};
    RotationMatrix mat(shape);
  
    double ca = cos(alpha);
    double sa = sin(alpha);
    double sq = sqrt(2.0);

    return combineInMatrix( Vec3(  -1 / sq,    1 / sq,   0),
			    Vec3(  sa / sq,   sa / sq,  ca),
			    Vec3(  ca / sq,   ca / sq,  -sa));  
  }

  RotationMatrix nearXMatrix(double alpha) {
    boost::array<RotationMatrix::index, 2> shape = {{3,3}};
    RotationMatrix mat(shape);
  
    double ca = cos(alpha);
    double sa = sin(alpha);

    return combineInMatrix( Vec3( -sa, ca, 0), 
			    Vec3(   0,  0, 1),
			    Vec3(  ca, sa, 0));
  }

  Vec3c calcPol(const ComplexMatrix& dat) {
#ifdef DEBUG
    cerr << "calcPol" << endl;
    cerr << "matrix:" << endl << dat << endl;
#endif

    Vec3c a1(dat[0][0], dat[0][1], dat[0][2]);
    Vec3c a2(dat[1][0], dat[1][1], dat[1][2]);
    Vec3c a3(dat[2][0], dat[2][1], dat[2][2]);

    Vec3c var[]={a1, a2, a3};

    for(int p = 0; p < 3; p++) {
      if (var[p].getSuperModule() > 1e-5 && 
	  var[(p + 1) % 3].getSuperModule() < 1e-5 &&
	  var[(p + 2) % 3].getSuperModule() < 1e-5) {
	Vec3c nz = var[p];
	Vec3c c1(-nz.y(), nz.x(), 0);
	Vec3c c2(0, -nz.z(), nz.y());
	if (c1.getSuperModule() > 1e-5) {
	  c1.normalize();
	  return c1;
	} else {
	  c2.normalize();
	  return c2;
	}
      } 
    }

    //null vector case
    double l1 = a1.abs();
    double l2 = a2.abs();
    double l3 = a3.abs();

    double desc_eps = 1e-7;
    if( min(min(l1,l2), l3) < max(max(l1,l2),l3) * desc_eps){
#ifdef DEBUG
      cerr << "\t a1 = " << a1 << endl;
      cerr << "\t a2 = " << a2 << endl;
      cerr << "\t a3 = " << a3 << endl << endl;
      cerr << "\t l1 = " << l1 << endl;
      cerr << "\t l2 = " << l2 << endl;
      cerr << "\t l3 = " << l3 << endl << endl;
#endif

      if(l1 <= min(l2, l3)) {
#ifdef DEBUG
	cerr << "l1 case: a2, a3 = " << (a2 & a3) << endl;
#endif
	return (a2 & a3).norm();
      }

      if(l2 <= min(l1, l3)) {
#ifdef DEBUG
	cerr << "l2 case: a1, a3 = " << (a1 & a3) << endl;
#endif
	return (a1 & a3).norm();
      }
#ifdef DEBUG
      cerr << "l3 case: a1, a2 = " << (a1 & a2) << endl;
#endif

      return (a1 & a2).norm();
    }

    //not null vector case
    a1.normalize();
    a2.normalize();
    a3.normalize();

#ifdef DEBUG
    cerr << "\t not null case" << endl;
    cerr << "\t a1 = " << a1 << endl;
    cerr << "\t a2 = " << a2 << endl;
    cerr << "\t a3 = " << a3 << endl << endl;
#endif

    Vec3c r1 = (a1 & a2);
    Vec3c r2 = (a2 & a3);
    Vec3c r3 = (a3 & a1);

#ifdef DEBUG
    cerr << "\t r1 = " << r1 << endl;
    cerr << "\t r2 = " << r2 << endl;
    cerr << "\t r3 = " << r3 << endl << endl;
#endif

    double v1 = r1.abs();
    double v2 = r2.abs();
    double v3 = r3.abs();

#ifdef DEBUG
    cerr << "\tv1 = " << v1 << endl;
    cerr << "\tv2 = " << v2 << endl;
    cerr << "\tv3 = " << v3 << endl << endl;
#endif

    if (v1 >= max(v2, v3)) {
#ifdef DEBUG
      cerr << "\t a1 a2 selected" << endl;
      cerr << "\t r1 = " << r1 << endl;
#endif
      return r1.norm();
    }

    if (v2 >= max(v1, v3)) {
#ifdef DEBUG
      cerr << "\t a2 a3 selected" << endl;
      cerr << "\t r2 = " << r2 << endl;
#endif
      return r2.norm();
    }

#ifdef DEBUG
    cerr << "\t a3 a1 selected" << endl;
    cerr << "\t r3 = " << r3 << endl;
#endif
    return r3.norm();
  }

  Mat3 strainFromKQ(const Vec3c& waveVector, const Vec3c& polarization){
    Mat3 ret;
    
    for(int p = 0; p < 3; ++p){
      for(int q = 0; q < 3; ++q){
	ret[p][q] = 1./2 * (waveVector[p] * polarization[q] + polarization[p] * waveVector[q]);
      }
    }
    return ret;
  }

  Mat3 stressFromCS(const MaterialTensor& tens, const Mat3& S){
    Mat3 ret;

    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q){
	CD val = 0;
	for(int r = 0; r < 3; ++r)
	  for(int s = 0; s < 3; ++s){
	    double ct = tens[p][q][r][s];
	    complex<double> cs = S[r][s];
	    val += cs * ct;
	  }
	ret[p][q]=val;
      }

    return ret;
  }

  Mat3 columns2matrix(const Vec3c& v1, const Vec3c& v2, const Vec3c& v3){
    Mat3 ret;

    ret[0][0] = v1[0]; ret[0][1] = v2[0]; ret[0][2] = v3[0];
    ret[1][0] = v1[1]; ret[1][1] = v2[1]; ret[1][2] = v3[1];
    ret[2][0] = v1[2]; ret[2][1] = v2[2]; ret[2][2] = v3[2];


    return ret;
  }

  Poly makeCharacterPoly(const RotationMatrix& a){
    vector<complex<double> > ret(4);

    ret[3] = 1;
    ret[2] = -a[0][0] - a[1][1] - a[2][2];
    ret[1] = a[1][1] * a[2][2] - a[1][2] * a[2][1] + a[0][0] * a[2][2] - a[0][2] * a[2][0] + a[0][0]*a[1][1] - a[0][1] * a[1][0];
    ret[0] = -(a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
    
    return Poly(ret);
  }

#ifdef FFTW_WORK
  MatrixFFTW loadFromPicture(string flnm, double low, double hi){
    QImage sour(flnm.c_str());

    int h = sour.height();
    int w = sour.width();

    if (h == 0 || w == 0) {
      throw(string("loadFromPicture: can't load file " + flnm));
    }

    MatrixFFTW ret(h, w);
    for(int p = 0; p < h; ++p)
      for(int q = 0; q < w; ++q)
	ret(p,q) = qGray(sour.pixel(q,p))/255.0 * (hi - low) + low;
    
    return ret;
  }

  void saveAsPicture(const MatrixFFTW& a, string flnm) {
    cerr << "saveAsPicture: flnm = " << flnm << " h = " << a.height() << " w = " << a.width() << endl;

    QImage dest(a.width(), a.height(), QImage::Format_ARGB32_Premultiplied);

    int w = a.width();
    int h = a.height();
    double minv = 0, maxv = 1;
    minv = maxv = abs(a(0,0));
    for(int p = 0; p < h; ++p) {
      for(int q = 0; q < w; ++q) {
	double v = abs(a(p,q));
	minv = min(minv, v);
	maxv = max(maxv, v);
      }
    }

    for(int p = 0; p < h; ++p)
      for(int q = 0; q < w; ++q){
	double v = abs(a(p,q));
	int val = 0;
	//if ((v-minv) < 0.001 * (maxv - minv)) {
	val = int(255 * (v - minv) / (maxv - minv));
	/*} else {
	  val = 255;
	  }*/
	dest.setPixel(q, p, qRgb(val, val, val));
      }

    if (dest.save(flnm.c_str())) {
      cerr << "success" << endl;
    } else {
      cerr << "it is impossible to save '" << flnm << "'" << endl;
    }
  } 


#endif

#ifdef BEAM_STRUCT_WORK
  void saveAsPictures(const Storage& dat, string flnmBase){
    MatrixFFTW mf(dat.depth(), dat.width());
 
    cerr << "storage h d w = " << dat.height() << " " << dat.depth() << " " << dat.width() << endl;
    cerr << "mf h w = " << mf.height() << " " << mf.width() << endl;
   
    for(int t = 0; t < dat.height(); ++t){
      stringstream flnm;
      flnm << flnmBase << t << ".png";

      copySliceIntoFftw(t, dat, mf);
      saveAsPicture(mf, flnm.str());
    }
  }

  mutex layerTransformMutex;

  void threadTransform(stack<int>* ptasks, Storage* ret, const Storage* dat, MatrixFFTW* sour, MatrixFFTW* dest, const PlanFFTW* plan) {
    while(ptasks->size()) {
      int task;
      layerTransformMutex.lock();
      task = ptasks->top();
      ptasks->pop();
      layerTransformMutex.unlock();

      copySliceIntoFftw(task, *dat, *sour);
      plan -> execute();
      copyFftwIntoSlice(task, *dest, *ret);
    }
  }

  Storage layerTransformThreaded(const Storage& dat, MatrixFFTW& sour, MatrixFFTW& dest, const PlanFFTW& plan){
    int h = dat.height();
    int d = dat.depth();
    int w = dat.width();

    Storage ret(h, d, w);
    stack<int> tasks;

    for(int t = 0; t < h; ++t){
      tasks.push(t);
    }

    int ns = 4;
    thread threads[ns];
    for(int p = 0; p < ns; ++p) {
      threads[p] = thread(threadTransform, &tasks, &ret, &dat, &sour, &dest, &plan);
    }
    
    for(int p = 0; p < ns; ++p) {
      threads[p].join();
    }

    return ret;
  }

  Storage layerTransform(const Storage& dat, MatrixFFTW& sour, MatrixFFTW& dest, const PlanFFTW& plan){
    int h = dat.height();
    int d = dat.depth();
    int w = dat.width();

    Storage ret(h, d, w);

    for(int t = 0; t < h; ++t) {
      copySliceIntoFftw(t, dat, sour);
      plan.execute();
      copyFftwIntoSlice(t, dest, ret);
    }

    return ret;
  }

  void copySliceIntoFftw(int t, const Storage& sour, MatrixFFTW& dest){
    if (dest.width() != sour.width()){
      throw(string("copySliceIntoFftw: width of storage and fftwMatrix mismatched"));
    }

    if (dest.height() != sour.depth()){
      throw(string("copySliceIntoFftw: height of fftwMatrix not equal to storage depth"));
    }

    if (t >= sour.height()){
      throw(string("copySliceIntoFftw: height of storage exceeded"));
    }

    int h = dest.height();
    int w = dest.width();

    for(int p = 0; p < h; ++p){
      for(int q = 0; q < w; ++q){
	dest(p, q) = sour(t, p, q);
      }
    }
  }

  void copyFftwIntoSlice(int t, const MatrixFFTW& sour, Storage& dest){
    if (sour.width() != dest.width()){
      throw(string("copyFftwIntoSlice: width of storage and fftwMatrix mismatched"));
    }

    if (sour.height() != dest.depth()){
      throw(string("copyFftwIntoSlice: height of fftwMatrix not equal to storage depth"));
    }

    if (t >= dest.height()){
      throw(string("copyFftwIntoSlice: height of storage exceeded"));
    }

    int h = sour.height();
    int w = sour.width();

    double mult = 1./(h * w);

    for(int p = 0; p < h; ++p){
      for(int q = 0; q < w; ++q){
	dest(t, p, q) = sour(p, q) * mult;
      }
    }
  }

 
  SpacialMatrix getSpacialMatrixFromStorage(const Storage& stor){
    int h = stor.depth();
    int w = stor.width();

    SpacialMatrix ret(h, w);
    
    for(int p = 0; p < h; ++p)
      for(int q = 0; q < w; ++q) 
	ret(p, q).takeFromStorage(stor, p, q);
    
    return ret;
  }
#endif

  bool rightPointing(double s1, 
		     double s2, 
		     double s3, 
		     const PolyMatrix& poly_mat, 
		     const MaterialTensor& tens, 
		     double omega) {
    PlaneWave pw(s1, s2, s3, poly_mat, tens, omega);

    return real(pw.getPointing().z()) < 0; 
  }

  Vec3c operator*(const ComplexMatrix&l, const Vec3c& r) {
    Vec3c ret;

    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q)
	ret[p] += l[p][q] * r[q];

    return ret;
  }

  ostream& operator << (ostream& os, const ComplexMatrix& cm) {
    for(int p = 0; p < 3; ++p) {
      for(int q = 0; q < 3; ++q) {
	os << cm[p][q] << " ";
      }
      os << endl;
    }
    return os;
  }

  RotationMatrix minusGamma(const RotationMatrix& r, double gamma) {
    RotationMatrix ret(r);
    
    for(int p = 0; p < 3; ++p) {
      ret[p][p] -= gamma;
    }

    return ret;
  }

  Vec3 takeAbs(const Vec3c& r) {
    return Vec3(abs(r.x()), abs(r.y()), abs(r.z()));
  }

  Vec3 takeReal(const Vec3c& r) {
    return Vec3(real(r.x()), real(r.y()), real(r.z()));
  }

  VolumeData
  newSphericalVolumeData(int n) {
    double ap = 18;

    VolumeData ret(Storage(n, n, n), 
		   Mat3(Vec3c(1, 0, 0), Vec3c(0, 1, 0), Vec3c(0, 0, 1)), 
		   Vec3(ap, ap, ap),
		   Vec3(-ap/2, -ap/2, -ap/2));
  
    for(int r=0; r<n; r++) {
      for(int p=0; p<n; p++) {
	for(int q=0; q<n; q++) {
	  double x = -ap/2 + ap/(n-1)*q;
	  double y = -ap/2 + ap/(n-1)*p;
	  double z = -ap/2 + ap/(n-1)*r;

	  ret.vd(r, p, q) = sin(sqrt(x*x + y*y + z*z))*2 + 5;
	}
      }
    }

    return ret;
  }

  VolumeData
  newBarZVolumeData(int n, const Mat3& cs) {
    double ap = 18;

    VolumeData ret(Storage(n, n, n), 
		   cs, 
		   Vec3(ap, ap, ap),
		   Vec3(-ap/2, -ap/2, -ap/2));
  
    double x0 = -5;
    double y0 = 5; 
    //double z0 = 5;

    for(int r=0; r<n; r++) {
      for(int p=0; p<n; p++) {
	for(int q=0; q<n; q++) {
	  double x = -ap/2 + ap/(n-1)*q;
	  double y = -ap/2 + ap/(n-1)*p;
	  //double z = -ap/2 + ap/(n-1)*r;


	  double dx = (x-x0);
	  double dy = (y-y0);
	  //double dz = (z-z0);

	  double val = 0;
	  if (dx*dx + dy*dy < 3*3) {
	    val = 1;
	  } else {
	    val = 0;
	  }

	  ret.vd(r, p, q) = val;
	}
      }
    }

    return ret;
  }

  void drawVolumePovray(const std::string& flnm, const VolumeData& vd, double thr) {
    PovrayExport pe(flnm, NO_HEADER);

    cerr << "cs : " << endl << vd.cs << endl;

    size_t w=vd.vd.width();
    size_t d=vd.vd.depth();
    size_t h=vd.vd.height();

    double rad = min(vd.geometry.x()/w, min(vd.geometry.y()/d, vd.geometry.z()/h))/2 * 1.5;
  
    double rthr = vd.vd.calcRealThreshold(thr);

    for(int r=0; r<int(h); r++) {
      for(int p=0; p<int(d); p++) {
	for(int q=0; q<int(w); q++) {
	  double x=q*vd.geometry.x()/w;
	  double y=p*vd.geometry.y()/d;
	  double z=r*vd.geometry.z()/h;

	  bool skip = true;
	  for(int dr = -1; dr <= 1 && skip; dr++) {
	    for(int dp = -1; dp <= 1 && skip; dp++) {
	      for(int dq = -1; dq <= 1 && skip; dq++) {
		if (dr!=0 && dp!=0 && dq!=0) {
		  if (r+dr < 0 || r+dr >= vd.vd.height()) {
		    skip = false;
		  } else if (p+dp<0 || p+dp >= vd.vd.depth()) {
		    skip = false;
		  } else if (q+dq<0 || q+dq >= vd.vd.width()) {
		    skip = false;
		  } else if (abs(vd.vd(r+dr, p+dp, q+dq))<rthr) {
		    skip = false;
		  }
		}
	      }
	    }
	  }

	  Vec3 currPoint = vd.origin + takeReal(vd.cs*Vec3c(x,y,z));
	  if (abs(vd.vd(r, p, q))>rthr && !skip) {
	    pe.sphere(currPoint.x(), currPoint.y(), currPoint.z(), rad);
	  }
	}
      }
    }
  }


double brezenDeviat(const Vec3& a, const Vec3& dir) {
  double r = a.z();
  double p = a.y();
  double q = a.x();

  double dx = p*dir.z() - q*dir.y();
  double dy = dir.x()*q - r*dir.z();
  double dz = r*dir.y() - dir.x() * p;

  return dx*dx + dy*dy + dz*dz;
}

}
