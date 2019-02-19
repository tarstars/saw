#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

#include "utils.h"
#include "vec3.h"
#include "coeff.h"
#include "mat3.h"
#include "poly3.h"
#include "poly.h"
#include "veccomp.h"


Mat3 n2crist(const Coeff& m, const Vec3 &n)
{
  Mat3 ret;

  double n1=n.dat[0];
  double n2=n.dat[1];
  double n3=n.dat[2];


  double g11 = m.c11 * n1 * n1 + m.c66 * n2 * n2 + m.c44 * n3 * n3 + 2 * m.c16 * n1 * n2;
  double g12 = m.c16 * ( n1 * n1 - n2 * n2 ) + (m.c12 + m.c66) * n1 * n2;
  double g13 = (m.c13 + m.c44) * n1 * n3;
  double g22 = m.c66 * n1 * n1 + m.c11 * n2 * n2 + m.c44 * n3 * n3 - 2 * m.c16 * n1 * n2;
  double g23 = (m.c13 + m.c44) * n2 * n3;
  double g33 = m.c44 * (n1 * n1 + n2 * n2) + m.c33 * n3 * n3;
  double g21 = g12;
  double g31 = g13;
  double g32 = g23;

  ret[0][0] = g11; ret[0][1] = g12; ret[0][2] = g13;
  ret[1][0] = g21; ret[1][1] = g22; ret[1][2] = g23;
  ret[2][0] = g31; ret[2][1] = g32; ret[2][2] = g33;

  return ret;
}

CPoly3 mat2poly(const Mat3& m)
{
  return CPoly3( 1,

		 -m[0][0] - m[1][1] - m[2][2],

		 m[1][1] * m[2][2] + m[0][0] * m[2][2] + m[0][0] * m[1][1] - 
		 m[1][2] * m[1][2] - m[0][2] * m[0][2] - m[0][1] * m[0][1],

		 m[1][1] * m[0][2] * m[2][0] + m[2][2] * m[0][1] * m[1][0] + m[0][0] * m[1][2] * m[2][1] -
		 m[0][0] * m[1][1] * m[2][2] - m[0][1] * m[1][2] * m[2][0] - m[0][2] * m[1][0] * m[2][1]);

}

vector<CEigElement> vel_polar(const Vec3 &n)
{
  Mat3 crist(n2crist(Coeff::paratellurite(),n)); 

  CPoly3 poly(mat2poly(crist));

  CD x1 = poly.find_root();
  CD x2;
  CD x3;
    
  poly.find_two_root(&x2, &x3, x1);

  //roots, Mat3 -> vector<EigElements>
  Vec3 p1(crist.row(0));
  Vec3 p2(crist.row(1));

  Vec3 d1=polariz(p1, p2, x1);
  Vec3 d2=polariz(p1, p2, x2);
  Vec3 d3=polariz(p1, p2, x3);

  vector<CEigElement> ret;
  ret.push_back(CEigElement(d1, n, x1));
  ret.push_back(CEigElement(d2, n, x2));
  ret.push_back(CEigElement(d3, n, x3));
  sort(ret.begin(), ret.end());

  return ret;
}

void one_stage(double go,ostream& dest3d)
{
    stringstream sflnmdest;
    sflnmdest<<int(go * 100);
    string flnmdest=sflnmdest.str();
    flnmdest="txt/"+string(4-flnmdest.size(),'0')+flnmdest+".txt";

    cerr<<flnmdest<<endl;
    ofstream dest(flnmdest.c_str());

    Vec3 e1(1, 0, 0);
    Vec3 e2(0, 1, 0);
    Vec3 e3(0, 0, 1);

    //    double phi = 0;
    Vec3 plane(sin(go),0,cos(go));

    Vec3 p1(e1),p2;

    if( (plane & p1).abs() < eps)
        p1=e2;
    
    p2=(p1 & plane).norm();
    p1=(plane & p2).norm();

    for(double alpha = 0; alpha < 6.28; alpha += 0.001)
    {
        Vec3 n=p1 * cos(alpha) + p2 * sin(alpha);

        vector<CEigElement> dat;
        dat=vel_polar(n);

        dest << setprecision(15) << alpha << " " 
	     << dat[0].slow << " " << dat[1].slow 
	     << " " << dat[2].slow << endl;

	dest3d << (dat[0].slow * n) << endl;
    }
}

void generate3d()
{
    ofstream dest3d("3d.txt");
    for(double go=0; go < M_PI/2; go += 0.01)
        one_stage(go,dest3d);
}

Poly n1n2g_to_poly(long double n1, long double n2, long double g)
{
  long double c11 = 5.6;
  long double c12 = 5.1;
  long double c13 = 2.2;
  long double c33 = 10.6;
  long double c44 = 2.65;
  long double c66 = 6.6;
  //long double c16 = 0;

  Poly g11( c11 * n1 * n1 + c66 * n2 * n2  - g, 0. , c44);
  Poly g12( (c12 + c66) * n1 * n2);
  Poly g13( 0, (c13 + c44) * n1 );
  Poly g22( c66 * n1 * n1 + c11 * n2 * n2 - g, 0 , c44);
  Poly g23( 0, (c13 + c44) * n2);
  Poly g33( c44 * (n1 * n1 + n2 * n2) - g, 0,  c33);
  Poly g21( g12);
  Poly g31( g13);
  Poly g32( g23);

  return g11 * g22 * g33 + g12 * g23 * g31 + g13 * g21 * g32 -
    g31 * g13 * g22 - g32 * g23 * g11 - g33 * g21 * g12;
}

poly_type::value_type
det_val(long double  n1, long double n2, long double g)
{
  cerr << "determinant calculation" << endl;
  cerr << "\tn1 n2 g = " << n1 << " " << n2 << " " << g << endl;

  long double c11 = 5.6;
  long double c12 = 5.1;
  long double c13 = 2.2;
  long double c33 = 10.6;
  long double c44 = 2.65;
  long double c66 = 6.6;
  //long double c16 = 0;

  Poly g11( c11 * n1 * n1 + c66 * n2 * n2  - g, 0. , c44);
  Poly g12( (c12 + c66) * n1 * n2);
  Poly g13( 0, (c13 + c44) * n1 );
  Poly g22( c66 * n1 * n1 + c11 * n2 * n2 - g, 0 , c44);
  Poly g23( 0, (c13 + c44) * n2);
  Poly g33( c44 * (n1 * n1 + n2 * n2) - g, 0,  c33);
  Poly g21( g12);
  Poly g31( g13);
  Poly g32( g23);

  Poly pol = g11 * g22 * g33 + g12 * g23 * g31 + g13 * g21 * g32 -
    g31 * g13 * g22 - g32 * g23 * g11 - g33 * g21 * g12;

  cerr << "\tpoly = " << pol << endl;

  poly_type ar = all_roots(pol);
  poly_type arf;

  for(int t = 0; t < (int)ar.size(); ++t)
    if(ar[t].imag() > 0)
      arf.push_back(ar[t]);

  vector<VecComp3> pols;
  for(int t = 0; t < (int)arf.size(); ++t)
    {
      poly_type::value_type root = arf[t];

      VecComp3 a1(g11(root), g12(root), g13(root));
      VecComp3 a2(g21(root), g22(root), g23(root));
      VecComp3 a3(g31(root), g32(root), g33(root));

      VecComp3 p1(a1 & a2);
      VecComp3 p2(a2 & a3);
      VecComp3 p3(a3 & a1);

      p1 = p1.norm();
      p2 = p2.norm();
      p3 = p3.norm(); 

      pols.push_back(p1);
    }

  poly_type::value_type n31 = arf[0];
  poly_type::value_type n32 = arf[1];
  poly_type::value_type n33 = arf[2];

  cerr << "\t polarizations " << endl;
  for(int r = 0; r < 3; ++r)
    cerr << "\t\t" << pols[r] << endl;

  //arf = {n31, n32, n33}
  //pols = {u1, u2, u3}

  if(arf.size() < 3) return 0;

  vector<poly_type::value_type> z(9);
  for(int r = 0; r < 3; ++r)
    {
      z[3 * 0 + r] = arf[r] * pols[r](1) + n1 * pols[r](3);
      z[3 * 1 + r] = arf[r] * pols[r](2) + n2 * pols[r](3);
      z[3 * 2 + r] = (arf[r] * pols[r](1) + n1 * pols[r](3) + arf[r] * pols[r](2) + n2 * pols[r](3)) * c13 + 
	(arf[r] * pols[r](3) ) * c33;
    }

  cerr << "\tsecond matrix " << endl;
  for(int p = 0; p < 3; ++p)
    {
      cerr << "\t\t";
      for(int q = 0; q < 3; ++q)
	cerr << z[p * 3 + q] << " ";
      cerr << endl;
    }
  cerr << endl;

  int di[][4]={{0,4,8,1},{1,5,6,1},{2,3,7,1},{2,4,6,-1},{0,5,7,-1},{1,3,8,-1}};

  poly_type::value_type det = 0;

  for(int t = 0; t < 6; ++t)
    det += z[di[t][0]] * z[di[t][1]] * z[di[t][2]] * (long double)(di[t][3]);

  return det;
}
