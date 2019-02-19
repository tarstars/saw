#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <complex>

#include "poly.h"
#include "coeff.h"
#include "tensor.h"
#include "vec3.h"
#include "veccomp.h"
#include "eigelement.h"
#include "mat3.h"
#include "poly3.h"
#include "utils.h"

using namespace std;


vector<double> 
roots_of_3(const CPoly3& poly)
{
  CD x1 = poly.find_root();

  vector<double> ret;

  if(imag(x1) < 1e-15)
    {
      CD x2;
      CD x3; 

      try{
	poly.find_two_root(&x2, &x3, x1);

	ret.push_back(real(x1));
	ret.push_back(real(x2));
	ret.push_back(real(x3));
      }
      catch(string msg)
	{
	  cerr << msg << endl;
	}
    }

  return ret;
}



vector<double> dist(double x, double y, double z)
{
  Tensor tenb(Coeff::paratellurite());
  
  Tensor ten(tenb.rotate(Mat3(M_PI / 4,0)));

  /*
  Mat3 rt(M_PI / 5, 0);
  for(int t = 0; t < 10; ++t)
    ten = ten.rotate(rt);
  */

  Vec3 vec(x, y, z);
  Mat3 cr = ten.crist(vec);

  vector<double> ret = roots_of_3(mat2poly(cr));
  
  for(int t = 0; t < (int)ret.size(); ++t)
    {
      double v = sqrt(ret[t] / real(rho));
      ret[t] = 1 / v;
    }
  
  return ret;
}

void 
d3()
{
  ofstream dest("sphere.txt");

  int n = 3000;
  double dz = 2.0 / n;
  //double dphi = 2 * M_PI / n;
  
  int meter = 0;
  stringstream med;
  for(double z = -1; z <= 1; z+=dz)
    {
      //for(double phi = 0; phi < 2 * M_PI; phi += dphi)
      double phi = 2 * M_PI * n * n * asin(z);
	{
	  double theta = asin(z);

	  double x = cos(theta) * cos(phi);
	  double y = cos(theta) * sin(phi);
	  double z = sin(theta);

	  vector<double> rv = dist(x, y, z);
	  
	  for(int t = 0; t < (int)rv.size(); ++t)
	    {
	      ++meter;
	      double r = rv[t];
	      med << x * r << " " << y * r << " " << z * r << endl;
	    }
	}
    }

  dest << "object 1 class array type float rank 1 shape 3 items " << meter << " data follows" << endl;
  dest << med.str();
  dest << "object \"my scalars\" class field" << endl;
  dest << "component \"positions\" value 1" << endl;
}

void testset()
{
   Tensor ten(Coeff::paratellurite());

   Mat3 rotst(0.3, 0), rotinv(-0.3, 0);

   Vec3 n(Vec3::theta_phi2vec(0.1,0.1));

   Mat3 one(ten.rotate(rotst).crist(n));
   Mat3 two(ten.crist(rotst * n));

   cout << one << endl;
   cout << endl;
   cout << two << endl;

   vector<double> oned = roots_of_3(mat2poly(one));
   vector<double> twod = roots_of_3(mat2poly(two));

   for(int t = 0; t < (int)oned.size(); ++t)
     cout << oned[t] << " ";
   cout << endl;

   for(int t = 0; t < (int)oned.size(); ++t)
     cout << twod[t] << " ";
   cout << endl;

   cout << "from full tensor" << endl;
   Poly pn3 = ten.poly_by_n3(cos(0.6), sin(0.6), 1);
   cout << pn3 << endl;

   cout << "from delesan formulation" << endl;
   Poly pdel = n1n2g_to_poly(cos(0.6), sin(0.6), 1);
   cout << pdel << endl;

   cout << "roots are " << endl;
   poly_type roots = all_roots(pn3);
   for(int t = 0; t < (int)roots.size(); ++t)
     cout << roots[t] << " ";
   cout << endl;

   Mat3 td;
   td[0][0] = 1; td[0][1] = 7; td[0][2] = 25;
   td[1][0] = 43; td[1][1] = 217; td[1][2] = 426;
   td[2][0] = 33; td[2][1] = 15; td[2][2] = 705;
   
   cout << "matrix " << endl << td << endl;
   cout << endl;
   cout << "det = " << td.det() << endl;

   Tensor ten_source(Coeff::paratellurite());
   ofstream dest("cpp_tensor.txt");
   (ten_source.rotate(Mat3(0.3,0) * Mat3(0.3, 1))).mupad_style_output(dest);

}


void test_matrix()
{
  cout << "x" << endl << Mat3(0.1, 0) << endl;
  cout << endl;
  cout << "y" << endl << Mat3(0.1, 1) << endl;
  cout << endl;
  cout << "z" << endl << Mat3(0.1, 2) << endl;
  cout << endl;
  cout << "x * y" << endl << Mat3(0.3, 0) * Mat3(0.3, 1) << endl;
  cout << endl;
}


int main()
{

   Tensor ten_source(Coeff::paratellurite());

   Mat3 coord_transform(Mat3(0.3, 0) * Mat3(0.3, 1));

   Tensor ten(ten_source.rotate(coord_transform));

   LD phi = 0.4;
   Vec3 vec_n(cos(phi), sin(phi), 0);

   CPoly3 gamma_poly = mat2poly(ten.crist(vec_n));
   VD gammas = roots_of_3(gamma_poly);

   cout << "gamma polynome " << gamma_poly << endl; 

   cout << "gamma values ";
   for(int t = 0; t < (int)gammas.size(); ++t)
     cout << setprecision(10) << gammas[t] << " ";
   cout << endl;

   double gamma_guess = gammas[0];

   Poly crist_poly = ten.poly_by_n3(vec_n[0], vec_n[1], gamma_guess);

   cout << "poly by n3 " << crist_poly << endl;

   poly_type roots = all_roots(crist_poly);
   
   cout << "n3 values ";
   for(int t = 0; t < (int)roots.size(); ++t)
     cout << roots[t] << " ";
   cout << endl;
   

   return 0;
}

