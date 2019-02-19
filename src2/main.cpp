/// @file main.cpp
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#include "errors.h"
#include "util.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace farn;
using namespace std;

// m = multiply(makeAxisRotation(0.3, 0), makeAxisRotation(0.4, 1))
// \gamma_1 = 2.182 - 2.186 
// \gamma_2 = 2.646 - 2.650
// \gamma_3 = 5.47 - 5.48
// \gamma_4 = 

void graph_out(string const &flnm, double ang){
  ofstream dest(flnm.c_str());

  double gamma;
  double a =  1;
  double b =  4;
  int n = 500;
  double d = (b - a) / n;

  MaterialTensor  t = makeTetragonalMaterialTensor(5.6, 5.145, 2.2, 10.6, 2.65, 6.6);
  RotationMatrix  m = multiply(multiply(
					makeAxisRotation(0.00001, 0), 
					makeAxisRotation(0.00001, 1)), 
			       makeAxisRotation(ang,2));


  for(gamma = a ; gamma < b; gamma += d)
    {
      try{
        complex<double> dv = calcDet(t, m, gamma);
        if(true || abs(dv) < 10000)
	  dest << gamma << "\t" << abs(dv) << endl;
      }catch(Error){}
    }    
}

/*
0.01 1.42 1.89 2.18

*/

void test_roots_number(){
    MaterialTensor  t = makeTetragonalMaterialTensor(5.6, 5.145, 2.2, 10.6, 2.65, 6.6);
    RotationMatrix m = multiply(makeAxisRotation(0.1, 0), makeAxisRotation(0.1, 1));

    try{
    double a = 1.0;
    double b = 2.3;
    cout << "one point  : " << a << " ; " << abs(calcDet(t, m, a)) << endl;
    cout << "other point: " << b << " ; " << abs(calcDet(t, m, b)) << endl;
    }catch(Error e)
      {
	cerr << e.what() << endl;
      }
}

void angle_dep(){
  cout << setprecision(10);
  cerr << setprecision(10);

  for(double an = 0; an < 1.57; an += 0.1)
    {
      stringstream flnm;
      flnm << "cut_" << an << ".txt";
      graph_out(flnm.str(), an);
      cerr << an << endl;
    }
}

void test_main(){
  MaterialTensor t = makeTetragonalMaterialTensor(5.6, 5.145, 2.2, 10.6, 6.6, 2.65);

  double x = 1 / sqrt(2);
  double y = x;
  double z = 0;

  PolyMatrix pm = makeChristoffelPolyMatrix(t, x, y, z);

  Poly christ_poly = det(pm);

  Poly::RootVec r = christ_poly.all_roots();
  
  for(int t = 0; t < (int)r.size(); ++t)
    cerr << sqrt(r[t] * 10. / 6.) * 1000. << " ";
  cerr << endl;

}

void test_laguerr(){
  MaterialTensor t = makeTetragonalMaterialTensor(5.6, 5.145, 2.2, 10.6, 6.6, 2.65);

}

int main()
{
  test_laguerr();

  return 0;
}

