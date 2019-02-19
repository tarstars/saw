#include <iomanip>
#include <iostream>

#include "util.h"
#include "piezo_tensor.h"
#include "vec3.h"

using namespace farn;
using namespace std;

/*
   0    0    0    0   d15  -2d22
-d22  d22    0  d15     0      0
 d31  d31  d33    0     0      0

d15 = 69.2e-12
d22 = 20.8e-12
d31 = -0.85e-12
d33 = 6.0e-12

   0    0    0    0   e15   2e21
 e21 -e21    0  e15     0      0
 e31  e32  e33    0     0      0

Дьелесан:
   0    0    0    0   e15   -e21
-e21  e21    0  e15     0      0
 e31  e32  e33    0     0      0

e11 = 0
e14 = 0
e15 = 3.7
e21 = e22 = 2.5
e31 = e32 = 0.2 
e33 = 1.3

 c11  c12  c13  c14    0   0
 c12  c11  c13 -c14    0   0
 c13  c13  c33    0    0   0
 c14 -c14    0  c44    0   0
   0    0    0    0  c44 c14
   0    0    0    0  c14 (c11 - c12) / 2

c11 = 20.3e10
c12 = 5.3e10
c13 = 7.5e10
c33 = 24.5e10
c44 = 6e10
c14 = 0.9e10
c66 = (c11 - c12) / 2

rho = 4.7e3

eps11 = 44
epss11 = 38.9
eps33 = 29
epss33 = 25.7

eps0 = 8.8542e-12

gamma_i = e_{kij} n_j n_k
eps = epss_{jk} n_j n_k

Gamma_{il} = c_{ijkl} n_j n_k + gamma_i * gamma_l / eps


1999_kushibiki
ce11 = 19.886e10
ce12 =  5.467e10
ce13 =  6.799e10
ce14 =  0.783e10
ce33 = 23.418e10
ce44 =  5.985e10
ce66 =  7.209e10

e15 = 3.655
e22 = 2.407
e31 = 0.328
e33 = 1.894

epss11 = 44.9*eps0
epss33 = 26.7*eps0

rho = 4642.8

LiNbO3 velocities
Ultrasonic transducers for simultaneous generation pf longitudinal and shear waves

36Y-Cut 7340 4001 4084
Z     

  7316 3573 3573
163-Y   6707 4528 3823
x       6572 4795 4079
10Y     7063 4271 4064
 */

void graph(){

  MaterialTensor tens = makeTrigonalMaterialTensor(19.886e10, 5.467e10, 6.799e10, 0.783e10, 23.418e10, 5.985e10, 7.209e10);
  PiezoTensor pieso = makeTrigonal3mMaterialTensor(3.65, 2.407, 0.328, 1.894);

  double rho = 4642.8;
  double eps0 = 8.8542e-12;

  //double phi = iphi / 180. * M_PI;
  //double theta = itheta / 180. * M_PI;

  for(double phi = 0; phi < 2 * M_PI; phi += 0.001){

    Vec3 n(cos(phi), sin(phi), 0);
    //Vec3 n(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    //cout << "n = " << n << endl;

    RotationMatrix christ = makePiezoChristoffelMatrix(tens, pieso, n, 44.9 * eps0, 26.7 * eps0);
    //cout << "christ = " << christ << endl;
    //RotationMatrix christ = makeChristoffelMatrix(tens, n);

    Poly pol = makeCharacterPoly(christ);
    Poly::RootVec rv = pol.all_roots();
    for(int t = 0; t < 3; ++t){
      double v = sqrt(real(rv[t]) / rho);
      //cout << v << " ";
      double s = 1/v;
      cout << s * n.x() << " " << s * n.y() << endl;
    }
  }

  cout << endl;

}

void testDir(int argc, char* argv[]){
  cout << scientific;
  cout.precision(10);

  /*if(argc < 3){
    return;
  }

  stringstream sour1(argv[1]), sour2(argv[2]);
  double iphi;
  double itheta;
  sour1 >> iphi;
  sour2 >> itheta;*/
  MaterialTensor tens = makeTrigonalMaterialTensor(19.886e10, 5.467e10, 6.799e10, 0.783e10, 23.418e10, 5.985e10, 7.209e10);
  PiezoTensor pieso = makeTrigonal3mMaterialTensor(3.65, 2.407, 0.328, 1.894);

  double rho = 4642.8;
  double eps0 = 8.8542e-12;

  double phi = 13 / 180. * M_PI;
  //double theta = itheta / 180. * M_PI;

  //Vec3 n(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  Vec3 n(0, sin(phi), cos(phi));
  cout << "n = " << n << endl;

  RotationMatrix christ = makePiezoChristoffelMatrix(tens, pieso, n, 44.9 * eps0, 26.7 * eps0);
  cout << "christ = " << christ << endl;

  Poly pol = makeCharacterPoly(christ);
  Poly::RootVec rv = pol.all_roots();
  for(int t = 0; t < 3; ++t){
    double v = sqrt(real(rv[t]) / rho);
    cout << v << " ";

    ComplexMatrix cm = makeZeroComplexMatrix();
    for(int p = 0; p < 3; ++p)
      for(int q = 0; q < 3; ++q){
	cm[p][q] = christ[p][q];
	if(p == q)
	  cm[p][q] -= rv[t];
      }

    Vec3 pol = calcPol(cm);
    cout << "pol = " << pol[0] << " " << pol[1] << " " << pol[2] << endl;
  }
  cout << endl;

}

int main(int argc, char* argv[]){
  testDir(argc, argv);
  //graph();
}
