#include "height_calculator.h"
#include "util.h"

#include <iostream>
#include <fstream>

#define DEBUG

#include <cmath>

using namespace farn;
using namespace std;

const double rho = 5.96e3;

HeightCalculator::HeightCalculator():
  ten(makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10)){
  
}

vector<pair<double, double> >
HeightCalculator::getAllHeights(double x1, double y1, double step) {
	cout << "getAllHeights: x1 = " << x1 << ",y1=" << y1 << 
		",step=" << step << endl;
	double dx = step*x1/hypot(x1,y1);
	double dy = step*y1/hypot(x1,y1);
	vector<pair<double, double> > retVal;

	for (double x = 0., y = 0.; abs(x) < abs(x1); x+=dx, y+=dy) {
		PolyMatrix polyMat = makePolyMatrix(ten, rho, x, y);
		Poly pol = det(polyMat);
		Poly::RootVec roots = pol.all_roots();
		for (unsigned int i = 0; i < roots.size(); ++i) {
    		if(abs(roots[i].imag()) < 1e-6 && roots[i].real() > 0){
      			//dat.push_back(roots[t].real());
				retVal.push_back(make_pair(hypot(x,y), roots[i].real()));
    		}
			
		}
	}
	cout << "getAllHeights: result size: " << retVal.size() << endl;
	return retVal;
}

double
HeightCalculator::getHeight(double x, double y, int modeCode){
  static ofstream heightLog("hl.log"); //use cerr instead in linux
  #ifdef DEBUG
  heightLog << "getHeight(" << x << ", " << y << ", " << modeCode << ")" << endl;
  #endif

  PolyMatrix polyMat = makePolyMatrix(ten, rho, x, y);
  Poly pol = det(polyMat);

  #ifdef DEBUG
  heightLog << "\tpolynome is " << pol << endl;
  #endif

  Poly::RootVec roots = pol.all_roots();

  #ifdef DEBUG
  heightLog << "\tall roots: ";
  for(int t = 0; t < int(roots.size()); ++t)
    heightLog << roots[t] << " ";
  heightLog << endl;
  #endif

  vector<double> dat;
  for(int t = 0; t < int(roots.size()); ++t){
    if(abs(roots[t].imag()) < 1e-6 && roots[t].real() > 0){
      dat.push_back(roots[t].real());
    }
  }

  sort(dat.begin(), dat.end());

  #ifdef DEBUG
  heightLog << "\tselected roots: ";
  for(int t = 0; t < int(dat.size()); ++t){
    heightLog << dat[t] << " ";
  }
  heightLog << endl << endl;
  #endif

  if (int(dat.size()) > modeCode){
    double val = dat[dat.size() - 1 - modeCode];
    return val;
  } 
  
  return 0;
}


