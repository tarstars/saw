#include "dx_surface_calculator.h"

#include "util.h"

#include <cmath>
#include <fstream>

using namespace farn;
using namespace std;

DxSurfaceCalculator::DxSurfaceCalculator(int np, const MaterialTensor& mt, double rho) : n(np) {
  #ifdef DEBUG
  cerr << "DxSurfaceCalculator: np = " << np << endl;
  #endif

  dat.push_back(OneDirResult(Vec3(0, 0, 1), mt, rho));

  double dz = 2.0 / n;
  for(int p = 1; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      double z = 1 - dz * p;
      double r = sqrt(1 - z * z);
      double phi = 2 * M_PI / n * q;
      dat.push_back(OneDirResult(Vec3(cos(phi) * r, sin(phi) * r, z), mt, rho));
    }
  }

  dat.push_back(OneDirResult(Vec3(0, 0, -1), mt, rho));
}

void
DxSurfaceCalculator::save(const char* pflnm, int ind) {
  char *ppref = getenv("DXDIR");
  string pref;
  if (ppref) {
    pref = string(ppref);
  }

  #ifdef DEBUG
  cerr << "pref = " << pref << endl;
  #endif

  string fullFilename = pref + pflnm + ".dx";

  #ifdef DEBUG
  cerr << "fullFilename = " << fullFilename << endl;
  #endif

  ofstream dest(fullFilename);

  dest << "object \"slowness positions\" class array type float rank 1 shape 3 items " << dat.size() << " data follows" << endl;
  for(OneDirResult odr : dat) {
    odr.exportPosition(dest, ind);
  }
  dest << endl;

  dest << "object \"slowness data\" class array type float rank 1 shape 3 items " << dat.size() << " data follows" << endl;
  for(OneDirResult odr : dat) {
    odr.exportPolarization(dest, ind);
  }
  dest << endl << endl;

  stringstream bodystream;
  int facesCounter = 0;

  
  for(int q = 0; q < n; ++q) {
    bodystream << 0 << " " << (1 + q) << " " << (1 + (q + 1) % n) << endl;
    ++facesCounter;
  }

  for(int p = 0; p < n - 2; ++p) {
    for(int q = 0; q < n; ++q) {
      bodystream << (1 + q + p * n) << " " << (1 + n + q + p * n) << " " << (1 + (q + 1) % n + p * n) << endl;
      ++facesCounter;
      bodystream << (1 + q + n + p * n) << " " << (1 + n + (q + 1) % n + p * n) << " " << (1 + (q + 1) % n + p * n) << endl;
      ++facesCounter;
    }
  }
  
  int ti = (1 + n * (n - 1));
  for(int q = 0; q < n; ++q) {
    bodystream << ti << " " << (ti - n + (q + 1) % n) << " " << (ti - n + q) << endl;
    ++facesCounter;
  }
  
  dest << "object \"slowness connections\" class array type int rank 1 shape 3 items " << facesCounter << " data follows" << endl;
  dest << bodystream.str();
  dest << "attribute \"element type\" string \"triangles\"" << endl;
  dest << "attribute \"ref\" string \"positions\"" << endl << endl;

  dest << "object \"slowness\" class field" << endl;
  dest << "component \"positions\" \"slowness positions\"" << endl;
  dest << "component \"data\" \"slowness data\"" << endl;
  dest << "component \"connections\" \"slowness connections\"" << endl; 
  dest << endl;
  dest << "end" << endl;
  dest << endl;
}

ostream& farn::operator<<(ostream& os, const DxSurfaceCalculator& r) {
  os << "Surface calculator: " << endl;

  for(OneDirResult val : r.dat) {
    os << val << endl;
  }

  return os;
}
