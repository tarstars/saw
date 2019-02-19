#include "storage3d.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace std;

void
test_print_file(const char* flnm) {
  Storage3d dat(flnm);
  cout << dat << endl;
}

int
test_copy(const char* pfrom, const char* pto) {
  Storage3d dat(pfrom);
  dat.save(pto);

  return 0;
}

void
test_fill(const char* flnm) {
  int h=1;
  int d=50;
  int w=50;

  double sigma = max(max(h, d), w)/ 3;
  double number = 10;

  Storage3d dat(h, d, w);

  for(int z = 0; z < h; ++z) {
    for(int y = 0; y < d; ++y) {
      for(int x = 0; x < w; ++x) {
	double val = exp(-((x-w/2)*(x-w/2)+(y-d/2)*(y-d/2)+(z-h/2)*(z-h/2)) / (sigma * sigma));
	//double val = double(x) / w;
	dat.at(z, y, x) = Storage3d::getMaxVal() * val;
      }
    }
  }

  dat.save(flnm);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    cout << "usage: cpp_conv command ..." << endl;
    cout << "command = create => cpp_conv create flnm" << endl;
    cout << "command = print => cpp_conv print flnm" << endl;
    cout << "command = copy => cpp_conv copy flnm1 flnm2" << endl;
    return 1;
  }

  string mode(argv[1]);
  if (mode == "create") {
    test_fill(argv[2]);
    cout << "creation done" << endl;
    return 0;
  }

  if (mode == "print") {
    test_print_file(argv[2]);
    cout << "printing is done" << endl;
    return 0;
  }

  if (mode == "copy") {
    test_copy(argv[2], argv[3]);
    cout << "copied" << endl;
  }

  return 0;
}
