#include "storage3d.h"


#include <fstream>
#include <limits>

using namespace std;

unsigned int read_integer(istream& is, int sz) {
  unsigned int ret = 0;
  
  for(int sh = sz - 1; sh >= 0; --sh) {
    is.read((char*)&ret + sh, 1);
  }

  return ret;
}

void write_integer(ostream& is, unsigned int val, int sz) {
  for(int sh = sz - 1; sh >= 0; --sh) {
    is.write((char*)&val + sh, 1);
  }
}

int read_dimension(istream& is) {
  return read_integer(is, 2);
}

Storage3d::Storage3d(const char* flnm) {
  ifstream is;
  is.open(flnm, ios::binary);
  
  w = read_dimension(is);
  d = read_dimension(is);
  h = read_dimension(is);

  dat.resize(w * d * h);

  streampos cp = is.tellg();
  is.seekg(0, ios::end);
  streampos ep = is.tellg();
  is.seekg(cp);

  int sz = (int)(ep - cp);
  int ds = sz / (w * d * h);

  for(int z = 0; z < h; ++z) {
    for(int y = 0; y < d; ++y) {
      for(int x = 0; x < w; ++x) {
	at(z, y, x) = read_integer(is, ds);
      }
    }
  }
}

ostream& operator<<(ostream& os, const Storage3d& r) {
  for(int z = 0; z < r.h; ++z) {
    for(int y = 0; y < r.w; ++y) {
      for(int x = 0; x < r.d; ++x) {
	os << r.at(z, y, x) << " ";
      }
      os << endl;
    }
    os << endl;
  }

  return os;
}

void
Storage3d::save(const char* flnm) {
  ofstream os;
  os.open(flnm, ios::binary);

  write_integer(os, w, 2);
  write_integer(os, d, 2);
  write_integer(os, h, 2);

  for(int z = 0; z < h; ++z) {
    for(int y = 0; y < d; ++y) {
      for(int x = 0; x < w; ++x) {
	write_integer(os, at(z, y, x), sizeof(int));
      }
    }
  }
}

unsigned int 
Storage3d::getMaxVal() {
  return numeric_limits<unsigned int>::max();
}
