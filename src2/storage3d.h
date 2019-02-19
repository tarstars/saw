#ifndef STORAGE3D
#define STORAGE3D

#include <iostream>
#include <vector>


class Storage3d {
  int h, d, w;
  std::vector<unsigned int> dat;

 public:

 Storage3d() : h(0), d(0), w(0) {}
 Storage3d(int hh, int dd, int ww) : h(hh), d(dd), w(ww), dat(hh * dd * ww) {}
  Storage3d(const char*);
  unsigned int& at(int z, int y, int x) {return dat[z * w * d + y * w + x];}
  const unsigned int& at(int z, int y, int x) const {return dat[z * w * d + y * w + x];}
  void save(const char*);
  static unsigned int getMaxVal();

  friend std::ostream& operator<<(std::ostream&, const Storage3d&);
};


#endif
