#ifndef WORLD
#define WORLD

#include "data_sphere.h"

class GlView;
class MainDialog;

class World{
 public:

  std::vector<DataSphere> spheres;

  World();

  GlView *pView;
  MainDialog * pMainDialog;

  void paint(bool);
  void registerWorld(GlView*);
  void findValue();
};
#endif
