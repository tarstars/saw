#include <QGLWidget>

#include <iostream>

#include "gl_view.h"
#include "main_dialog.h"
#include "world.h"
#include "gl_utils.h"

using namespace std;
using namespace Eigen;

World::World(){
  DataSphere sphere(150);
  for(int t = 0; t < 3; ++t){
    spheres.push_back(DataSphere(0));
  }

  sphere.initSlowness(spheres[0], spheres[1], spheres[2]);
}

void
World::paint(bool isSurface){
  glPushAttrib(GL_CURRENT_BIT);{
    glColor3d(1.0, 0.5, 0.0);
    if(pMainDialog->rbSelectSurface1->isChecked()){
      spheres[0].paint(1,isSurface);
    }
    
    glColor3d(0.3, 1.0, 0.0);
    if(pMainDialog->rbSelectSurface2->isChecked())
      spheres[1].paint(1,isSurface);

    glColor3d(0.0, 0.6, 1.0);
    if(pMainDialog->rbSelectSurface3->isChecked())
      spheres[2].paint(1,isSurface);
  }glPopAttrib();


  
  // glPushAttrib(GL_CURRENT_BIT);
  // glBegin(GL_QUADS);{
  //   glColor3f(1.0, 1.0, 0.);

  //   glNormal3d(0,1,0);
  //   glVertex3f(-.5, 0, -.5);
  //   glVertex3f( .5, 0, -.5);
  //   glVertex3f( .5, 0,  .5);
  //   glVertex3f(-.5, 0,  .5);

  //   glNormal3d(0,-1,0);
  //   glVertex3f(-.5, 0, -.5);
  //   glVertex3f(-.5, 0,  .5);
  //   glVertex3f( .5, 0,  .5);
  //   glVertex3f( .5, 0, -.5);

  // }glEnd();
  // glPopAttrib();  
  
}

void
World::registerWorld(GlView* pView_){
  pView = pView_;
  pView -> pWorld = this;
}

void
World::findValue(){
  Vector3d direction = pMainDialog -> getDirection();
  int as = pMainDialog -> getActiveSurface();

  Vector3d nearest = spheres[as].findNearest(direction);
}
