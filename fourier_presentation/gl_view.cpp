#include <iostream>
#include <cmath>
#include <Eigen/Geometry>

#include "gl_view.h"

using namespace std;
using namespace Eigen;

GlView::GlView(QWidget* parent):QGLWidget(parent){
  xang = 0;
  yang = 0;
  zang = 0;
}

void
GlView::initializeGL(){
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glEnable(GL_SMOOTH);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);

  GLfloat specularMat[] = {0,0,1,1};
  GLfloat emissionMat[] = {0,0,0,1};

  glMaterialfv(GL_FRONT, GL_SPECULAR, specularMat);
  glMaterialfv(GL_FRONT, GL_EMISSION, emissionMat);

  glEnable(GL_LIGHT0);
  {
    GLfloat ambientLight[] = {0.2, 0.2, 0.2, 1};
    GLfloat diffuseReflection[] = {0.5, 0.5, 0.5, 1};
    GLfloat specularLight[] = {0.5, 0.5, 0.5, 1};

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  }
}

void
GlView::resizeGL(int w, int h){
  int s = min(w,h);
  glViewport(( w - s) / 2, (w - s) / 2, s, s); 

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-1, 1, -1, 1, 4, 12);

  glMatrixMode(GL_MODELVIEW);

}

// double f(double x, double z){
//   double r = 1 - (x * x + z * z);
//   if (r < 0) throw(string("no point"));
//   return sqrt(r);
// }

double f(double x, double y){
  double r = sqrt(x * x + y * y);
  return cos(r * 10) * (1 - r);
}

void
glTriangle(const Vector3d& r1, const Vector3d& r2, const Vector3d& r3){
  Vector3d norm = (r2 - r1).cross(r3 - r1).normalized();
  
  glNormal3d(norm.x(), norm.y(), norm.z());
  glVertex3d(r1.x(), r1.y(), r1.z());
  glVertex3d(r2.x(), r2.y(), r2.z());
  glVertex3d(r3.x(), r3.y(), r3.z());
}

void
GlView::paintGL(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glColor3d(0,1,0);
  
  glLoadIdentity();

  glTranslated(0, 0, -7);

  glRotated(xang, 1, 0, 0);
  glRotated(yang, 0, 1, 0);
  glRotated(zang, 0, 0, 1);

  glBegin(GL_TRIANGLES);{

    double dx = 0.02, dz = 0.02;

    for(double x = -1; x < 1; x += dx)
      for(double z = -1; z < 1; z += dz){
	try{
	  double y1 = f(x, z);
	  double y2 = f(x + dx, z);
	  double y3 = f(x, z + dz);
	  double y4 = f(x + dx, z + dz);

	  Vector3d v11(x, y1, z);
	  Vector3d v12(x + dx, y2, z);
	  Vector3d v13(x, y3, z + dz);

	  Vector3d v21(x + dx, y2, z);
	  Vector3d v22(x + dx, y4, z + dz);
	  Vector3d v23(x, y3, z + dz);

	  glTriangle(v11, v12, v13);
	  glTriangle(v21, v22, v23);
	}catch(...){};
      }
  }glEnd();
}

void
GlView::setRotx(int xang){
  this -> xang  = xang;
  updateGL();
}

void
GlView::setRoty(int yang){
  this -> yang  = yang;
  updateGL();
}

void
GlView::setRotz(int zang){
  this -> zang  = zang;
  updateGL();
}

