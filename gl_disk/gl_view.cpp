#include "gl_utils.h"
#include "gl_view.h"

#include <cmath>
#include <Eigen/Geometry>
#include <iostream>

using namespace std;
using namespace Eigen;

GlView::GlView(QWidget* parent):
  QGLWidget(parent){

}

void
GlView::initializeGL(){
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glEnable(GL_SMOOTH);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_LIGHTING);

  {
    glEnable(GL_LIGHT0);
    
    GLfloat ambientLight[] = {0, 0, 0, 1};
    GLfloat diffuseReflection[] = {.2, .7, .7, 1};
    GLfloat specularLight[] = {.5, .5, .5, 1.0};
    GLfloat pos[]={2, 2, 8, 1};
    GLfloat dir[] = {0, 0, -1};

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);

  }
}

void
GlView::resizeGL(int w, int h){
  int s = qMin(w, h);
  glViewport((w - h) / 2, (h - s) / 2, s, s);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-.5, .5, -0.5, 0.5, 4, 12);
  glMatrixMode(GL_MODELVIEW);
}

void
setTriangle(const Vector3d& r1, const Vector3d& r2, const Vector3d& r3) {
  Vector3d v1 = (r2 - r1);
  Vector3d v2 = (r3 - r1);
  Vector3d norm = v1.cross(v2).normalized();

  glNormal3d(norm.x(), norm.y(), norm.z());
  
  glVertex3d(r1.x(), r1.y(), r1.z());
  glVertex3d(r2.x(), r2.y(), r2.z());
  glVertex3d(r3.x(), r3.y(), r3.z());
}


void
GlView::paintGL(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  double rtheta = theta * M_PI / 180;
  double rphi = phi * M_PI / 180;

  Vector3d dir(cos(rtheta) * cos(rphi), cos(rtheta) * sin(rphi), sin(rtheta));
  Vector3d dir_norm(cos(rtheta + M_PI / 2) * cos(rphi), cos(rtheta + M_PI / 2) * sin(rphi), sin(rtheta + M_PI / 2));

  glLoadIdentity();
  glTranslated(0, 0, -7);
  glRotated(rotx, 1.0, 0.0, 0.0);
  glRotated(roty, 0.0, 1.0, 0.0);
  glRotated(rotz, 0.0, 0.0, 1.0);

  /*  
    glColor3f(0, 1, 1);
    ring(dir * 0.3, dir, 0.1, 0.25);

    glColor3f(1, 1, 1);
    disk(dir * 0.2, dir, 0.2);*/

  arrow(dir * 0.3, dir, 0.5);

  glPushAttrib(GL_CURRENT_BIT);
  glBegin(GL_LINES);{
    glNormal3d(0,0,1);
    glColor3f(3.0, .0, 0.0);
    glVertex3f(-0.3, 0.0, 0.0);
    glVertex3f( 1.0, 0.0, 0.0);

    glNormal3d(1.0, 0, 0);
    glColor3f(0.0, 3.0, 0.0);
    glVertex3f(0.0, -0.3, 0.0);
    glVertex3f(0.0,  1.0, 0.0);

    glNormal3d(0, 1, 0);
    glColor3f(0.0, 0.0, 3.0);
    glVertex3f( 0.0, 0.0, -0.3);
    glVertex3f( 0.0, 0.0,  1.0);

    glNormal3d(dir_norm.x(), dir_norm.y(), dir_norm.z());
    glColor3f(3.0, 3.0, 3.0);
    glVertex3f( 0.0, 0.0, 0.0);
    glVertex3f( dir.x(), dir.y(),  dir.z());
  };glEnd();
  glPopAttrib();
}

void
GlView::setPhi(int phi){
  this->phi = phi;
  updateGL();
}

void
GlView::setTheta(int theta){
  this->theta = theta;
  updateGL();
}


void
GlView::setRotx(int rotx){
  this->rotx = rotx;
  updateGL();
}

void
GlView::setRoty(int roty){
  this->roty = roty;
  updateGL();
}

void
GlView::setRotz(int rotz){
  this->rotz = rotz;
  updateGL();
}
