#include <cmath>
#include <iostream>

#include <eigen3/Eigen/Geometry>

#include "surface_view.h"
#include "util_3d.h"

using namespace std;
using namespace Eigen;

SurfaceView::SurfaceView(QWidget* parent):
  QGLWidget(parent){

}

void
SurfaceView::initializeGL(){
  glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE);
  glEnable(GL_BLEND);

  glEnable(GL_SMOOTH);
  
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_LIGHTING);


  /*GLfloat specularMat[] = {0.8, 0.8, 0.8, 1.0};
  GLfloat emissionMat[] = {0.3, 0.3, 0.3, 1.0};

  glMaterialfv(GL_FRONT, GL_SPECULAR, specularMat);
  glMaterialfv(GL_FRONT, GL_EMISSION, emissionMat);
  */
  {
    glEnable(GL_LIGHT0);

    GLfloat ambientLight[] = {0.3, 0.3, 0.3, 1.0};
    GLfloat diffuseReflection[] = {0.3, 0.3, 0.3, 1.0};
    GLfloat specularLight[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat pos[] = {0., 3., -20., 1};
    GLfloat dir[] = {0., 0., -1.};

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
  }

}

void
SurfaceView::resizeGL(int w, int h){
  int s = qMin(w, h);

  glViewport((w - s) / 2, (h - s) / 2, s, s);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-1.3, 1.3, -1.3, 1.3, 6 , 14);

  glMatrixMode(GL_MODELVIEW);
}

Vector3d
sphere2vector(double theta, double phi){
  return Vector3d(cos(phi) * sin(theta), 
		  sin(phi) * sin(theta),
		  cos(theta));
}


Vector3d
pqPoint(int p, int q, int n){
  double theta = acos(-1 + p * 2. / n);
  double phi = 2 * M_PI / n * q;
  return sphere2vector(theta, phi);
}

void
setSphere(double r){
  int n = 20;

  glBegin(GL_TRIANGLES);{

    for(int p = 1; p < n - 1; ++p)
      for(int q = 0; q < n; ++q){
	setTriangle(r * pqPoint(p + 1, q, n), 
		    r * pqPoint(p, q, n), 
		    r * pqPoint(p + 1, q + 1, n));

	setTriangle(r * pqPoint(p, q, n), 
		    r * pqPoint(p, q + 1, n), 
		    r * pqPoint(p + 1, q + 1, n));

      }

    for(int q = 0; q < n; q++){
      setTriangle(r * pqPoint(0, 0, n),
		  r * pqPoint(1, q + 1, n),
		  r * pqPoint(1, q, n));
  
      setTriangle(r * pqPoint(n - 1, q, n),
		  r * pqPoint(n - 1, q + 1, n),
		  r * pqPoint(n, 0, n));
    }
  }glEnd();
}



double
spc(int ind, double A, int n){
  return -A + 2 * A * ind / (n - 1);
}

double
f(double x, double y, const Vector3d& dir){
  Vector3d cur(x, y, 0);
  double d = cur.dot(dir);
  return  0.2 * sin(3 * 2 * M_PI * d); 
}

void
setzxy(const Vector3d& dir){
  int n = 100;
  double A = 1;
  glBegin(GL_TRIANGLES);{
    for(int p = 0; p < n; ++p)
      for(int q = 0; q < n; ++q){
	double cx = spc(q, A, n);
	double px = spc(q + 1, A, n);
	double cy = spc(p, A, n);
	double py = spc(p + 1, A, n);

	double zcc = f(cx, cy, dir);
	double zpc = f(px, cy, dir);
	double zcp = f(cx, py, dir);
	double zpp = f(px, py, dir);

	if( (cx * cx + cy * cy < 1) &&
	    (px * px + py * py < 1)){
	  setTriangle(Vector3d(cx, cy, zcc),
		      Vector3d(px, cy, zpc), 
		      Vector3d(px, py, zpp));
 
	  setTriangle(Vector3d(cx, cy, zcc),
		      Vector3d(px, py, zpp), 
		      Vector3d(cx, py, zcp));
	}
     
      }
  };glEnd();
}

void
SurfaceView::paintGL(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity();
  glTranslated(0, 0, -10);

  //rotate scene
  glRotated(rotx, 1, 0, 0);
  glRotated(roty, 0, 1, 0);
  glRotated(rotz, 0, 0, 1);

  //coordinate system
  glColor3d(3, 0, 0); setCylinder(Vector3d(1, 0, 0), 1.5, 0.02);
  glColor3d(0, 3, 0); setCylinder(Vector3d(0, 1, 0), 1.5, 0.02);
  glColor3d(0, 0, 3); setCylinder(Vector3d(0, 0, 1), 1.5, 0.02);

  glColor3d(0.6, 0.6, 0.6);
  setSphere(0.25);

  glPushMatrix();{
    glTranslated(1.5, 0, 0);
    glColor3d(0.9, 0.4, 0.4);
    setSphere(0.25);
  }glPopMatrix();

  glPushMatrix();{
    glTranslated(0, 1.5, 0);
    glColor3d(0.4, 0.9, 0.4);
    setSphere(0.25);
  }glPopMatrix();

  glPushMatrix();{
    glTranslated(0, 0, 1.5);
    glColor3d(0.4, 0.4, 0.9);
    setSphere(0.25);
  }glPopMatrix();

  //surface 
  glColor3d(1.1, 1.1, 0);
  setzxy(dir);

  //direction arrow
  glColor3d(1.0, 0.0, 0);
  setCylinder(dir, .9, 0.025);
  setRing(dir * .9, dir, 0.1, .025);
  setCone(dir * .9, dir, 0.1, 0.1);

  //vertical projection
  Vector3d e3 = Vector3d(0,0,1).normalized();
  Vector3d pos = e3 * e3.dot(dir);
  glColor3d(0, 1, 1);
  setRing(pos, e3, 0, 0.5);
  
  
}

void
SurfaceView::onTimer(){
  update();
}

void
SurfaceView::setTheta(int xTheta){
  theta = 2 * M_PI * xTheta / 360;
  dir = sphere2vector(theta, phi);
  update();
}

void
SurfaceView::setPhi(int xPhi){
  phi = 2 * M_PI * xPhi / 360;
  dir = sphere2vector(theta, phi);
  update();
}

void
SurfaceView::setRotx(int xRotx){
  rotx = xRotx;
  update();
}

void
SurfaceView::setRoty(int xRoty){
  roty = xRoty;
  update();
}

void
SurfaceView::setRotz(int xRotz){
  rotz = xRotz;
  update();
}

void
SurfaceView::showAngles(){
  cerr << "rotx roty rotz = " << rotx << " " << roty << " " << rotz << endl;
}
