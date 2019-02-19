#include <QGLWidget>

#include "util_3d.h"

#include <iostream>

using namespace Eigen;
using namespace std;

Vector3d
getNormal(const Vector3d& xdir){
  Vector3d dir(xdir.normalized());

  Vector3d cands[]={Vector3d(-dir.z(), dir.y(), dir.x()),
		    Vector3d(-dir.y(), dir.x(), dir.z()),
		    Vector3d(dir.x(), -dir.z(), dir.y())};

  int ind = 0;
  double minVal = cands[ind].dot(dir);
  
  for(int t = 1; t < 3; ++t){
    double val = cands[t].dot(dir);
    if(val < minVal){
      minVal = val;
      ind = t;
    }
  }

  Vector3d cand = cands[ind];
  Vector3d ret = cand - dir * dir.dot(cand);

  return ret.normalized();
}

ostream& operator<<(ostream& os, Vector3d & r){
  return os << r.x() << " " << r.y() << " " << r.z();
}

void
setCylinder(const Vector3d& xdir,
	    double len,
	    double r){
  Vector3d dir(xdir.normalized());
  Vector3d e1 = getNormal(dir);
  Vector3d e2 = dir.cross(e1);

  int n = 20;
  glBegin(GL_TRIANGLES);{
    for(int t = 0; t < n; ++t){
      double phi11 = M_PI * 2 * t / n;
      double phi12 = M_PI * 2 * (t + 1) / n;
      double phi21 = M_PI * 2 * (t - 1./2) / n;
      double phi22 = M_PI * 2 * (t + 1./2) / n;

      Vector3d r11 = r * (e1 * cos(phi11) + e2 * sin(phi11));
      Vector3d r12 = r * (e1 * cos(phi12) + e2 * sin(phi12));
      Vector3d r21 = r * (e1 * cos(phi21) + e2 * sin(phi21)) + len * dir;
      Vector3d r22 = r * (e1 * cos(phi22) + e2 * sin(phi22)) + len * dir;

      setTriangle(r11, r22, r21);
      setTriangle(r11, r12, r22);
    }
  };glEnd();
}

#define OV(a) (a).x(),(a).y(),(a).z()

void
setTriangle(const Vector3d& a, const Vector3d& b, const Vector3d& c){
  Vector3d r1 = c - a;
  Vector3d norm((c - a).cross(b - a).normalized());
  glNormal3d(OV(norm));

  glVertex3d(OV(a));
  glVertex3d(OV(b));
  glVertex3d(OV(c));
}

void
setRing(const Vector3d& shift,
	const Vector3d& xdir,
	double r1, 
	double r2){
  Vector3d dir(xdir.normalized());
  Vector3d e1 = getNormal(dir);
  Vector3d e2 = dir.cross(e1);

  int n = 20;
  glBegin(GL_TRIANGLES);{
    for(int t = 0; t < n; ++t){
      double phi11 = M_PI * 2 * t / n;
      double phi12 = M_PI * 2 * (t + 1) / n;
      double phi21 = M_PI * 2 * (t - 1./2) / n;
      double phi22 = M_PI * 2 * (t + 1./2) / n;

      Vector3d r11 = r1 * (e1 * cos(phi11) + e2 * sin(phi11)) + shift;
      Vector3d r12 = r1 * (e1 * cos(phi12) + e2 * sin(phi12)) + shift;
      Vector3d r21 = r2 * (e1 * cos(phi21) + e2 * sin(phi21)) + shift;
      Vector3d r22 = r2 * (e1 * cos(phi22) + e2 * sin(phi22)) + shift;

      setTriangle(r11, r22, r21);
      setTriangle(r11, r12, r22);
    }
  }glEnd();  
}

void
setCone(const Vector3d& shift,
	const Vector3d& xdir,
	double h, 
	double r){
  Vector3d dir(xdir.normalized());
  Vector3d e1 = getNormal(dir);
  Vector3d e2 = dir.cross(e1);

  int n = 20;
  glBegin(GL_TRIANGLES);{
    for(int t = 0; t < n; ++t){
      double phi11 = M_PI * 2 * t / n;
      double phi12 = M_PI * 2 * (t + 1) / n;
 
      Vector3d r11 = r * (e1 * cos(phi11) + e2 * sin(phi11)) + shift;
      Vector3d r12 = r * (e1 * cos(phi12) + e2 * sin(phi12)) + shift;
      Vector3d r21 = shift + h * dir;
 
      setTriangle(r11, r12, r21);
    }
  }glEnd();  
}
