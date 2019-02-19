#include <QGLWidget>
#include <Eigen/Geometry>
#include <iostream>

#include "gl_utils.h"
#include "data_sphere.h"

//using namespace Eigen;
using namespace std;

/*
g++ triang_ind.cpp data_sphere.cpp test_utils.cpp utils.cpp -o test_utils -I/usr/include/qt4/QtOpenGL -I/usr/include/qt4 -L/usr/lib/qt4 -lQtOpenGL -I../src2 ../src2/util.cpp ../src2/poly.cpp
*/

ostream& pp(const Eigen::Vector3d & r){
  return cerr << r.x() << " " << r.y() << " " << r.z();
}

void
draw_sphere(double x, double y, double z, double r, 
	    double cr, double cg, double cb){
  glPushMatrix();{
    glTranslated(x,y,z);
    DataSphere sp(10);
    glPushAttrib(GL_CURRENT_BIT);{
      glColor3d(cr,cg,cb);
      sp.paint_surface(r);
    }glPopAttrib();
  }glPopMatrix();
}

Eigen::Vector3d 
findNormal(const Eigen::Vector3d& r){
  Eigen::Vector3d cand1(r.x(), -r.z(), r.y());
  Eigen::Vector3d cand2(-r.y(), r.x(), r.z());


  Eigen::Vector3d cand = cand1;
  if (cand1.dot(r) > cand2.dot(r))
    cand = cand2;
  
  
  Eigen::Vector3d unitr = r.normalized();
  Eigen::Vector3d ret(cand - unitr * cand.dot(unitr));

  return ret.normalized();
}

void rgvrt(const Eigen::Vector3d &v){
  glVertex3d(v.x(), v.y(), v.z());
}

void rgnorm(const Eigen::Vector3d &v){
  glNormal3d(v.x(), v.y(), v.z());
}

void triangle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3){
  Eigen::Vector3d n = ((v2 - v1).cross(v3 - v1)).normalized();
 
  rgnorm(n);
  rgvrt(v1);
  rgvrt(v2);
  rgvrt(v3);
}

void 
disk(const Eigen::Vector3d& r, const Eigen::Vector3d& norm1, double radius){
  Eigen::Vector3d norm(norm1.normalized());
  Eigen::Vector3d locNormx = findNormal(norm).normalized();
  Eigen::Vector3d locNormy = norm.normalized().cross(locNormx);

  //int error;
  //error = r;

  int n = 30;
  double dphi = 2 * M_PI / n;
  glBegin(GL_TRIANGLES);{
    for(int t = 0; t <= n; ++t){
      double phi1 = t * dphi;
      double phi2 = phi1 + dphi;
      Eigen::Vector3d p1 = r;
      Eigen::Vector3d p2 = r + radius * cos(phi1) * locNormx + radius * sin(phi1) * locNormy;
      Eigen::Vector3d p3 = r + radius * cos(phi2) * locNormx + radius * sin(phi2) * locNormy;

      triangle(p1, p2, p3);
    }
  }glEnd();
}


void cylinder(const Eigen::Vector3d & base, const Eigen::Vector3d &dir, double len, double radius){
  Eigen::Vector3d h0(dir.normalized());
  Eigen::Vector3d h1 = findNormal(h0).normalized();
  Eigen::Vector3d h2 = h0.cross(h1);

  glBegin(GL_TRIANGLES);{
    double dphi = 0.3;
    for(double phi = 0; phi < 6.3; phi += dphi){
      Eigen::Vector3d p1 = base + radius * (h1 * cos(phi) + h2 * sin(phi));
      Eigen::Vector3d p2 = base + radius * ((h1 * cos(phi + dphi) + h2 * sin(phi + dphi)));
      Eigen::Vector3d p3 = base + len * h0 + radius * ((h1 * cos(phi + dphi / 2) + h2 * sin(phi + dphi / 2)));
      Eigen::Vector3d p4 = base + len * h0 + radius * ((h1 * cos(phi + dphi / 2 * 3) + h2 * sin(phi + dphi / 2 * 3)));

      triangle(p1,p2,p3);
      triangle(p3,p2,p4);
    }
  }glEnd();
}

void cone(const Eigen::Vector3d &base, const Eigen::Vector3d &direction,  double radius, double height){
  Eigen::Vector3d h0(direction.normalized());
  Eigen::Vector3d h1(findNormal(direction));
  Eigen::Vector3d h2(h0.cross(h1));

  Eigen::Vector3d p1 = base + direction * height;
  double dphi = 0.3;

  glBegin(GL_TRIANGLES);{
    for(double phi = 0; phi < 6.3; phi += dphi){
      Eigen::Vector3d p2 = base + radius * (h1 * cos(phi) + h2 * sin(phi));
      Eigen::Vector3d p3 = base + radius * (h1 * cos(phi + dphi) + h2 * sin(phi + dphi));

      triangle(p1,p2,p3);
    }
  }glEnd();
}

void ring(const Eigen::Vector3d &base, const Eigen::Vector3d &dir, double in_rad, double out_rad){
  Eigen::Vector3d h0(dir.normalized());
  Eigen::Vector3d h1(findNormal(h0));
  Eigen::Vector3d h2(h0.cross(h1));

  double dphi = 0.3;
  glBegin(GL_TRIANGLES);{
    for(double phi = 0; phi < 6.3; phi += dphi){
      Eigen::Vector3d p11 = base + in_rad * (h1 * cos(phi) + h2 * sin(phi));
      Eigen::Vector3d p12 = base + in_rad * (h1 * cos(phi + dphi) + h2 * sin(phi + dphi));
      Eigen::Vector3d p21 = base + out_rad * (h1 * cos(phi + dphi / 2) + h2 * sin(phi + dphi / 2)); 
      Eigen::Vector3d p22 = base + out_rad * (h1 * cos(phi + dphi / 2 * 3) + h2 * sin(phi + dphi / 2 * 3));

      triangle(p11, p21, p12);
      triangle(p12, p21, p22);
    }
  }glEnd();
}
  
void arrow(const Eigen::Vector3d& r, const Eigen::Vector3d& norm, double tl) {
  double ld = 0.6;
  double thick = 0.1;
  double head = 0.3;

  disk(r, -norm, tl * thick);
  cylinder(r, norm, ld * tl, thick * tl);
  ring(r + norm * ld * tl, -norm, thick * tl, head * tl);
  cone(r + norm * ld * tl, norm, head * tl, (1-ld) * tl);
}
