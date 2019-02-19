#include <QGLWidget>

#include <cmath>
#include <iostream>

#include <Eigen/Geometry>

#include "data_sphere.h"
#include "util.h"

using namespace Eigen;
using namespace std;
using namespace farn;

DataSphere::DataSphere(int n_):n(n_),unitVectors(n_ * (n_ - 1) + 2){
  double dz = 2. / n;
  unitVectors[0] = Vector3d(0,0,1);
  unitVectors[n * (n - 1) + 1] = Vector3d(0,0, -1);

  graph = vector<set<int> >(unitVectors.size());

  slowness = vector<double> (unitVectors.size(), 1);

  for(int p = 1; p < n; ++p)
    for(int q = 0; q < n; ++q){
      double z = 1 - dz * p;
      double r = sqrt(1 - z * z);
      double phi = M_PI * 2. / n * q;

      Vector3d vec(r * cos(phi), r * sin(phi), z);
      unitVectors[index_of_pq(p,q)] =vec;
    }


  for(int t = 0; t < n ; ++t){
    registerTriangle(TriangInd(index_of_pq(0,0),
				  index_of_pq(1,t),
				  index_of_pq(1,ic(t+1))));
    

    registerTriangle(TriangInd(index_of_pq(n - 1, t),
				  index_of_pq(n,0),
				  index_of_pq(n - 1, ic(t+1))));
    

    for(int p = 1; p < n - 1; ++p){
      registerTriangle(TriangInd(index_of_pq(p,t),
				    index_of_pq(p+1,t),
				    index_of_pq(p+1,ic(t+1))));

      registerTriangle(TriangInd(index_of_pq(p+1,ic(t + 1)),
				    index_of_pq(p,ic(t + 1)),
				    index_of_pq(p,t)));
				    }
    
  }
}

void
DataSphere::registerTriangle(const TriangInd& ti){
  triangles.push_back(ti);
  
  graph[ti.i1].insert(ti.i2);
  graph[ti.i1].insert(ti.i3);
  graph[ti.i2].insert(ti.i1);
  graph[ti.i2].insert(ti.i3);
  graph[ti.i3].insert(ti.i1);
  graph[ti.i3].insert(ti.i2);
}

#ifdef SLOWNESS_WORK
void
DataSphere::initSlowness(DataSphere &ds1, DataSphere &ds2, DataSphere &ds3){
  ds1 = DataSphere(*this);
  ds2 = DataSphere(*this);
  ds3 = DataSphere(*this);

  DataSphere *pds[] = {&ds1, &ds2, &ds3};

  MaterialTensor tensor = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);

  for(int t = 0; t < (int)unitVectors.size(); ++t){
    Vector3d v = unitVectors[t];

    PolyMatrix pm = makeChristoffelPolyMatrix(tensor, v.x(), v.y(), v.z());
    Poly christ_poly = det(pm);
    Poly::RootVec r = christ_poly.all_roots();


    for(int ind = 0; ind < 3; ++ind){
      double gamma = abs(r[ind]);
      double velocity = sqrt(gamma / 5.96e3);
      double m = 1000. / velocity;

      DataSphere * pDest = pds[ind];
      pDest -> slowness[t] = m;
    }
  }
}
#endif 

int
DataSphere::index_of_pq(int p, int q)const{
  if (0 == p)
    return 0;

  if (n == p)
    return n * (n - 1) + 1;

  return 1 + q + (p - 1) * n;
}

void
DataSphere::paint(double m, bool isSurface)const{
  if(isSurface)
    paint_surface(m);
  else
    paint_wire(m);
}

void
DataSphere::paint_surface(double m)const{
  glBegin(GL_TRIANGLES);{
    for(int t = 0; t < (int)triangles.size(); ++t){
      TriangInd const &ti = triangles[t];
      int i1 = ti.i1, i2 = ti.i2, i3 = ti.i3;

      Vector3d n1 = unitVectors[i1];
      Vector3d n2 = unitVectors[i2];
      Vector3d n3 = unitVectors[i3];
      
      Vector3d a = n2 - n1;
      Vector3d b = n3 - n1;
      
      Vector3d no =  -a.cross(b);
      no.normalize();

      glNormal3d(no.x(), no.y(), no.z());

      double m1 = m * slowness[i1];
      double m2 = m * slowness[i2]; 
      double m3 = m * slowness[i3];

      glVertex3d(m1*unitVectors[i1](0),m1*unitVectors[i1](1),m1*unitVectors[i1](2));
      glVertex3d(m2*unitVectors[i2](0),m2*unitVectors[i2](1),m2*unitVectors[i2](2));
      glVertex3d(m3*unitVectors[i3](0),m3*unitVectors[i3](1),m3*unitVectors[i3](2));
    }
  } glEnd();
}

void
DataSphere::paint_wire(double m)const{
  glBegin(GL_LINES);{
    for(int t = 0; t < (int)triangles.size(); ++t){
      TriangInd const &ti = triangles[t];
      int i1 = ti.i1, i2 = ti.i2, i3 = ti.i3;

      Vector3d n1 = unitVectors[i1];
      Vector3d n2 = unitVectors[i2];
      Vector3d n3 = unitVectors[i3];
      
      Vector3d a = n2 - n1;
      Vector3d b = n3 - n1;
      
      Vector3d no =  -a.cross(b);
      no.normalize();

      glNormal3d(no.x(), no.y(), no.z());

      double m1 = m * slowness[i1];
      double m2 = m * slowness[i2]; 
      double m3 = m * slowness[i3];

      glVertex3d(m1*unitVectors[i1](0),m1*unitVectors[i1](1),m1*unitVectors[i1](2));
      glVertex3d(m2*unitVectors[i2](0),m2*unitVectors[i2](1),m2*unitVectors[i2](2));

      glVertex3d(m2*unitVectors[i2](0),m2*unitVectors[i2](1),m2*unitVectors[i2](2));
      glVertex3d(m3*unitVectors[i3](0),m3*unitVectors[i3](1),m3*unitVectors[i3](2));

      glVertex3d(m3*unitVectors[i3](0),m3*unitVectors[i3](1),m3*unitVectors[i3](2));
      glVertex3d(m1*unitVectors[i1](0),m1*unitVectors[i1](1),m1*unitVectors[i1](2));
    }
  } glEnd();
}

Vector3d
DataSphere::findNearest(const Vector3d & aim){
  int index = 0;
  
  bool flag = true;
  bool todo = true;
  double globalMax = 0;
  while(todo){
    todo = false;
    bool localflag = true;
    double localMax = 0;
    int localIndex = 0;
    for(set<int>::iterator it = graph[index].begin(); it!=graph[index].end(); ++it){
      double val = unitVectors[*it].dot(aim);
      if(localflag){
	localflag = false;
	localMax = val;
	localIndex = *it;
      } else {
	if(val > localMax){
	  localMax = val;
	  localIndex = *it;
	}
      }
    }
    if(flag){
      globalMax = localMax;
      flag = false;
      todo = true;
      index = localIndex;
    }else{
      if(localMax > globalMax){
	globalMax = localMax;
	index = localIndex;
	todo = true;
      }
    }
  }
  return unitVectors[index];
}
