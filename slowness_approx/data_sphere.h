#ifndef DATA_SPHERE
#define DATA_SPHERE

#include <Eigen/Core>
#include <vector>
#include <set>

#include "triang_ind.h"

class DataSphere{
 public:
  int n;

  std::vector<Eigen::Vector3d> unitVectors;
  std::vector<TriangInd> triangles;
  std::vector<double> slowness;
  std::vector<std::set<int> > graph;

  DataSphere(int);
  void initSlowness(DataSphere&, DataSphere&, DataSphere&);

  int index_of_pq(int, int)const;
  void paint(double, bool)const;
  void paint_surface(double m = 1)const;
  void paint_wire(double m = 1)const;
  int ic(int a)const{return (a%n);}
  void registerTriangle(const TriangInd&);
  Eigen::Vector3d findNearest(const Eigen::Vector3d &);
};

#endif
