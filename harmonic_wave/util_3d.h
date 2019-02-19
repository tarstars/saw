#ifndef UTILS
#define UTILS

#include <eigen3/Eigen/Geometry>

Eigen::Vector3d getNormal(const Eigen::Vector3d& xdir);

void setCylinder(const Eigen::Vector3d& xdir,
		 double len,
		 double r);

void setTriangle(const Eigen::Vector3d& a, 
		 const Eigen::Vector3d& b, 
		 const Eigen::Vector3d& c);

void setRing(const Eigen::Vector3d&,
	     const Eigen::Vector3d&, 
	     double, 
	     double);


void setCone(const Eigen::Vector3d&, 
	     const Eigen::Vector3d&,
	     double, 
	     double);
#endif
