#ifndef UTILS
#define UTILS

#include <Eigen/Geometry> 

void draw_sphere(double x, double y, double z, double r=1, 
	    double cr=1, double cg=1, double cb=1);

void cylinder(const Eigen::Vector3d &base, 
	      const Eigen::Vector3d &dir, 
	      double len = 5, 
	      double radius = 1.0);

void cone(const Eigen::Vector3d &base, 
	  const Eigen::Vector3d &direction,  
	  double radius, 
	  double height);

void disk(const Eigen::Vector3d& r, 
	  const Eigen::Vector3d& norm, 
	  double radius);

void ring(const Eigen::Vector3d &base,
	  const Eigen::Vector3d &dir,
	  double in_rad,
	  double out_rad);

void arrow(const Eigen::Vector3d& r, 
	   const Eigen::Vector3d& norm, 
	   double tl);

#endif
