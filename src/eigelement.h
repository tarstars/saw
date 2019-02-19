#ifndef _EIGELEMENT_
#define _EIGELEMENT_

#include "types.h"

#include "vec3.h"

class CEigElement
{
public:
    CD gam;
    Vec3 vec;
    CD proj;
    CD vel;
    CD slow;

public:
    CEigElement(const Vec3 & xvec, const Vec3 &n, CD xgam):
                                        gam(xgam), 
                                        vec(xvec), 
                                        proj(n * xvec),
                                        vel( sqrt(gam / rho) ),
					  slow(CD(1.0) / vel){}

    bool operator<(const CEigElement& r)const {return abs(proj) < abs(r.proj);}
};

#endif
