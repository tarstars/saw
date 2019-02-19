#ifndef _POLY3_
#define _POLY3_

#include "types.h"

class Vec3;

class CPoly3
{
public:
    CD a3, a2, a1, a0;
    CPoly3(CD ax3 = 0, CD ax2 = 0, CD ax1 = 0, CD ax0 = 0): a3(ax3), a2(ax2), a1(ax1), a0(ax0){}
    CD val(CD x)const {return a0 + x * (a1 + x * (a2 + x * a3));}
    CD der_val(CD x)const {return a1 + x * (CD(2) * a2 + x * CD(3) * a3);}
    CD find_root() const;
    void   find_two_root(CD *, CD *, CD) const; 
};

Vec3 polariz(Vec3 r1, Vec3 r2, CD x);
std::ostream& operator<<(std::ostream& os, const CPoly3& r);

#endif
