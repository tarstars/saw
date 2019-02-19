/// @file test_util.h
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#ifndef __TEST_UTIL_H__INCLUDED__
#define __TEST_UTIL_H__INCLUDED__

#include "poly.h"

#include <complex>
#include <limits>
#include <vector>

namespace std {
    bool operator<(const complex<double>& a, const complex<double>& b)
    {
        if(a.real() < b.real())
            return true;
        if(a.real() > b.real())
            return false;

        if(a.imag() < b.imag())
            return true;
        if(a.imag() > b.imag())
            return false;


        return true;
    }
};


namespace farn {

    using namespace std;


    bool eqWithEps(double a, double b, double eps = 8. * numeric_limits<double>::epsilon()) {
        if (a == 0.)
            return abs(b) < eps;
        else
            return abs(a - b) / (abs(a) + abs(b)) < eps;
    }

    bool eqWithEps(complex<double> a, complex<double> b,
                   double eps = 8. * numeric_limits<double>::epsilon()) {
        // return eqWithEps(a.real(), b.real()) && eqWithEps(a.imag(), b.imag());
        if (a == complex<double>(0., 0.))
            return abs(b) < eps;
        else
            return abs(a - b) / (abs(a) + abs(b)) < eps;
    }

    bool polyEqWithEps(Poly const& p1, Poly const& p2) {
        if (p1.c.size() > p2.c.size())
            return false;
        Poly::CoeffVec::const_iterator i, j;
        for (i = p1.c.begin(), j = p2.c.begin(); i != p1.c.end(); ++i, ++j) {
            if (!eqWithEps(*i, *j)) return false;
        }
        for (; j != p2.c.end(); ++j) {
            if (!eqWithEps(0., *j)) return false;
        }
        return true;
    }

    bool rootsEquality(const vector<complex<double> >& a, const vector<complex<double> >& b)
    {
        if(a.size() != b.size()) 
            return false;

        vector<complex<double> > aa(a);
        vector<complex<double> > bb(b);

        sort(aa.begin(), aa.end());
        sort(bb.begin(), bb.end());
        vector<complex<double> >::const_iterator ita, itb;

        for (ita = aa.begin(), itb = bb.begin(); ita != aa.end(); ++ita, ++itb)
            if (!eqWithEps(*ita, *itb)) {
                return false;
            }
        return true;
    }

};

#endif // __TEST_UTIL_H__INCLUDED__
