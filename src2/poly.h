/// @file poly.h
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#ifndef __POLY_H__INCLUDED__
#define __POLY_H__INCLUDED__

#include <complex>
#include <ostream>
#include <vector>

#include <boost/assign/list_of.hpp>

namespace farn {

    using namespace boost;
    using namespace boost::assign;
    using namespace std;

    class Poly
    {
    public:
        typedef complex<double> CoeffType;
        typedef vector<CoeffType> CoeffVec;
        typedef complex<double> RootType;
        typedef vector<RootType> RootVec;

    public:
        Poly() {}
        Poly(CoeffVec const& poly) :
            c(poly.begin(), poly.end())
        {}

        Poly(complex<double> c0) :
            c(1, c0)
        {}

        Poly(complex<double> c0, complex<double> c1) :
            c(list_of(c0)(c1).convert_to_container<Poly::CoeffVec>())
        {}

        Poly(complex<double> c0, complex<double> c1, complex<double> c2) :
            c(list_of(c0)(c1)(c2).convert_to_container<Poly::CoeffVec>())
        {}

        Poly(complex<double> c0, complex<double> c1, complex<double> c2, complex<double> c3) :
            c(list_of(c0)(c1)(c2)(c3).convert_to_container<Poly::CoeffVec>())
        {}

        Poly(complex<double> c0, complex<double> c1, complex<double> c2, complex<double> c3, complex<double> c4) :
            c(list_of(c0)(c1)(c2)(c3)(c4).convert_to_container<Poly::CoeffVec>())
        {}

        Poly operator/(const Poly&) const;
        Poly operator*(const Poly&) const;
        Poly operator+(const Poly&) const;
        Poly operator-(const Poly&) const;

        complex<double> operator()(complex<double> x) const;
  
        template<typename T>
	  void fst(T x, T *pz, T *pf, T *ps, double *) const;

        CoeffVec::size_type deg() const {
            return c.size() - 1;
        }

        RootVec all_roots() const;

        friend ostream& operator<<(ostream&, Poly const&);

    // private:
        CoeffVec c;
    };

    ostream& operator<<(ostream& os, Poly const& p);
    ostream& operator<<(ostream& os, Poly::RootVec const& p);


    inline complex<double> Poly::operator()(complex<double> x) const
    {
        complex<double> ret = 0;
        CoeffVec::const_reverse_iterator t;
        for (t = c.rbegin(); t != c.rend(); ++t) {
            ret = ret * x + (*t);
        }
        return ret;
    }

    //#define _TEST1_

}

#endif // __POLY_H__INCLUDED__
