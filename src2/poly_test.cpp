/// @file poly_test.cpp
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#include "poly.h"
#include "test_util.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <boost/assign/list_of.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE FarnPolyTest
#include <boost/test/unit_test.hpp>

using namespace boost;
using namespace boost::assign;
using namespace farn;
using namespace std;

BOOST_AUTO_TEST_CASE(test_polyRoots)
{
    
    Poly sam1((Poly::CoeffVec const&)(list_of(-5)(1)));
    Poly sam2((Poly::CoeffVec const&)(list_of(-2)(1)));
    Poly sam = sam1 * sam2;

    vector<complex<double> > allRoots = list_of<complex<double> >(5)(2);
    
    Poly unity((Poly::CoeffVec const&)(list_of(-1)(0)(0)(0)(0)(1)));

    BOOST_CHECK(rootsEquality(allRoots, sam.all_roots()));
}

BOOST_AUTO_TEST_CASE(test_polyOps)
{
    double c11 = 2, c12 = 3, c13 = 5;
    double c21 = 7, c22 = 11, c23 = 13;

    Poly p1((Poly::CoeffVec const&)(list_of(c11)(c12)(c13)));
    Poly p2((Poly::CoeffVec const&)(list_of(c21)(c22)(c23)));
    Poly p3 = list_of(c11 * c21)
        (c11 * c22 + c12 * c21)
        (c11 * c23 + c12 * c22 + c13 * c21)
        (c12 * c23 + c13 * c22)
        (c13 * c23).convert_to_container<Poly::CoeffVec>();
    Poly p4 = list_of(c11 + c21)(c12 + c22)(c13 + c23)
        .convert_to_container<Poly::CoeffVec>();

    BOOST_CHECK(polyEqWithEps(p3, p1 * p2));
    BOOST_CHECK(polyEqWithEps(p2, p3 / p1));
    BOOST_CHECK(polyEqWithEps(p4, p1 + p2));
    BOOST_CHECK(polyEqWithEps(p2, p4 - p1));
    BOOST_CHECK(polyEqWithEps(p2, p2 - Poly()));
}
