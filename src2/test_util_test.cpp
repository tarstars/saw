/// @file test_util_test.cpp
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#include "test_util.h"

#include <cmath>
#include <iostream>
#include <limits>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE FarnPolyTest
#include <boost/test/unit_test.hpp>

// using namespace boost;
using namespace farn;
using namespace std;

BOOST_AUTO_TEST_CASE(test_eqWithEps)
{
    BOOST_CHECK(eqWithEps(1 + 1e-15, 1 + 1e-16));
    BOOST_CHECK(!eqWithEps(1 + 1e-14, 1 + 1e-15));
    BOOST_CHECK(eqWithEps(M_PI, 3.1415926535897932384));

    BOOST_CHECK(eqWithEps(complex<double>(1, 1e-15), complex<double>(1, 1e-16)));
    BOOST_CHECK(!eqWithEps(complex<double>(1, 1e-14), complex<double>(1, 1e-15)));
}
