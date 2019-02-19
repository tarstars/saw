/// @file util_test.cpp
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>


#include "test_util.h"
#include "util.h"

#include <iostream>
#include <limits>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE farn
#include <boost/test/unit_test.hpp>

using namespace farn;
using namespace std;

BOOST_AUTO_TEST_CASE( test_tensorPrint )
{
    using namespace farn;
    using namespace std;

    typedef boost::multi_array<double, 4> array_type;
    boost::array<array_type::index, 4> shape = {{3, 3, 3, 3}};
    array_type tensor(shape);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    tensor[i][j][k][l] = l + 10 * (k + 10 * (j + 10 * i)) + 10000;
                }
            }
        }
    }

    ostringstream actual;
    actual << tensor;

    char const* expected =
        "10000 10001 10002  | 10100 10101 10102  | 10200 10201 10202  | \n"
        "10010 10011 10012  | 10110 10111 10112  | 10210 10211 10212  | \n"
        "10020 10021 10022  | 10120 10121 10122  | 10220 10221 10222  | \n\n"
        "11000 11001 11002  | 11100 11101 11102  | 11200 11201 11202  | \n"
        "11010 11011 11012  | 11110 11111 11112  | 11210 11211 11212  | \n"
        "11020 11021 11022  | 11120 11121 11122  | 11220 11221 11222  | \n\n"
        "12000 12001 12002  | 12100 12101 12102  | 12200 12201 12202  | \n"
        "12010 12011 12012  | 12110 12111 12112  | 12210 12211 12212  | \n"
        "12020 12021 12022  | 12120 12121 12122  | 12220 12221 12222  | \n\n";

    BOOST_CHECK_EQUAL(actual.str(), string(expected));
}

BOOST_AUTO_TEST_CASE( test_makeTetragonalMaterialTensor )
{
    MaterialTensor t = makeTetragonalMaterialTensor(1, 2, 3, 4, 5, 6);
    
    int i, j, k, l;

    for (i = 0; i < 4; ++i)
        BOOST_CHECK_EQUAL(3, t.shape()[i]);

    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            for (k = 0; k < 3; ++k)
                for (l = 0; l < 3; ++l) {
                    BOOST_CHECK_EQUAL(t[i][j][k][l], t[j][i][k][l]);
                    BOOST_CHECK_EQUAL(t[i][j][k][l], t[i][j][l][k]);
                    BOOST_CHECK_EQUAL(t[i][j][k][l], t[k][l][i][j]);
                }
    BOOST_CHECK_EQUAL(1, t[0][0][0][0]);
    BOOST_CHECK_EQUAL(1, t[1][1][1][1]);
    BOOST_CHECK_EQUAL(2, t[0][0][1][1]);
    BOOST_CHECK_EQUAL(3, t[0][0][2][2]);
    BOOST_CHECK_EQUAL(3, t[1][1][2][2]);
    BOOST_CHECK_EQUAL(4, t[2][2][2][2]);
    BOOST_CHECK_EQUAL(5, t[1][2][1][2]);
    BOOST_CHECK_EQUAL(5, t[0][2][0][2]);
    BOOST_CHECK_EQUAL(6, t[0][1][0][1]);
}

bool matrixEq(const RotationMatrix &a /*expected*/, const RotationMatrix& b /*actual*/) {
    double eps = 8 * numeric_limits<double>::epsilon();
    int i, j;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            if (a[i][j] == 0.) {
                if (fabs(b[i][j]) >= eps) return false;
            } else {
                if (fabs(a[i][j] - b[i][j]) / (fabs(a[i][j]) + fabs(b[i][j]))  > eps)
                    return false;
            }

    return true;
}

bool tensorEq(const MaterialTensor &a /*expected*/, const MaterialTensor& b /*actual*/) {
    double eps = 8 * numeric_limits<double>::epsilon();

    MaterialTensor::const_iterator ita, itb;

    int i,j,k,l;
    for(i = 0; i < 3; ++i)
        for(j = 0; j < 3; ++j)
            for( k = 0; k < 3; ++k)
                for( l = 0; l < 3; ++l)
                {
                    double vala = a[i][j][k][l];
                    double valb = b[i][j][k][l];
                    if (vala == 0.) {
                        if (fabs(valb) >= eps) return false;
                    } else {
                        if (fabs(vala - valb) / (fabs(vala) + fabs(valb))  > eps)
                            return false;
                    }
                }

    return true;
}

BOOST_AUTO_TEST_CASE( test_makeAxisRotation ) {
    {  
        RotationMatrix mat1 = makeAxisRotation(0.3, 0);
        RotationMatrix mat2 = makeAxisRotation(0.6, 0);
        RotationMatrix mat1inv = makeAxisRotation(-0.3, 0);
        RotationMatrix matEq = makeAxisRotation(0, 0);

        RotationMatrix mat2alt = multiply(mat1, mat1);

        BOOST_CHECK(matrixEq(mat2, mat2alt));
        BOOST_CHECK(matrixEq(matEq, multiply(mat1, mat1inv)));
    }

    {
        RotationMatrix mat1 = makeAxisRotation(0.3, 1);
        RotationMatrix mat2 = makeAxisRotation(0.6, 1);
        RotationMatrix mat1inv = makeAxisRotation(-0.3, 1);
        RotationMatrix matEq = makeAxisRotation(0, 1);
    

        RotationMatrix mat2alt = multiply(mat1, mat1);

        BOOST_CHECK(matrixEq(mat2, mat2alt));
        BOOST_CHECK(matrixEq(matEq, multiply(mat1, mat1inv)));
    }

    {
        RotationMatrix mat1 = makeAxisRotation(0.3, 2);
        RotationMatrix mat2 = makeAxisRotation(0.6, 2);
        RotationMatrix mat1inv = makeAxisRotation(-0.3, 2);
        RotationMatrix matEq = makeAxisRotation(0, 2);
    

        RotationMatrix mat2alt = multiply(mat1, mat1);

        BOOST_CHECK(matrixEq(mat2, mat2alt));
        BOOST_CHECK(matrixEq(matEq, multiply(mat1, mat1inv)));
    }

    {
        RotationMatrix mat1 = makeAxisRotation(0.3, 0);
        RotationMatrix mat2 = makeAxisRotation(0.3, 1);
        RotationMatrix mat3 = makeAxisRotation(0.3, 2);
        RotationMatrix mat1inv = makeAxisRotation(-0.3, 0);
        RotationMatrix mat2inv = makeAxisRotation(-0.3, 1);
        RotationMatrix matEq = makeAxisRotation(0, 0);

        BOOST_CHECK(matrixEq(
                             multiply(multiply(mat1, mat2), mat3),
                             multiply(mat1, multiply(mat2, mat3))));

        BOOST_CHECK(matrixEq(matEq, multiply(
                                             multiply(mat1, mat2),
                                             multiply(mat2inv, mat1inv))));
    }
}

BOOST_AUTO_TEST_CASE(test_TensorRotation)
{
    RotationMatrix mat1 = makeAxisRotation(0.3, 0);
    RotationMatrix mat2 = makeAxisRotation(0.6, 0);
    RotationMatrix mat1inv = makeAxisRotation(-0.3, 0);

    MaterialTensor tensor = makeTetragonalMaterialTensor(1, 2, 3, 4, 5, 6);

    BOOST_CHECK(tensorEq(tensor, rotateTensor(rotateTensor(tensor, mat1), mat1inv)));
    BOOST_CHECK(tensorEq(rotateTensor(tensor, mat2), rotateTensor(rotateTensor(tensor, mat1), mat1)));
}


BOOST_AUTO_TEST_CASE( test_polyMatrix )
{

    double c11 = 2, c12 = 3, c13 = 5, c33 = 7, c44 = 11, c66 = 13;
    MaterialTensor t = makeTetragonalMaterialTensor(c11, c12, c13, c33, c44, c66);
    double n1 = 17, n2 = 19, n3 = 23;

    PolyMatrix tstMatrix = makePolyMatrix(t, 0, n1, n2);
    ComplexMatrix numbersMatrix = evaluatePolyMatrix(tstMatrix, n3);

    BOOST_CHECK(eqWithEps(c11 * n1 * n1 + c66 * n2 * n2 + c44 * n3 * n3, numbersMatrix[0][0]));
    BOOST_CHECK(eqWithEps((c12 + c66) * n1 * n2, numbersMatrix[0][1]));
    BOOST_CHECK(eqWithEps((c13 + c44) * n1 * n3, numbersMatrix[0][2]));
    BOOST_CHECK(eqWithEps(c66 * n1 * n1 + c11 * n2 * n2 + c44 * n3 * n3, numbersMatrix[1][1]));
    BOOST_CHECK(eqWithEps((c13 + c44) * n2 * n3, numbersMatrix[1][2]));
    BOOST_CHECK(eqWithEps(c44 * (n1 * n1 + n2 * n2) + c33 * n3 * n3, numbersMatrix[2][2]));

    int i, j;
    for(i = 0; i < 3; ++i)
        for(j = i + 1; j < 3; ++j)
            BOOST_CHECK(eqWithEps(numbersMatrix[i][j], numbersMatrix[j][i]));
    
}

BOOST_AUTO_TEST_CASE( test_det )
{
    PolyMatrix m1 = makeZeroPolyMatrix();
    m1[0][0] = Poly(1);
    m1[0][1] = Poly(2);
    m1[0][2] = Poly(4);
    m1[1][0] = Poly(1);
    m1[1][1] = Poly(3);
    m1[1][2] = Poly(9);
    m1[2][0] = Poly(1);
    m1[2][1] = Poly(4);
    m1[2][2] = Poly(16);
    Poly p1 = det(m1);

    BOOST_CHECK(polyEqWithEps(p1, Poly(2)));

    PolyMatrix m2 = makeZeroPolyMatrix();
    m2[0][0] = Poly(0, 1);
    m2[0][1] = Poly(1, 1);
    m2[0][2] = Poly(1, 2);
    m2[1][0] = Poly(1, 3);
    m2[1][1] = Poly(1, 5);
    m2[1][2] = Poly(1, 6);
    m2[2][0] = Poly(1, 8);
    m2[2][1] = Poly(1, 9);
    m2[2][2] = Poly(1, 11);
    Poly p2 = det(m2);

    BOOST_CHECK(eqWithEps(complex<double>(-78.), p2(2.)));
}


template<typename T>
ostream& operator<<(ostream &os, vector<T> const &dat)
{
    for(typename vector<T>::const_iterator it = dat.begin(); it != dat.end(); ++it)
        os << (*it) << " ";
    return os;
}

BOOST_AUTO_TEST_CASE( test_calcGammas )
{
    typedef complex<double> CD;

    CD g1 = CD(0.5, -1.);
    CD g2 = CD(-.6, -1.);
    CD g3 = CD(-0.9, -2.5);
    CD g4 = CD(g1.real(), -g1.imag());
    CD g5 = CD(g2.real(), -g2.imag());
    CD g6 = CD(g3.real(), -g3.imag());

    Poly p = Poly( -g1, 1) * Poly( -g4, 1) *
        Poly( -g2, 1) * Poly( -g5, 1) *
        Poly( -g3, 1) * Poly( -g6, 1);

    Poly::RootVec gammas = calcGammas(p);

    BOOST_CHECK_EQUAL( gammas.size(), 3);
    
    Poly::RootVec expected = list_of(g1) (g2) (g3);

    BOOST_CHECK(rootsEquality(expected, gammas));
}

BOOST_AUTO_TEST_CASE( test_kernel )
{
    ComplexMatrix m = makeZeroComplexMatrix();

    m[0][0] = 1;
    m[0][1] = 2;
    m[0][2] = 4;

    m[1][0] = 1;
    m[1][1] = 3;
    m[1][2] = 9;

    m[2][0] = 2;
    m[2][1] = 5;
    m[2][2] = 13;

    ComplexVec v = kernel(m);
    ComplexVec zero = list_of(0)(0)(0);
    ComplexVec actual = multiply(m, v);

    BOOST_CHECK(rootsEquality(zero, actual));
}


BOOST_AUTO_TEST_CASE( test_norm )
{
    typedef complex<double> CD;
    {
        ComplexVec v = list_of(3) (4) (0);
        ComplexVec actual = norm(v), expected = list_of(0.6)(0.8)(0);
        BOOST_CHECK(rootsEquality(expected, actual));
    }
    {
        ComplexVec v = list_of(CD(0,3)) (CD(0,4)) (CD(0));
        ComplexVec actual = norm(v), expected = list_of(CD(0,0.6))(CD(0,0.8))(0);
        BOOST_CHECK(rootsEquality(expected, actual));
    }
}

BOOST_AUTO_TEST_CASE( test_matrixOfRows )
{
    typedef complex<double> CD;
    typedef vector<CD> VCD;

    VCD r0 = list_of(1)(2)(3);
    VCD r1 = list_of(4)(5)(6);
    VCD r2 = list_of(7)(8)(9);

    ComplexMatrix m = matrixOfRows(r0, r1, r2);

    BOOST_CHECK(eqWithEps(m[0][0], r0[0]));
    BOOST_CHECK(eqWithEps(m[0][1], r0[1]));
    BOOST_CHECK(eqWithEps(m[0][2], r0[2]));
    BOOST_CHECK(eqWithEps(m[1][0], r1[0]));
    BOOST_CHECK(eqWithEps(m[1][1], r1[1]));
    BOOST_CHECK(eqWithEps(m[1][2], r1[2]));
    BOOST_CHECK(eqWithEps(m[2][0], r2[0]));
    BOOST_CHECK(eqWithEps(m[2][1], r2[1]));
    BOOST_CHECK(eqWithEps(m[2][2], r2[2]));
}

BOOST_AUTO_TEST_CASE( test_makeAmplitudeMatrix )
{
    // ComplexMatrix makeAmplitudeMatrix(MaterialTensor const& c,
    //                                   ComplexMatrix const& q, /* Matrix of kernels */
    //                                   ComplexMatrix const& n /* Matrix of directions */);

    MaterialTensor t = makeZeroTensor();
    ComplexMatrix q = makeZeroComplexMatrix();
    ComplexMatrix n = makeZeroComplexMatrix();

    t[1][2][2][0] = 3;
    t[1][2][0][1] = 2;
    q[2][0] = 5;
    q[0][2] = 11;
    n[2][1] = 7;
    n[0][0] = 13;

    ComplexMatrix actual = makeAmplitudeMatrix(t, q, n);

    int i, j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            if (i == 1) {
                if (j == 0) {
                    BOOST_CHECK(eqWithEps(429, actual[i][j]));
                    continue;
                } else if (j == 2) {
                    BOOST_CHECK(eqWithEps(70, actual[i][j]));
                    continue;
                }
            }
            BOOST_CHECK(eqWithEps(0, actual[i][j]));
        }
    }
}
