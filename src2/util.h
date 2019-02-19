/// @file util.h
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#ifndef __UTIL_H__INCLUDED__
#define __UTIL_H__INCLUDED__


#include "poly.h"

#include <ostream>

#include <boost/multi_array.hpp>


namespace farn {

    class Mat3;
    class MatrixFFTW;
    class PiezoTensor;
    class PlanFFTW;
    class SpacialMatrix;
    class Storage;
    class Vec3;
    class Vec3c;
    class VolumeData;
    class WaveMatrix;

    enum CalculationType {FixedForce, FixedDisplacement};

    struct TransmissPoint {
      double t, s1, s2;
      bool operator<(const TransmissPoint& r) const {return t < r.t;}
    };

    std::ostream& operator<<(std::ostream&, const TransmissPoint& r);

    typedef boost::multi_array<double, 4> MaterialTensor;
    typedef boost::multi_array<double, 2> RotationMatrix;

    MaterialTensor makeZeroTensor();
    MaterialTensor makeTetragonalMaterialTensor(double c11,
                                              double c12,
                                              double c13,
                                              double c33,
                                              double c44,
                                              double c66);

    MaterialTensor makeTrigonalMaterialTensor(double c11,
					      double c12,
					      double c13,
					      double c14,
					      double c33,
					      double c44, 
					      double c66);

    MaterialTensor makeParatelluriteMaterialTensor();
    MaterialTensor makeTelluriumMaterialTensor();

    PiezoTensor makeTrigonal3mMaterialTensor(double e15, 
					     double e21, 
					     double e31, 
					     double e33);


    RotationMatrix makePiezoChristoffelMatrix(MaterialTensor const& t, const PiezoTensor& pt, const Vec3& v, double eps11, double eps33);
					      
    MaterialTensor rotateTensor(MaterialTensor const& t, RotationMatrix const& m);

    typedef boost::multi_array<Poly, 2> PolyMatrix;
    PolyMatrix makePolyMatrix(MaterialTensor const&, double gamma = 0, double n1 = 1, double n2 = 0);
    PolyMatrix makeZeroPolyMatrix();
    
    RotationMatrix makeChristoffelMatrix(MaterialTensor const & r, Vec3 const & v);
    PolyMatrix makeChristoffelPolyMatrix(MaterialTensor const & r, Vec3 const& v);
    PolyMatrix makeChristoffelPolyMatrix(MaterialTensor const & t, double n1, double n2, double n3);

    typedef boost::multi_array<complex<double>, 2> ComplexMatrix;
    ComplexMatrix evaluatePolyMatrix(PolyMatrix const& m, complex<double> x);
    ComplexMatrix makeZeroComplexMatrix();

    Poly det(PolyMatrix const& m);
    complex<double> det(ComplexMatrix const& m);

    typedef vector<complex<double> > ComplexVec;
    ComplexVec kernel(ComplexMatrix const&);
    ComplexVec norm(ComplexVec const& v);
    ComplexVec multiply(ComplexMatrix const &, ComplexVec const&);
    double vecAbs(ComplexVec const& v);

    RotationMatrix multiply(RotationMatrix const& a, RotationMatrix const& b);

    Poly::RootVec calcGammas(Poly const&);
    ComplexMatrix makeAmplitudeMatrix(MaterialTensor const& c,
                                      ComplexMatrix const& q, /* Matrix of kernels */
                                      ComplexMatrix const& n /* Matrix of directions */);

    ComplexMatrix matrixOfRows(Poly::RootVec const& r0,
                               Poly::RootVec const& r1,
                               Poly::RootVec const& r2);

    RotationMatrix makeAxisRotation(double phi, int axis);
    RotationMatrix makeXYSlice(double phi);
    RotationMatrix minusGamma(const RotationMatrix& r, double gamma);

    complex<double> calcDet(MaterialTensor const &t, 
                            RotationMatrix const &m,
                            double gamma, 
                            double n1 = 1, 
                            double n2 = 0);

    RotationMatrix combineInMatrix(const Vec3&, const Vec3&, const Vec3&);
    RotationMatrix transpose(const RotationMatrix&);

    WaveMatrix create_wave_matrix(int n, double aperture, double freq, const MaterialTensor&, double rho, const Vec3c& force, CalculationType);
    WaveMatrix create_framed_wave_matrix(int n, double, double, double, double, double aperture, double freq, const MaterialTensor&, double rho, const Vec3c& force, CalculationType);

    RotationMatrix genChungMatrix(double alpha);
    RotationMatrix nearXMatrix(double alpha);

    Vec3c calcPol(const ComplexMatrix&);
    Mat3 strainFromKQ(const Vec3c& k, const Vec3c& q);
    Mat3 stressFromCS(const MaterialTensor& c, const Mat3& S);
    Mat3 columns2matrix(const Vec3c& v1, const Vec3c& v2, const Vec3c& v3);
    Poly makeCharacterPoly(const RotationMatrix& m);

    MatrixFFTW loadFromPicture(std::string flnm, double lower, double hi);
    void saveAsPicture(const MatrixFFTW&, std::string flnm);
    Storage layerTransform(const Storage&, MatrixFFTW&, MatrixFFTW&, const PlanFFTW&);

    bool rightPointing(double s1, 
		       double s2, 
		       double s3, 
		       const PolyMatrix& poly_mat, 
		       const MaterialTensor& tens, 
		       double omega);

    Vec3 operator*(const RotationMatrix&, const Vec3&);


    void saveAsPictures(const Storage& dat, std::string flnmBase);

    void copySliceIntoFftw(int, const Storage&, MatrixFFTW&);
    void copyFftwIntoSlice(int, const MatrixFFTW&, Storage&);
    SpacialMatrix getSpacialMatrixFromStorage(const Storage&);

    ostream& operator<<(ostream& os, const RotationMatrix& mat);
    ostream& operator<<(ostream& os, const PolyMatrix& m);

    ostream& getLog();

    Vec3 takeAbs(const Vec3c&);
    Vec3 takeReal(const Vec3c&);

    void drawVolumePovray(const std::string& flnm, const VolumeData& cs, double thr);
    VolumeData newSphericalVolumeData(int n);
    VolumeData newBarZVolumeData(int n, const Mat3&);

    template <typename T>
    std::ostream& operator<<(std::ostream& os, boost::multi_array<T, 4> const& tensor) {
        //TODO(dmitry): Make this generic for any dimensions vector.

        using namespace boost;
        using namespace std;

        int p, q, r, s;
        
        for (p = 0; p < 3; ++p) {
            for (q = 0; q < 3; ++q) {
                for (r = 0; r < 3; ++r) {
                    for (s = 0; s < 3; ++s) {
		      double val = tensor[p][q][r][s];
		      //if (abs(val) > 1) {
		      //os << p << q << r << s << " " << val << endl;
			//}
		      os << val << " ";
                    }
                    os << " | ";
                }
                os << endl;
            }
            os << endl;
        }
        return os;
    }

    template<typename T>
      ostream& operator<<(ostream &os, vector<T> const &dat)
      {
	for(typename vector<T>::const_iterator it = dat.begin(); it != dat.end(); ++it)
	  os << (*it) << " ";
	return os;
      }

    Vec3c operator*(const ComplexMatrix&, const Vec3c& );
    std::ostream& operator<<(std::ostream&, const ComplexMatrix&);
    double brezenDeviat(const Vec3& a, const Vec3& dir);
}

#endif // __UTIL_H__INCLUDED__
