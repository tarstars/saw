package bsl

import (
	"math"
	"math/cmplx"
)

//complex 3d vector
type Vec3 [3]complex128

//Vec3 from real numbers
func NewVec3(x, y, z float64) Vec3 {
	return Vec3{complex(x, 0), complex(y, 0), complex(z, 0)}
}

//scalar product
func DotProduct(x, y Vec3) complex128 {
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]
}

//cross product
func CrossProduct(x, y Vec3) (ret Vec3) {
	ret[0] =  x[1]*y[2] - x[2]*y[1]
	ret[1] = -x[0]*y[2] + x[2]*y[0]
	ret[2] =  x[0]*y[1] - x[1]*y[0]

	return
}

//absolute value of vector
func Abs(x Vec3) (s float64) {
	for _, v := range(x) {
		s += real(v * cmplx.Conj(v))
	}
	s = math.Sqrt(s)
	return
}

//return normalized vector
func Normalized(x Vec3) Vec3 {
	r := Abs(x)
	for p := 0; p < 3; p++ {
		x[p] /= complex(r, 0)
	}
	return x
}

//vector and complex product
func VectorComplexProduct(x Vec3, r complex128) (y Vec3) {
	for i, _ := range x {
		y[i] = x[i] * r
	}
	return
}

//vector and vector summ
func VectorSumm(x, y Vec3) (z Vec3) {
	for t := range x {
		z[t] = x[t] + y[t]
	}
	return
}

//vector and vector diff
func VectorDiff(x, y Vec3) (z Vec3) {
	for t := range x {
		z[t] = x[t] - y[t]
	}
	return
}

//equality of two vectors
func VecEq(x, y Vec3) bool {
	return Abs(VectorDiff(x, y)) < 1e-10
}

//rotate vec3 by rotation matrix
func Vec3Rotation(vec Vec3, m Mat3) (ret Vec3) {
	for k := 0; k < 3; k++ {
		for p := 0; p < 3; p++ {
			ret[k] += vec[p] * m[k][p]
		}
	}
	return 
}
