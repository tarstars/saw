package bsl

import ("fmt"
	"math")

//3x3 matrix
type Mat3 [3][3]complex128

//transpose matrix
func Transpose(m Mat3) (ret Mat3) {
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			ret[p][q] = m[q][p]
		}
	}
	return
}

//calculate determinant for complex matrix
func MatrixDet(m Mat3) complex128 {
	return  m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - 
		m[0][1] * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) + 
		m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1])
}

//calculate product of matrix and vector
func MatrixVectorProduct(m Mat3, v Vec3) (ret Vec3) {
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			ret[p] += m[p][q] * v[q]
		}
	}

	return
}

//pretty matrix print
func MatrixPPrint(m Mat3, pref string) {
	for p := 0; p < 3; p++ {
		fmt.Print(pref)
		for q := 0; q < 3; q ++ {
			val := real(m[p][q])
			if math.Abs(val) < 1e-6 {
				val = 0
			}
			fmt.Print(val, " ")
		}
		fmt.Println()
	}
}

//unity matrix
func NewUnityMatrix() (ret Mat3) {
	for p := 0; p < 3; p++ {
		ret[p][p] = 1
	}
	return
}
