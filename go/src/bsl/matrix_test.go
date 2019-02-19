package bsl

import "testing"

//Test basic matrix operations
func TestMatrixSetGet(t *testing.T) {
	a := &Mat3{}

	if a[0][0] != 0 {
		t.Error("new give not zero matrix")
	}

	v := complex(3.14, 0)
	a[1][1] = v
	if rv := a[1][1]; rv != v {
		t.Error("expected: ", v, " received: ", rv)
	}
}

//Test calculation of matrix determinant
func TestMatrixDet(t *testing.T) {
	a := Mat3{{7, 13, 25}, {25, 17, 19}, {22, 20, 20}}

	v := MatrixDet(a)
	expected := 1804 + 0i

	if v != expected {
		t.Error("expected: ", expected, " actual: ", v)
	}
}

//Test of matrix and vector product
func TestMatrixVectorProduct(t *testing.T) {
	a := Mat3{{11, 12, 13}, {21, 22, 23}, {31, 32, 33}}

	v := Vec3{1, 0, 0}
	v1 := MatrixVectorProduct(a, v)
	expected := Vec3{11, 21, 31}
	if v1 != expected {
		t.Error("expected: ", expected, " received:", v1)
	}

	v = Vec3{0, 1, 0}
	v1 = MatrixVectorProduct(a, v)
	expected = Vec3{12, 22, 32}
	if v1 != expected {
		t.Error("expected: ", expected, " received:", v1)
	}

	v = Vec3{0, 0, 1}
	v1 = MatrixVectorProduct(a, v)
	expected = Vec3{13, 23, 33}
	if v1 != expected {
		t.Error("expected: ", expected, " received:", v1)
	}

}
