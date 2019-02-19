package bsl

import (
	"testing"
)

//test function for determinant of polynomial matrix
func TestDet(t *testing.T) {
	a := PolyMat{
		{Poly{1}, Poly{1}, Poly{1}}, 
		{Poly{1}, Poly{2}, Poly{3}},
		{Poly{1}, Poly{4}, Poly{9}}}

	v1 := PolyMatDet(a)
	v2 := Poly{2}
	if !Eq(v1, v2) {
		t.Error("expected: ", v2, " received: ", v1)
	}
}

//test function for generation of poly matrix by sx, sy and material tensor
func TestGenPolyMat(t *testing.T) {
	mt := Tensor4{}
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			mt[2][p][2][q] = complex(float64((p + 1) * 10 + (q + 1)), 0)
		}
	}

	pm := GenPolyMat(mt, 1, 1, 1)

	if pm[0][0][0] != -1 {
		t.Error("must be -1")
	}
}
