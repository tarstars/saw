package bsl

import (
	"testing"
	//"fmt"
	"math/cmplx"
)

//test calculation of polynome degree function
func TestPolyDeg(t *testing.T) {
	p:=Poly{1, 2, 3}
	
	d:=2
	if v:=p.Deg(); v != d {
		t.Error("expected: ", d, " real: ", v)
	}
}

//test of elimination of one root from polynome by bezu
func TestPolyEliminateRoot(t *testing.T) {
	p := Poly{1, 3, 3, 1}

	p1 := p.EliminateRoot(-1)
	nd := 2
	if v := p1.Deg(); v!=nd {
		t.Error("expected: ", nd, " received: ", v)
	}

	np := Poly{1, 2, 1}
	if !Eq(np, p1) {
		t.Error("expected: ", np, " received: ", p1)
	}
}

//test of equation between two polynome
func TestPolyEq(t *testing.T) {
	if Eq(Poly{1, 2, 3, 4}, Poly{1, 2, 3, 5}) {
		t.Error("must be not equal")
	}

	if !Eq(Poly{1, 2, 3, 4}, Poly{1, 2, 3, 4}) {
		t.Error("must be equal")
	}
}

//test of evaluation of polynome, derivation of polynome and second degree of polynome
//in given point
func TestPolyVal(t *testing.T) {
	v, dv, ddv := (&Poly{1, 3, 3, 1}).Val(-1)

	if v != 0 {
		t.Error("v != 0")
	}

	if dv != 0 {
		t.Error("dv != 0")
	}

	if ddv != 0 {
		t.Error("ddv != 0")
	}

	if v, _, _ := (&Poly{1, 3, 3, 1}).Val(0); v != 1 {
		t.Error("must be ", 1, " indeed ", v)
	}

	if _, dv, _ := (&Poly{1, 3, 3, 1}).Val(0); dv != 3 {
		t.Error("must be ", 3, " indeed ", dv)
	}

	if _, _, ddv := (&Poly{1, 3, 3, 1}).Val(0); ddv != 6 {
		t.Error("must be ", 6, " indeed ", ddv)
	}
}

//test of one root calculation for polynome
func TestPolyOneRoot(t *testing.T) {
	p := &Poly{-10, 17, -8, 1}
	root := p.OneRoot()
	
	if v, _, _ := p.Val(root); v != 0 {
		t.Error("non zero value: ", v, " for root ", root)
	}	
}

//test for value, derivation and second derivation of polynome
func TestPolyValue(t *testing.T) {
	if v, dv, ddv := (&Poly{-10, 17, -8, 1}).Val(5); v != 0 {
		t.Error("expected 0, received ", v, " ", dv, " ", ddv)
	}
}

//test of polynomial coefficients order 
func TestPolyOrder(t *testing.T) {
	p := Poly{1, 2, 3, 4}
	if p[0] != 1 {
		t.Error("reverse sequence expected")
	}
}

//test of all roots of given polynome
func TestAllRoots( *testing.T) {
	p := Poly{1, 3, 3, 1}
	//fmt.Println("poly before: ", p)
	p.AllRoots()
	//fmt.Println("poly after: ", p)
}

//test of all roots of given polynome
func TestAllRoots1(tst *testing.T) {
	p := Poly{1, 0, 0, 1}
	roots := p.AllRoots()

	for _, t := range roots  {
		v, _, _ := p.Val(t)
		if (cmplx.Abs(v) > 1e-15) {
			tst.Error("to high poly value for root: v = ", v, " root = ", t)
		}
	}
}

//test of two polynome multiplication
func TestMult(t *testing.T) {
	x := Poly{1, 2, 1}
	y := Poly{1, 3, 3, 1}
	z := Poly{1, 5, 10, 10, 5, 1}

	if v := Mult(x, y); !Eq(v, z) {
		t.Error("expected: ", z, " obtained: ", v)
	}
}

//test of two polynome addition
func TestAdd1(t *testing.T) {
	v := Add(Poly{1,2},Poly{1,2,3})
	
	rec := v.Deg()
	if expect := 2; rec != expect {
		t.Error("expected: ", expect, " received: ", rec)
	}
	
	expect1 := Poly{2, 4, 3}
	if !Eq(expect1, v) {
		t.Error("expected: ", expect1, " received: ", v)
	}	
}

//test of two polynome subtraction
func TestSubst1(t *testing.T) {
	v := Subtr(Poly{1, 2, 3},Poly{1, 2})
	
	rec := v.Deg()
	if expect := 2; rec != expect {
		t.Error("expected: ", expect, " received: ", rec)
	}
	
	expect1 := Poly{0, 0, 3}
	if !Eq(expect1, v) {
		t.Error("expected: ", expect1, " received: ", v)
	}	
}

//test of two polynome subtraction
func TestSubst2(t *testing.T) {
	v := Subtr(Poly{1}, Poly{0})
	
	rec := v.Deg()
	if expect := 0; rec != expect {
		t.Error("expected: ", expect, " received: ", rec)
	}
	
	expect1 := Poly{1}
	if !Eq(expect1, v) {
		t.Error("expected: ", expect1, " received: ", v)
	}	
}

//test of two polynome subtraction
func TestSubst3(t *testing.T) {
	v := Subtr(Poly{0}, Poly{1})
	
	rec := v.Deg()
	if expect := 0; rec != expect {
		t.Error("expected: ", expect, " received: ", rec)
	}
	
	expect1 := Poly{-1}
	if !Eq(expect1, v) {
		t.Error("expected: ", expect1, " received: ", v)
	}	
}

//test of multiplication of polynome on complex number
func TestMultComplex(t *testing.T) {
	x := Poly{1,2,3,4}
	MulCmplx(&x, 3.0)
	expect := Poly{3, 6, 9, 12}

	if !Eq(x, expect) {
		t.Error("expected: ", expect, " received: ", x)
	}
}

//elimination of trailing zeros in polynome
func TestReduce(t *testing.T) {
	x := &Poly{1, 2, 1, 0, 0}

	x.Reduce()
	v := &Poly{1, 2, 1}

	if !Eq(*v, *x) {
		t.Error("expected: ", v, " received: ", x)
	}

	y := &Poly{0, 0, 0, 0}

	y.Reduce()
	if y.Deg() != 0 {
		t.Error("degree of {0, 0, 0} polynome must be zero")
	}

	if len(*y) != 1 {
		t.Error("one record must last in polynome after the reduce")
	}
	
}
