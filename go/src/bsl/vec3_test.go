package bsl

import "testing"

//test initializaion of vec3
func TestVec3Init(t *testing.T) {
	x := complex(1.2, 0)
	y := complex(3.4, 0)
	z := complex(5.6, 0)

	a := Vec3{x, y, z}
	
	if rv := a[0]; rv != x {
		t.Error("Expected: ", x, " received ", rv)
	} 

	if rv := a[1]; rv != y {
		t.Error("Expected: ", y, " received ", rv)
	} 

	if rv := a[2]; rv != z {
		t.Error("Expected: ", z, " received ", rv)
	} 
}

//test dot product
func TestVec3DotProduct(t *testing.T) {
	a := Vec3{1, 2, 3}
	b := Vec3{10, 20, 30}

	result := DotProduct(a, b)
	expect := 140 + 0i

	if result != expect {
		t.Error("expected: ", expect, "result: ", result)
	}
}

//test cross product
func TestVec3CrossProduct(t *testing.T) {
	a := Vec3{1, 2,  4}
	b := Vec3{1, 5, 25}

	result := CrossProduct(a, b)
	expect := Vec3{30, -21, 3}
	if result != expect {
		t.Error("expected: ", expect, "result: ", result)
	}
}

//test vector-complex product
func TestVecComplexProduct(t *testing.T) {
	a := Vec3{1, 2, 4}
	b := Vec3{10, 20, 40}
	c := VectorComplexProduct(a, 10)

	if c != b {
		t.Error("expected: ", b, " real:", c)
	}
}
