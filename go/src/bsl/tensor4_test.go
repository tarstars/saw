package bsl

import (
	"testing"
	//"fmt"
)

//test of set/get methods in Tensor4 
func TestTensor4SetGet(t *testing.T) {
	val := 5.0+0.0i
	mt := &Tensor4{}
	mt[1][1][2][2] = val
	if rv := mt[1][1][2][2]; val != rv {
		t.Error("Expected: ", val, " received ", rv)
	} 
}

//test of tetragonal material tensor creator
func TestMakeTetragonalMaterialTensor(t *testing.T) {
	c11 := 5.6e10 + 0.0i
	c12 := 5.145e10 + 0.0i
	c13 := 2.2e10 + 0.0i
	c33 := 10.6e10 + 0.0i
	c66 := 6.6e10 + 0.0i
	c44 := 2.65e10 + 0.0i

	mt := NewTetragonalMaterialTensor(c11, c12, c13, c33, c44, c66) 

	if rv := mt[0][0][0][0]; rv != c11 {
		t.Error("Expected: ", c11, " received ", rv)
	} 
	
	if rv := mt[0][0][1][1]; rv != c12 {
		t.Error("Expected: ", c12, " received ", rv)
	} 

	if rv := mt[2][2][2][2]; rv != c33 {
		t.Error("Expected: ", c33, " received ", rv)
	} 

	if rv := mt[1][2][2][1]; rv != c44 {
		t.Error("Expected: ", c44, " received ", rv)
	} 

	if rv := mt[1][0][0][1]; rv != c66 {
		t.Error("Expected: ", c66, " received ", rv)
	} 
}

//test of combInd function
func TestCombInd(t *testing.T) {
	val := 0
	if rv := combInd(0, 0); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 1
	if rv := combInd(1, 1); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 2
	if rv := combInd(2, 2); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 3
	if rv := combInd(1, 2); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 3
	if rv := combInd(2, 1); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 4
	if rv := combInd(2, 0); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 4
	if rv := combInd(0, 2); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 5
	if rv := combInd(0, 1); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 

	val = 5
	if rv := combInd(1, 0); rv != val {
		t.Error("Expected: ", val, " received ", rv)
	} 
}

func TestTrigonalPhotoelasticTensor(t *testing.T) {
	pt := NewTrigonalPhotoelasticTensor(11, 12, 13, 31, 33, 14, 41, 44)

	if pt[0][0][0][0] != 11 {
		t.Error("1111 problem")
	}

	if pt[0][0][1][1] != 12 {
		t.Error("1122 problem")
	}

	if pt[0][0][2][2] != 13 {
		t.Error("1133 problem")
	}
	
	if pt[0][0][1][2] != 14 {
		t.Error("1123 problem")
	}

	if pt[2][2][2][2] != 33 {
		t.Error("2222 problem")
	}

	if pt[1][1][1][2] != -14 {
		t.Error("1123 problem")
	}
}
