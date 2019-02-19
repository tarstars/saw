package bsl

import ("io"
	"fmt"
	"math")

//tensor of 4 rang
type Tensor4  [3][3][3][3]complex128

//combine two indices into one
func combInd(p, q int) int {
	return [...]int{
		0, 5, 4,
		5, 1, 3,
		4, 3, 2}[p*3 + q]
}

//Make material tensor for tetragonal crystals
func NewTetragonalMaterialTensor(c11, c12, c13, c33, c44, c66 complex128) Tensor4 {
	ret := Tensor4{}
	
	mat := [6][6]complex128{
		{c11, c12, c13,   0.0,   0.0,   0.0},
   		{c12, c11, c13,   0.0,   0.0,   0.0},
		{c13, c13, c33,   0.0,   0.0,   0.0},
 		{0.0, 0.0, 0.0,   c44,   0.0,   0.0},
		{0.0, 0.0, 0.0,   0.0,   c44,   0.0},
		{0.0, 0.0, 0.0,   0.0,   0.0,   c66}}

	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
				}
			}
		}
	}

	return ret
}

//Make material tensor for orthorombic crystals
func NewOrthorombicMaterialTensor(c11, c12, c13, c22, c23, c33, c44, c55, c66 complex128) Tensor4 {
	ret := Tensor4{}
	
	mat := [6][6]complex128{
		{c11, c12, c13,   0.0,   0.0,   0.0},
   		{c12, c22, c23,   0.0,   0.0,   0.0},
		{c13, c23, c33,   0.0,   0.0,   0.0},
 		{0.0, 0.0, 0.0,   c44,   0.0,   0.0},
		{0.0, 0.0, 0.0,   0.0,   c55,   0.0},
		{0.0, 0.0, 0.0,   0.0,   0.0,   c66}}

	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
				}
			}
		}
	}

	return ret
}

//Make material tensor for orthorombic crystals
func NewOrthorombicPhotoelasticTensor(p11, p12, p13,    p21, p22, p23,    p31, p32, p33, p44, p55, p66 complex128) Tensor4 {
	ret := Tensor4{}
	
	mat := [6][6]complex128{
		{p11, p12, p13,   0.0,   0.0,   0.0},
   		{p21, p22, p23,   0.0,   0.0,   0.0},
		{p31, p32, p33,   0.0,   0.0,   0.0},
 		{0.0, 0.0, 0.0,   p44,   0.0,   0.0},
		{0.0, 0.0, 0.0,   0.0,   p55,   0.0},
		{0.0, 0.0, 0.0,   0.0,   0.0,   p66}}

	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
				}
			}
		}
	}

	return ret
}

//Make trigonal photoelastic tensor
func NewTrigonalPhotoelasticTensor(p11, p12, p13, p31, p33, p14, p41, p44 complex128) Tensor4 {
	ret := Tensor4{}

	mat := [6][6]complex128{
		{p11, p12, p13,  p14,   0,   0},
		{p12, p11, p13, -p14,   0,   0},
		{p31, p31, p33,    0,   0,   0},
		{p41,-p41,   0,  p44,   0,   0},
		{  0,   0,   0,    0, p44, p41},
		{  0,   0,   0,    0, p14, (p11 - p12) / 2}}

	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
				}
			}
		}
	}

	return ret
}

//Make tetragonal photoelastic tensor
func NewTetragonalPhotoelasticTensor(p11, p12, p13, p31, p33, p44, p66 complex128) Tensor4 {
	ret := Tensor4{}

	mat := [6][6]complex128{
		{p11, p12, p13,    0,   0,   0},
		{p12, p11, p13,    0,   0,   0},
		{p31, p31, p33,    0,   0,   0},
		{  0,   0,   0,  p44,   0,   0},
		{  0,   0,   0,    0, p44,   0},
		{  0,   0,   0,    0,   0,  p66}}

	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
				}
			}
		}
	}

	return ret
}


//Make material tensor for tetragonal crystals
func NewTrigonalMaterialTensor(c11, c12, c13, c14, c33, c44, c66 complex128) Tensor4 {
	ret := Tensor4{}
	
	mat := [6][6]complex128{
		{c11,  c12, c13,  c14,   0,   0},
		{c12,  c11, c13, -c14,   0,   0}, 
		{c13,  c13, c33,    0,   0,   0},
		{c14, -c14,   0,  c44,   0,   0},
		{  0,    0,   0,    0, c44, c14},
		{  0,    0,   0,    0, c14, c66}};

	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
				}
			}
		}
	}

	return ret
}

//pretty print for tensor
func Tensor4PPrint(dest io.Writer, tens Tensor4, pref string) {
	for p:=0; p < 3; p++ {
		for q:=0; q < 3; q++ {
			for r:=0; r < 3; r++ {
				for s:=0; s < 3; s++ {
					val := real(tens[p][q][r][s])
					if math.Abs(val) < 1e-8 {
						val = 0
					}

					fmt.Fprintln(dest, pref, p + 1, q + 1, r + 1, s + 1, val)
				}
			}
		}
	}
}

//decompose abbreviated subscript into the full form
func IndDecompose(ind int) (k, l int) {
	switch ind {
	case 0: return 0, 0
	case 1: return 1, 1
	case 2: return 2, 2
	case 3: return 2, 1
	case 4: return 2, 0
	case 5: return 1, 0
	}
	return 0, 0
}

//print tensor4 as 6x6 matrix
func Tensor4SixPrint(dest io.Writer, tens Tensor4) {
	for p := 0; p < 6; p++ {
		for q := 0; q < 6; q++ {
			i, j := IndDecompose(p)
			k, l := IndDecompose(q)
			fmt.Fprint(dest, real(tens[i][j][k][l]), " ")
		}
		fmt.Fprintln(dest) 
	}
}
