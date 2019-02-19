package bsl

//tensor of 3 rang
type Tensor3  [3][3][3]complex128


//Make material matrix for paratellurite crystal
func NewTrigonalPiezoTensor(e15, e21, e31, e33 complex128) (ret Tensor3) {
	mat := [3][6]complex128{
		{   0,   0,   0,   0,  e15,  -e21},
		{-e21, e21,   0, e15,    0,     0},
		{ e31, e31, e33,   0,    0,     0}}

	for p := 0; p < 3; p ++ {
		for q := 0; q < 3; q ++ {
			for r := 0; r < 3; r++ {
				ret[p][q][r] = mat[p][combInd(q, r)]
			}
		}
	}
	
	return 
}

// //Make material tensor for tetragonal crystals
// func MakeTetragonalMaterialTensor(c11, c12, c13, c33, c44, c66 float64) Tensor4 {
// 	ret := Tensor4{}
	
// 	mat := [6][6]float64{
// 		{c11, c12, c13,   0.0,   0.0,   0.0},
//    		{c12, c11, c13,   0.0,   0.0,   0.0},
// 		{c13, c13, c33,   0.0,   0.0,   0.0},
//  		{0.0, 0.0, 0.0,   c44,   0.0,   0.0},
// 		{0.0, 0.0, 0.0,   0.0,   c44,   0.0},
// 		{0.0, 0.0, 0.0,   0.0,   0.0,   c66}}

// 	for p:=0; p < 3; p++ {
// 		for q:=0; q < 3; q++ {
// 			for r:=0; r < 3; r++ {
// 				for s:=0; s < 3; s++ {
// 					ret[p][q][r][s] = mat[combInd(p, q)][combInd(r, s)]
// 				}
// 			}
// 		}
// 	}

// 	return ret
// }
