package bsl

type PolyMat [3][3]Poly

//Calculate determinant of polinomial matrix
//result is polinome
func PolyMatDet(pm PolyMat) Poly {
	t1 := Mult(pm[0][0], Subtr(Mult(pm[1][1], pm[2][2]), Mult(pm[2][1], pm[1][2])))
	t2 := Mult(pm[0][1], Subtr(Mult(pm[1][0], pm[2][2]), Mult(pm[2][0], pm[1][2])))
	t3 := Mult(pm[0][2], Subtr(Mult(pm[1][0], pm[2][1]), Mult(pm[2][0], pm[1][1])))
	
	return Add(Subtr(t1, t2), t3)
}

//generates polinomial matrix for given material tensor
//slowness sx, sy and density of material
func GenPolyMat(c Tensor4, sx, sy, rho complex128) (ret PolyMat) {
       ssel := [2]complex128{sx, sy}
        for p := 0; p < 3; p++ {
                for q := 0; q < 3; q++ {
                       ret[p][q] = Poly{0, 0, 0}
                       for i := 0; i < 3; i++ {
                               for k := 0; k < 3; k++ {
                                       term := c[i][p][k][q]
                                       d := 0
                                       if i == 2 {d++} else {term *= ssel[i]}
                                       if k == 2 {d++} else {term *= ssel[k]}

                                       ret[p][q][d] += term
                               }
                       }
                       if p == q {
                               ret[p][q][0] -= rho
                       }
                }
        }

        return
}

//evaluate each polinome in matrix to value and return 
//resulted matrix
func PolyMatEval(mat PolyMat, x complex128) (ret Mat3) {
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			ret[p][q],_,_ = mat[p][q].Val(x)
		}
	}
	return
}
