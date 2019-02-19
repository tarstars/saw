package bsl

import (
	//"fmt"
	"math"
)

type IndexRet struct {
	N complex128
	QE, QD Vec3
}

//calculates index of refraction by direction of light propagation and
//dielectric permettivity tensor
func IndexByDirection(mo Vec3, eps Mat3) []IndexRet {
	var a PolyMat

	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			c1 := -mo[p] * mo[q]
			c0 := -eps[p][q]
			if p == q {
				c1 += 1
			}
			a[p][q] = Poly{c0, c1}
			(&a[p][q]).Reduce()
		}
	}

	charPoly := PolyMatDet(a)
	(&charPoly).Reduce()
	ar := charPoly.AllRoots()

	ret := make([]IndexRet, 0, 2)

	for _, gamm := range(ar) {
		if real(gamm) > 0 {
			zeroDetMat := PolyMatEval(a, gamm)
			qE, _ := CalcPol(zeroDetMat)
			qD := Normalized(MatrixVectorProduct(eps, qE))
			n := math.Sqrt(real(gamm))

			ret = append(ret, IndexRet{complex(n, 0), qE, qD})
		}
	}

	return ret
}

type IndexRetProj struct {
	N3 complex128
	N float64
	QE, QD Vec3
}

type Indexes []*IndexRetProj

func (s Indexes) Len() int {return len(s)}

func (s Indexes) Less(i, j int) bool {
	return real(s[i].N3) < real(s[j].N3)
}

func (s Indexes) Swap(i, j int) {s[i], s[j] = s[j], s[i]}

func IndexByProjection(n1, n2 complex128, eps Mat3) (ret Indexes) {
	polyMat := GenPolyMatrixForEpsilon(n1, n2, eps);
	//fmt.Println("polyMat = ", polyMat)
	pmd := PolyMatDet(polyMat)
	//fmt.Println("det = ", pmd)
	ar := pmd.AllRoots()
	for _, pdr := range(ar) {
		zm := PolyMatEval(polyMat, pdr)
		qE, _ := CalcPol(zm)
		qD := Normalized(MatrixVectorProduct(eps, qE))
		
		ret = append(ret, &IndexRetProj{pdr, 
			math.Sqrt(real(n1*n1 + n2*n2 + pdr*pdr)), 
			qE, 
			qD})
	}
	return 
}
