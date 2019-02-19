package bsl

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"
)

//test for christoffel matrix creation
func TestChristoffel(t *testing.T) {
	mt,_ := NewParatelluriteMaterialTensor() 
	cc := complex(1.0/math.Sqrt(2), 0)
	n := Vec3{cc, cc, 0}
	ch := Christoffel(n, mt)
	ch = ch
	fmt.Println("")
}

//test for calculation of polarization in real crystal
func TestPolarizationInCrystal(t *testing.T) {
	//fmt.Println("TestPolarization")
	mt,_ := NewParatelluriteMaterialTensor() 
	sx, sy := 0.0002 + 0i, 0.0003 + 0i
	pm := GenPolyMat(mt, sx, sy, 5.96e3)
	cp := PolyMatDet(pm)
	roots := cp.AllRoots()
	for _, v := range(roots) {
		cmat := PolyMatEval(pm, v)
		p, mat := CalcPol(cmat)
		val := Abs(MatrixVectorProduct(mat, p))
		if (val > 1e-15) {
			t.Error("p = ", p, val)
		}
	}
}

//test for CalcPol function for hypothetical matrices
func TestCalcPolarization(t *testing.T) {
	implTestOneMatrix(t, Mat3{
		{1, 1, 0}, 
		{1, 1, 0}, 
		{1, 0, 0}})

	implTestOneMatrix(t, Mat3{
		{1, 0, 0}, 
		{1, 1, 0}, 
		{1, 1, 0}})

	implTestOneMatrix(t, Mat3{
		{1, 1, 0}, 
		{1, 0, 0}, 
		{1, 1, 0}})

	implTestOneMatrix(t, Mat3{
		{(-3.814697265625e-05+0i), (0+0i), (0+0i)},
		{(0+0i), (3.5967239917487717e+09+0i), (3.5691281284236534e+10+0i)},
		{(0+0i), (3.569128128423653e+10+0i), (1.775372989816785e+11+0i)}})
}

//test of CalcPol for one matrix
func implTestOneMatrix(t *testing.T, mat Mat3) {
	p1, xmat := CalcPol(mat)
	m := Abs(MatrixVectorProduct(xmat, p1))

	if (!(math.Abs(1 - Abs(p1)) < 1e-14 && m < 1e-14)) {
		t.Error("Wrong polarization vector ", p1, " abs(mat * p1) = ", m)
	}
}

//test of CalcPol for real crystal in case of degeneracy
func TestCalcPolRealCrystalZeroDet(t *testing.T) {
	mt, _ := NewParatelluriteMaterialTensor()
	phi := 4.583662361046586
	n := Vec3{complex(math.Cos(phi), 0), complex(math.Sin(phi), 0), 0}
	_, polmat := NewChristoffelPoly(mt, n)
	mat := PolyMatEval(polmat, 2.65e10)
	q, nm := CalcPol(mat)
	if Abs(MatrixVectorProduct(nm, q)) > 1e-8 || math.Abs(1 - Abs(q)) > 1e-8 {
		t.Error("nm = ", nm, "\nq = ", q)
	}
}

//test for composite wave construction
func TestCompositeWave(t *testing.T) {
	mt, rho := NewParatelluriteMaterialTensor()
	f0 := Vec3{1, 0, 0}
	omega := complex(1e8, 0)
	cv := NewCompositeWave(mt, rho, 0.0001, 0.0002, omega, f0)
	ft := cv.testCalcForce()
	if !VecEq(ft, f0) {
		t.Error("test force = ", ft, " given force = ", f0)
	}
}

//test for algorithm of waves selection
func TestAlgorithmOfWavesSelection(t *testing.T) {
	h,w := 500, 500
	a := complex(5e-3, 0)
	freq := complex(1e8, 0)
	omega := 2 * math.Pi * freq
	mt, rho := NewParatelluriteMaterialTensor()
	mt = Tensor4Rotation(mt, NewChangMatrix(10 * math.Pi / 180))

	for p := 0; p < h; p++ {
		for q := 0; q < w; q++ {
			sx := SlownessByIndex(q, w, a, freq)
			sy := SlownessByIndex(p, h, a, freq)

			pm := GenPolyMat(mt, sx, sy, rho)
			pol := PolyMatDet(pm)
			ar := pol.AllRoots()

			num := 0
			realRoots := 0
			positiveRoots := 0

			for _, v := range(ar) {
				currMat := PolyMatEval(pm, v)
				q, _ := CalcPol(currMat)
				slow := Vec3{sx, sy, v}
				k := VectorComplexProduct(slow, omega)
				S := StrainByKQ(k, q)
				T := StressByStrainElasticity(S, mt)
				vgr := MatrixVectorProduct(T, q)
				if math.Abs(imag(v)) > 1e-10 {
					if imag(v) < -1e-10 {
						num ++;
					}
				} else {
					if imag(vgr[2]) < 0 {
						num ++;
					}
				}

				if math.Abs(imag(v)) < 1e-8 {
					realRoots ++
					if real(v) > 0 {
						positiveRoots ++
					}
				}
			}


			//if realRoots == 6 && positiveRoots != 3 {
				//fmt.Println("weird case of ", positiveRoots, " on ", realRoots - positiveRoots, " real roots sx = ", sx, " sy = ", sy)
			//}

			if num != 3 {
				t.Error("Number of selected roots is not 3. It is: ", num)
				return
			}
		}
	}

}

func TestMatrixByZ(t *testing.T) {
	vec := Normalized(Vec3{0.1, 0.2, 0.3})
	mat := MatrixByZ(vec)
	
	if Abs(VectorDiff(vec, mat[2])) > 1e-10 {
		t.Error("row 2 is not equal parameter")
	}

	if cmplx.Abs(DotProduct(vec, mat[0])) > 1e-10 {
		t.Error("row 2 is not equal parameter")
	}

	if cmplx.Abs(DotProduct(vec, mat[1])) > 1e-10 {
		t.Error("row 2 is not equal parameter")
	}
}

func TestCalcRealThetaPhi(t *testing.T) {
	theta, phi := CalcRealThetaPhi(NewVec3(1, 0, 0))
	if math.Abs(theta - math.Pi/2) > 1e-5 || math.Abs(phi - 0) > 1e-5 {
		t.Error("theta, phi = ", theta, phi)
	}

	for theta := 1.0; theta < 89.0; theta ++ {
		for phi := -179.0; phi < 179.0; phi ++ {
			theta0, phi0 := Rad(theta), Rad(phi)
			theta1, phi1 := CalcRealThetaPhi(NewVec3ByThetaPhi(theta0, phi0))
			if math.Abs(theta1 - theta0) > 1e-12 || math.Abs(phi1 - phi0) > 1e-12 {
				t.Error("theta, phi = ", theta, phi, Deg(theta1), Deg(phi1))
			}
		}
	}
}
