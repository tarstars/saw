package bsl

import (
	"io"
	"fmt"
	"log"
	//"os"
	"math"
)

//composite wave consists of three plane waves with equal projection
//of wave vector on XY plane 
type CompositeWave struct {
	pv [3]PlaneWave
	b  [3]complex128
	a     complex128
	farnMat Mat3
	farnDet complex128
}

//make composite wave from material tensor, material density, slowness projection
//frequency omega, force direction f0. Current problem: not accurate selection
//of appropriate plane waves. It is necessary to use pointing vector
func NewCompositeWave(mt Tensor4, rho, sx, sy, omega complex128, f0 Vec3) (cv CompositeWave) {
	//fmt.Println("NewCompositeWave: sx, sy = ", sx, ", ", sy)
	//generate matrix of polynomes, make matrix determinant zero
	pm := GenPolyMat(mt, sx, sy, rho)
	//fmt.Println("poly matrix = ", pm, "\n")
	pol := PolyMatDet(pm)
	//fmt.Println("\npolynome = ", pol, "\n")
	ar := pol.AllRoots()
	//fmt.Println("\nar = ", ar)
	
	//for each sz value calculate plane wave with all intrinsic values
	ind := 0
	for _, v := range ar {
		slowness := Vec3{sx, sy, v}
		polarization, _ := CalcPol(PolyMatEval(pm, v))
		cand := NewPlaneWave(mt, slowness, polarization, omega)
		//fmt.Println("cand = ")
		//cand.Print(os.Stdout)
		//fmt.Println()
		if cand.HasPhysicalMeaning() {
			if ind < 3 {
				cv.pv[ind] = cand
			}
			ind ++
		}
	}

	if ind != 3 {
		log.Panic("number of selected roots is not 3. it is ", ind)
	}

	//linear system equation solver
	for p := 0; p < 3; p++ {
		cv.farnMat[p] = cv.pv[p].GetRightColumn();
	}

	//fmt.Println("baseMatrix = ", baseMatrix)

	cv.farnDet = MatrixDet(cv.farnMat);
	///fmt.Println("baseDet = ", baseDet)

	for p := 0; p < 3; p++ {
		helperMatrix := cv.farnMat
		helperMatrix[p] = f0
		cv.b[p] = MatrixDet(helperMatrix) / cv.farnDet
	}

	/*
	for p := 0; p < 3; p++ {
		fmt.Println("b[", p, "] = ", cv.b[p])
	}
	fmt.Println()
        */

	return
} 


func (cv CompositeWave) testCalcForce() (tf Vec3) {
	for p := 0; p < 3; p++ {
		tf = VectorSumm(tf, VectorComplexProduct(cv.pv[p].GetRightColumn(), cv.b[p]))
	}
	return
}

func (cv CompositeWave) Print(dest io.Writer) {
	fmt.Fprintln(dest, "a = ", cv.a)
	for p := 0; p < 3; p++ {
		fmt.Fprintln(dest, "b = ", cv.b[p])
		cv.pv[p].Print(dest)
		fmt.Fprintln(dest)
	}
	fmt.Fprintln(dest, "farnell matrix = ", cv.farnMat)
	fmt.Fprintln(dest, "farnell det = ", cv.farnDet)
}

func (cv CompositeWave) IsGoodForTest() bool {
	for t := 0; t < 3; t++ {
		if cv.pv[t].IsInhomogenious() {
			return false
		}
	}
	return true
}

func (cv CompositeWave) GetTestAlpha() float64 {
	ind := 0
	val := real(cv.pv[0].GetSlowness()[2])
	for t := 1; t < 3; t++ {
		if real(cv.pv[t].GetSlowness()[2]) > val {
			val = real(cv.pv[t].GetSlowness()[2])
			ind = t
		}
	}

	mvv := Normalized(cv.pv[ind].GetSlowness())

	return math.Acos(real(mvv[2]))
}

func (cv CompositeWave) GetTestWalkoff(dest io.Writer) float64 {
	ind := 0
	val := real(cv.pv[0].GetSlowness()[2])
	for t := 1; t < 3; t++ {
		if real(cv.pv[t].GetSlowness()[2]) > val {
			val = real(cv.pv[t].GetSlowness()[2])
			ind = t
		}
	}

	s := cv.pv[ind].GetSlowness()
	fmt.Fprintln(dest, real(s[0]), real(s[2]))

	mvv := Normalized(s)
	vp := Normalized(cv.pv[ind].GetEnergyVector())

	return math.Acos(real(DotProduct(mvv, vp)))
}

func (cv *CompositeWave) SetAmplitude(amp complex128) {
	cv.a = amp
}

func (cv *CompositeWave) FillStorage(p, q int, st *Storage) {
	for t := 0; t < 3; t++ {
		cv.pv[t].FillStorage(p, q, cv.a*cv.b[t], st)
	}
}

func (cv *CompositeWave) ShiftZ(dz complex128) {
	for r := 0; r < 3; r++ {
		cv.b[r] = cv.b[r] * cv.pv[r].GetPhazeMultiplier(dz)
	}
}

func (cv *CompositeWave) GetFarnDet() complex128 {
	return cv.farnDet
}
