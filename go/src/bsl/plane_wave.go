package bsl

import (
	"fmt"
	"io"
	"math"
	"math/cmplx"
)

type PlaneWave struct {
	slow        Vec3
	waveVector  Vec3
	strain      Mat3
	stress      Mat3
	energy      Vec3
	pol         Vec3
}

//construct plane wave in the frequency domain by the material tensor, slowness vector,
//frequency omega
func NewPlaneWave(mt Tensor4, slowness, polarization Vec3, omega complex128) (pv PlaneWave) {
	pv.slow = slowness
	pv.waveVector = VectorComplexProduct(slowness, omega)
	pv.strain = StrainByKQ(pv.waveVector, polarization)
	pv.stress = StressByStrainElasticity(pv.strain, mt)
	pv.energy = VectorComplexProduct(MatrixVectorProduct(pv.stress, polarization), complex(0,1))
	pv.pol    = polarization
	return
}

//get right column of tension tensor
func (pv PlaneWave) GetRightColumn() Vec3 {
	return Vec3{pv.stress[0][2], pv.stress[1][2], pv.stress[2][2]}
}

//check of physical sense
func (pv PlaneWave) HasPhysicalMeaning() bool {
	if math.Abs(imag(pv.slow[2])) > 1e-10 {
		return imag(pv.slow[2]) < 0
	}

	return real(pv.energy[2]) < 0
}

//print plane wave
func (pv PlaneWave) Print(dest io.Writer) {
	fmt.Fprintln(dest, "slow = ", pv.slow)
	fmt.Fprintln(dest, "waveVector = ", pv.waveVector)
	fmt.Fprintln(dest, "strain = ", pv.strain)
	fmt.Fprintln(dest, "stress = ", pv.stress)
	fmt.Fprintln(dest, "energy = ", pv.energy)
	fmt.Fprintln(dest, "polarization = ", pv.pol)
}

//test on homogenious wave
func (pv PlaneWave) IsInhomogenious() bool {
	return math.Abs(imag(pv.slow[2])) > 1e-6
}

//get slowness vector
func (pv PlaneWave) GetSlowness() Vec3 {
	return pv.slow
}

//get energy vector
func (pv PlaneWave) GetEnergyVector() Vec3 {
	return pv.energy
}

//get number of components to marshall
func (pv PlaneWave) GetNumberOfComponents() int {
	return 15
}

//fill storage with components
func (pv PlaneWave) FillStorage(p, q int, a complex128, st *Storage) {
	for t := 0; t < 3; t++ {
		st.Set(0 + t, p, q, st.At(0 + t, p, q) + a * pv.pol[t])
	}

	trp := []int{1, 2, 3, 2, 1, 1}
	trq := []int{1, 2, 3, 3, 3, 2}

	for t := 0; t < 6; t++ {
		addVal := a * pv.strain[ trp[t]-1 ][ trq[t]-1]
		st.Set(3 + t, p, q, st.At(3+t, p, q) + addVal)
	}

	for t := 0; t < 6; t++ {
		addVal := a * pv.stress[ trp[t]-1 ][ trq[t]-1]
		st.Set(9 + t, p, q, st.At(9+t, p, q) + addVal)
	}
}

//load data from storage. Inverse to FillStorage
func (pv *PlaneWave) LoadFromStorage(p, q int, st *Storage) {
	for t := 0; t < 3; t++ {
		pv.pol[t] = st.At(0, p, q)
	}

	trp := []int{1, 2, 3, 2, 1, 1}
	trq := []int{1, 2, 3, 3, 3, 2}

	for t := 0; t < 6; t++ {
		val := st.At(3+t, p, q)
		pv.strain[ trp[t]-1 ][ trq[t]-1] = val
		pv.strain[ trq[t]-1 ][ trp[t]-1] = val
	}

	for t := 0; t < 6; t++ {
		val := st.At(9+t, p, q)
		pv.stress[ trp[t]-1 ][ trq[t]-1] = val
		pv.stress[ trq[t]-1 ][ trp[t]-1] = val
	}

}

//get displacement
func (pv *PlaneWave) GetDisplacement() (ret float64) {
	return Abs(pv.pol)
}

//get phaze multiplyer by dz
func (pv *PlaneWave) GetPhazeMultiplier(dz complex128) (complex128) {
	return cmplx.Exp(-complex(0, 1)*dz * pv.waveVector[2])
}
