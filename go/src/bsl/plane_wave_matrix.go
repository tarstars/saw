package bsl

import (
	"github.com/runningwild/go-fftw"
)

type PlaneWaveMatrix struct {
	h, w int
	dat []PlaneWave
}

//set method for matrix
func (mat *PlaneWaveMatrix) Set(p, q int, w PlaneWave) {
	mat.dat[mat.w * p + q] = w
}

//receive element by index	
func (mat *PlaneWaveMatrix) At(p, q int) (*PlaneWave) {
	return &mat.dat[mat.w*p + q]
}

//construce plane wave matrix from storage
func NewPlaneWaveMatrix(st *Storage) (ret *PlaneWaveMatrix) {
	nx, ny := st.Width(), st.Depth()
	ret = &PlaneWaveMatrix{ny, nx, make([]PlaneWave, ny*nx)}

	for p := 0; p < ny; p++ {
		for q := 0; q < nx; q++ {
			ret.At(p,q).LoadFromStorage(p, q, st)
		}
	}

	return
}

//get spacial distribution of displacement module
func (pwm *PlaneWaveMatrix) GetDisplacementSlice() (ret [][]complex128) {
	ret = fftw.Alloc2d(pwm.h, pwm.w)
	for p := 0; p < pwm.h; p++ {
		for q := 0; q < pwm.w; q++ {
			ret[p][q] = complex(pwm.At(p, q).GetDisplacement(), 0)
		}
	}
	return
}

