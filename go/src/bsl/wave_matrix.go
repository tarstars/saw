package bsl

import (
	"math"
)

//matrix of composite waves. each wave has own distribution of parameters
//on xy plane
type WaveMatrix struct {
	h, w int
	dat []CompositeWave
}

//slowness by position in array in fourier space
func SlownessByIndex(ind, limit int, aperture, freq complex128) complex128 {
	if ind > limit / 2 {
		ind -= limit
	}

	return complex(float64(ind), 0.0) * ( 1.0 / (freq*aperture) )
} 

//set method for matrix
func (mat *WaveMatrix) Set(p, q int, w CompositeWave) {
	mat.dat[mat.w * p + q] = w
}

//receive element by indexes	
func (mat *WaveMatrix) At(p, q int) (*CompositeWave) {
	return &mat.dat[mat.w*p + q]
}

//create wave matrix by number of discretes along x and y axes,
//material elasticity tensor, material density, frequency of external action
func NewWaveMatrix(h, w int, mt Tensor4, rho, ax, ay, freq complex128, f0 Vec3) (ret WaveMatrix) {
	ret.h = h
	ret.w = w
	ret.dat = make([]CompositeWave, h * w)

	for p := 0; p < h; p++ {
		for q := 0; q < w; q++ {
			sx := SlownessByIndex(q, w, ax, freq)
			sy := SlownessByIndex(p, h, ay, freq)

			ret.Set(p, q, NewCompositeWave(mt, rho, sx, sy, freq * 2.0 * math.Pi, f0))
		}
	}
	return
}

//take fourier amplitudes from a matrix and place them into composite waves in the wave matrix
func (mat *WaveMatrix) LoadProfileImage(pi [][]complex128) {
	//todo: panic if dimensions mismatches

	for p := 0; p < mat.h; p++ {
		for q := 0; q < mat.w; q++ {
			mat.At(p, q).SetAmplitude(pi[p][q])
		}
	}
}

//fill storage by data from the wave matrix
func (mat *WaveMatrix) FillStorage(st *Storage) {
	st.Clear()
	for p := 0; p < mat.h; p++ {
		for q := 0; q < mat.w; q++ {
			mat.At(p, q).FillStorage(p, q, st)
		}
	}
}

//takes into account shift along the z-axis
func (mat *WaveMatrix) ShiftZ(dz complex128) {
	for p := 0; p < mat.h; p++ {
		for q := 0; q < mat.w; q++ {
			mat.At(p, q).ShiftZ(dz)
		}
	}
}
