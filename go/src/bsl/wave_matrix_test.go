package bsl

import (
	"math"
	"testing"
)

//test of creation of waveMatrix
func TestWaveMatrix(t *testing.T) {
	nx := 4 //discretes along x axis
	ny := 5 //discretes along y axis
	f := complex(1e8, 0)  //frequency 100 MHz
	a := complex(5e-3, 0) //aperture 5mm
	f0 := Vec3{1, 0, 0}
	mt, rho := NewParatelluriteMaterialTensor()
	mt = Tensor4Rotation(mt, NewChangMatrix(10 * math.Pi / 180))

	//wm := 
	NewWaveMatrix(nx, ny, mt, rho, a, f, f0)  
}
