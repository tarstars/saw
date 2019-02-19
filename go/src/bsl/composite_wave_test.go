package bsl

import (
	"math"
	//"os"
	"testing"
)

func TestNewCompositeWave(t *testing.T) {
	mt, rho := NewParatelluriteMaterialTensor()
	mt = Tensor4Rotation(mt, NewChangMatrix(10 * math.Pi / 180))
	cv := NewCompositeWave(mt, rho, 0.000222, 0.000036, 2 * math.Pi * 1e8, Vec3{1,0,0})
	_ = cv
	//cv.Print(os.Stdout)
}

