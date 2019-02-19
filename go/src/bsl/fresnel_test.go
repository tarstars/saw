package bsl

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestFresnelIndex(t *testing.T) {
	phi := 13 * math.Pi / 180
	mo := Normalized(Vec3{0.0, complex(math.Sin(phi), 0), complex(math.Cos(phi), 0)})

	exx, ezz := 5.3973, 4.9890

	epss := Mat3{
		{complex(exx,0), 0, 0}, 
		{0, complex(exx, 0), 0}, 
		{0, 0, complex(ezz,0)}}

	indices :=IndexByDirection(mo, epss)

	if cmplx.Abs(indices[0].N - 2.31841326) > 1e-5 ||
		cmplx.Abs(indices[1].N - 2.323209) > 1e-5 {
		t.Error("indices", indices)	
	}
}
