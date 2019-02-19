package bsl

import (
	"fmt"
	"github.com/runningwild/go-fftw"
	"image"
	"image/color"
	"image/png"
	"log"
	"math/cmplx"
	"os"
)

//load picture as a matrix by filename
func LoadMatrix(flnm string) (ret [][]complex128) {
	source, err := os.Open(flnm)
	if err != nil {
		log.Fatal(err)
	}
	defer source.Close()

	pic, err := png.Decode(source)
	if err != nil {
		log.Fatal(err)
	}

	rc := pic.Bounds()
	minx, miny, maxx, maxy := rc.Min.X, rc.Min.Y, rc.Max.X, rc.Max.Y

	nx, ny := (maxx - minx), (maxy - miny)

	ret = fftw.Alloc2d(ny, nx)

	for p := miny; p < maxy; p++ {
		for q := minx; q < maxx; q++ {
			ret[p][q] = complex(float64(color.GrayModel.Convert(pic.At(q, p)).(color.Gray).Y)/255.0, 0)
		}
	}

	return
}

//save complex matrix to file with given filename
func SaveMatrix(dat [][]complex128, flnm string, isCentered bool) {
	var valMin, valMax float64
	flag := true

	var shp, shq int


	ny := len(dat)
	nx := len(dat[0])

	if isCentered {
		shp = ny/2
		shq = nx/2
	}

	fmt.Println("save matrix with nx = ", nx, " ny = ", ny, " name = ", flnm)

	for p := 0; p < ny; p++ {
		for q := 0; q < nx; q++ {
			if flag {
				flag = false
				valMin = cmplx.Abs(dat[p][q])
				valMax = valMin
			} else {
				val := cmplx.Abs(dat[p][q])
				if val < valMin {
					valMin = val
				}
				if val > valMax {
					valMax = val
				}
			}
		}
	}

	pts := image.NewGray(image.Rect(0, 0, nx, ny)) 
	for p := 0; p < ny; p++ {
		for q := 0; q < nx; q++ {
			val := cmplx.Abs(dat[(p+shp)%ny][(q+shq)%nx])
			pts.Set(q, p, color.Gray{uint8(255*(val - valMin)/(valMax-valMin))})
		}
	} 

	dest, err := os.Create(flnm)
	if err != nil {
		log.Fatal(err)
	}
	defer dest.Close()

	png.Encode(dest, pts)
}

//module of matrix
func MatrixAbs(dat [][]complex128) (ret float64) {
	for p := 0; p < len(dat); p++ {
		for q := 0; q < len(dat[0]); q++ {
			ret += cmplx.Abs(dat[p][q]*dat[p][q])
		}
	}
	return 
}
