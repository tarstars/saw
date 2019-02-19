package bsl

import (
	"fmt"
	"github.com/runningwild/go-fftw"
)


//type for storage 3D data
type Storage struct {
	h, d, w int 
	dat []complex128
	sdat [][]complex128
	ssdat [][][]complex128
}

//create storage by three dimensions
func NewStorage(h, d, w int) (ret *Storage) {
	fmt.Println("create storage h, d, w = ", h, d, w)

	ret = &Storage{h, d, w, make([]complex128, h*d*w), make([][]complex128, h*d), make([][][]complex128, h)}

	for p := 0; p < h*d; p++ {
		ret.sdat[p] = ret.dat[p*w:(p+1)*w]
	}

	for p := 0; p < h; p++ {
		ret.ssdat[p] = ret. sdat[p*d:(p+1)*d]
	}

	return 
}

//set element by indexes and value
func (pd *Storage) Set(r, p, q int, val complex128) {
	pd.dat[r * pd.w * pd.d + p * pd.w + q] = val
}

//get element by indexes
func (pd *Storage) At(r, p, q int) (complex128) {
	return pd.dat[r * pd.w * pd.d + p * pd.w + q]
} 

//layered Fourier transform for storage
func (pd *Storage) LayeredTransform(source, dest [][]complex128, plan *fftw.Plan) {
	for r := 0; r < pd.h; r++ {
		NonMoveCopy(source, pd.ssdat[r])
		plan.Execute()
		NonMoveCopy(pd.ssdat[r], dest)
	}

	ny, nx := len(dest), len(dest[0])
	for p := 0; p < ny; p++ {
		for q:=0; q < nx; q++ {
			dest[p][q] /= complex(float64(nx*ny), 0)
		}
	}
}

//get width
func (pd *Storage) Width() int {
	return pd.w
}

//get height
func (pd *Storage) Height() int {
	return pd.h
}

//get depth
func (pd *Storage) Depth() int {
	return pd.d
}

//get matrix by vertical coordinate
func (pd *Storage) GetMatrix(ind int) [][]complex128 {
	return pd.ssdat[ind]
}

//set values of the layer with given index by matrix
func (pd *Storage) SetZLayer(ind int, mat [][]complex128) {
	NonMoveCopy(pd.ssdat[ind], mat)
}

//get slice of storage with fixed depth
func (pd *Storage) SliceHW(d int) [][]complex128 {
	ret := fftw.Alloc2d(pd.h, pd.w)
	for p := 0; p < pd.h; p++ {
		for q := 0; q < pd.w; q++ {
			ret[p][q] = pd.At(p, d, q) 
		}
	}
	return ret
}

//get slice of storage with fixed width
func (pd *Storage) SliceHD(width int) [][]complex128 {
	ret := fftw.Alloc2d(pd.h, pd.d)
	for p := 0; p < pd.h; p++ {
		for q := 0; q < pd.d; q++ {
			ret[p][q] = pd.At(p, q, width) 
		}
	}
	return ret
}

//fill storage with 0
func (pd *Storage) Clear() {
	for t := 0; t < len(pd.dat); t++ {
		pd.dat[t] = 0
	}
}
