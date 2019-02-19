package bsl

import (
	//"math"
	"math/cmplx"
	//"fmt"
)

//polinome with complex coefficients
type Poly []complex128

//degree of polinome
func (p Poly) Deg() int {
	return len(p) - 1 
}

//decrease degree of polinome by divide on (x - x0)
func (p Poly) EliminateRoot(x complex128) Poly {
	n := len(p)
	ret := make([]complex128, n - 1)
	
	var z complex128
  	for q := n - 2; q >= 0; q-- {
		z = z * x + p[q + 1]
		ret[q] = z
	}
	
	return ret
}

//two polinomes are equal if their degrees are equal and
//coefficients are near too
func Eq(x, y Poly) bool {
	if x.Deg() != y.Deg() {
		return false
	}

	n := len(x)
	dev := 0.0
	for p := 0; p < n; p++ {
		dev += cmplx.Abs(x[p] - y[p])
	}

	return dev < 1e-10
}

//calculate value of polinome, it's derivate of first and second
//orders in given point 
func (p Poly) Val(x complex128) (v, dv, ddv complex128) {
	n := len(p)
	for t := n - 1; t >= 0; t-- {
		ddv = ddv * x + dv
		dv  =  dv * x +  v
		v   =   v * x + p[t]
	}

	ddv *= 2.0

	return 
}

//find one root of polinome by laguerre method
func (p Poly) OneRoot() (root complex128) {
	//fmt.Println("OneRoot:")
	root = 0.00001
	n := complex(float64(p.Deg()), 0.0)
	for nit := 0; nit < 100; nit++ {
		//fmt.Println("\troot = ", root)
		v, dv, ddv := p.Val(root)
		//fmt.Println("\tv, dv, ddv = ", v, dv, ddv)
		if cmplx.Abs(v) <= 1e-60 {
			return
		}
		g := dv / v
		g2 := g * g
		h := g2 - ddv / v
		sq := cmplx.Sqrt((n - 1) * (n * h - g2))
		gp := g + sq
		gm := g - sq
		if (cmplx.Abs(gp) < cmplx.Abs(gm)) {
			gp = gm
		}
		dx := n / gp
		root -= dx
	}
	return
}

//find all polinomial roots 
func (p Poly) AllRoots() (roots []complex128) {
	tmp := p
	n := p.Deg()
	roots = make([]complex128, n)

	for t := 0; t < n; t++ {
		root := tmp.OneRoot()
		tmp = tmp.EliminateRoot(root)
		roots[t] = root
	}

	return 
}

//multiplication of two polinomes
func Mult(x, y Poly) (ret Poly) {
	ret = make([]complex128, x.Deg() + y.Deg() + 1)
	for p := 0; p < len(x); p++ {
		for q := 0; q < len(y); q++ {
			ret[p + q] += x[p] * y[q]
		}
	}

	return
}

//addition of two polinomes
func Add(x, y Poly) (ret Poly) {
	n := len(x)
	m := len(y)

	if n > m {
		n, m = m, n
		x, y = y, x
	}

	ret = make(Poly, m)
	for t, v := range y {
		ret[t] = v
		if t < n {
			ret[t] += x[t]
		}
	}

	return
}

//subtraction of two polinomes
func Subtr(x, y Poly) (ret Poly) {
	n := len(x)
	m := len(y)

	var s complex128
	s = -1

	if n > m {
		n, m = m, n
		x, y = y, x
		s = 1
	}

	ret = make(Poly, m)
	for t, v := range y {
		ret[t] = v
		if t < n {
			ret[t] -= x[t]
		}
		ret[t] *= s
	}

	return
}
		
//inplace multiplication of polinome on complex number
func MulCmplx(x *Poly, r complex128) {
	for t := range *x {
		(*x)[t] *= r
	}
}

//eliminate leadiing zeros from polinome
func (p *Poly) Reduce() {
	t := 0
	for r := range *p {
		if cmplx.Abs((*p)[r]) > 1e-10 {
			t = r
		}
	}

	*p = (*p)[:(t + 1)]
}
