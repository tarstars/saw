//Package bsl implements library for acoustic calculations
//
//There are several types such as vectors, tensors and matrices.
//Polinomes with basic operations such as divide, multiplication
//summation, substraction and roots finding implemented as well.
package bsl

import (
	"bytes"
	"fmt"
	"flag"
	"image"
	"image/color"
	"image/png"
	"github.com/runningwild/go-fftw"
	"log"
	"math"
	"math/cmplx"
	"os"
	"runtime"
	"sort"
)

//test evanescent energy flow
func TestEvanescentEnergyFlow() {
	mt,rho := NewParatelluriteMaterialTensor()
	_ = NewCompositeWave(mt, rho, 0.01, 0.02, 1e8, Vec3{1, 0, 0})
}

//test laguerre method
func TestLaguerre() {
	p := Poly{1,2,3,4,5}
	fmt.Println("p = " , p)
	fmt.Println("all roots = ", p.AllRoots())
}

//test for slowness calculation
func TestSlownessCalc() {
	fmt.Println("test slowness")
	mt,rho := NewParatelluriteMaterialTensor()
	n := Normalized(NewVec3(1, 1, 0))
	chPoly, _ := NewChristoffelPoly(mt, n)
	dr := chPoly.AllRoots()
	fmt.Println("direction roots are", dr)

	sx, sy := 0.0005, 0.0003
	pm := GenPolyMat(mt, complex(sx, 0), complex(sy, 0), rho)
	pol := PolyMatDet(pm)
	ar := pol.AllRoots()

	fmt.Println("all roots are ", ar)
}

//Christoffel returns christoffel matrix for given direction n in material with
//tensor elasticity t
func Christoffel(n Vec3, t Tensor4) *Mat3 {
	ret := &Mat3{}

	for p:=0; p < 3; p ++ {
		for r:=0; r < 3; r++ {
			var v complex128
			for q:=0; q < 3; q++ {
				for s:=0; s < 3; s++ {
					v += t[p][q][r][s] * n[q] * n[s]
				}
			}
			ret[p][r] =v
		}
	}

	return ret
}

//function for creation of roots number color map
func setPoint(img *image.RGBA, mt Tensor4, sx, sy float64, p, q int) {
	/*pal := [...]color.RGBA{
		{  0,   0,   0, 255}, {255,   0,   0, 255}, {  0, 255,   0, 255},
		{  0,   0, 255, 255}, {255, 255,   0, 255}, {255,   0, 255, 255},
		{  0, 255, 255, 255}}*/

	c1, c2, c3, c4 := byte(255), byte(200), byte(140), byte(80)

	pal := [...]color.RGBA{
		{  c1,   c1,   c1, 255}, {c2,   c2,   c2, 255}, {  c3, c3,   c3, 255},
		{  c4,   c4, c4, 255}, {255, 255,   0, 255}, {255,   0, 255, 255},
		{  0, 255, 255, 255}}


	pm := GenPolyMat(mt, complex(sx, 0), complex(sy, 0), 5.96e3)
	pol := PolyMatDet(pm)
	ar := pol.AllRoots()

	v := 0
	for _, root := range ar {
		if (math.Abs(imag(root)) < 1e-10) {
			v++
		}
	}

	img.Set(q, p, pal[v/2])
}

//structure for parameters of setPoint call
type PointPoolTask struct {
	img *image.RGBA
	mt Tensor4
	sx, sy float64
	p, q int
}

//for the Runnable interface
func (pt *PointPoolTask) Run() {
	setPoint(pt.img, pt.mt, pt.sx, pt.sy, pt.p, pt.q)
}

//test indexes in rgba imaage
func TestImageDirections() {
	n := 500
	cx, cy := n/2, n/2
	w, h := 10, 100

	img := image.NewRGBA(image.Rect(0, 0, n, n))

	for p := 0; p < h; p++ {
		for q := 0; q < w; q++ {
			img.Set(cx - q, cy - p, color.RGBA{0, 255, 0, 255})
		}
	}

	dest, _ := os.Create("test_direction.png")
	png.Encode(dest, img)
}

//creation of color map of roots number
func CreateRootsMap() {
	runtime.GOMAXPROCS(4)

	mt_main_axis,_ := NewParatelluriteMaterialTensor()
	mt := Tensor4Rotation(mt_main_axis,  NewChangMatrix(10*math.Pi/180))

	ap := 0.0017
	n := 1200

	img := image.NewRGBA(image.Rect(0, 0, n, n))

	pool := NewPool(5)

	for p := 0; p < n; p++ {
		fmt.Println(p)
		for q := 0; q < n; q++ {
			sx := -ap + 2.0 * ap * float64(q) / float64(n)
			sy := -ap + 2.0 * ap * float64(p) / float64(n)
			
			pool.AddTask(&PointPoolTask{img, mt, sx, sy, p, q})
		}
	}
	pool.WaitAll()

	step := 0.0002
	for sx := -float64(int(ap / step)) * step; sx < ap; sx += step {
		p := n / 2
		q := int((sx + ap) / (2.0 * ap) * float64(n))
		img.Set(q, p, color.RGBA{255, 0, 0, 255})
	}

	for sy := -float64(int(ap / step)) * step; sy < ap; sy += step {
		q := n / 2
		p := int((sy + ap) / (2.0 * ap) * float64(n))
		img.Set(q, p, color.RGBA{255, 0, 0, 255})
	}

	dest, _ := os.Create("a.png")
	png.Encode(dest, img)
}

//set sphere into the povray scene for all roots including imaginary
func setSphereWithImaginary(ch chan []byte, mt Tensor4, rho, sx, sy float64) {
	pm := GenPolyMat(mt, complex(sx, 0), complex(sy, 0), complex(rho, 0))
	pol := PolyMatDet(pm)
	ar := pol.AllRoots()

	for _, v := range ar {
		x, y, z := sx, sy, real(v)
		x, y, z = x * 1000, y * 1000, z * 1000

		var lb bytes.Buffer

		if (math.Abs(imag(v)) > 1e-10) {
			fmt.Fprintf(&lb, "#declare customColor = <1, 0, 0>;")
		} else {
			fmt.Fprintf(&lb, "#declare customColor = <0, 1, 0>;")
		}

		if math.Abs(z) > 1e-5 {
			fmt.Fprintln(&lb, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{ color customColor }}")
			ch <- lb.Bytes()
		}

		lb = bytes.Buffer{}
		z = imag(v) * 1000
		if math.Abs(z) > 1e-5 {
			fmt.Fprintln(&lb, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{ color <0, 0, 1> }}")
			ch <- lb.Bytes()
		}			

	}
}

//set sphere into the povray scene for all roots including imaginary
func setSphereReal(ch chan []byte, mt Tensor4, rho, sx, sy complex128) {
	pm := GenPolyMat(mt, sx, sy, rho)
	pol := PolyMatDet(pm)
	ar := pol.AllRoots()

	for _, v := range ar {
		x, y, z := sx, sy, real(v)
		x, y, z = x * 1000, y * 1000, z * 1000

		var lb bytes.Buffer

		if math.Abs(imag(v)) < 1e-10 {
			fmt.Fprintln(&lb, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{ color <0,1,0> }}")
			ch <- lb.Bytes()
		}
	}
}

func setSphere(ch chan []byte, mt Tensor4, rho, sx, sy complex128) {
	setSphereReal(ch, mt, rho, sx, sy)
}

//structure for setSphere parameters
type SpherePoolTask struct {
	ch chan []byte
	mt Tensor4
	rho complex128
	sx, sy complex128
}

//Runnable interface implementation
func (spt *SpherePoolTask) Run() {
	setSphere(spt.ch, spt.mt, spt.rho, spt.sx, spt.sy)
}

//Create povray scene by calculation of six-degree equation roots
func CreatePovrayScene() {
	runtime.GOMAXPROCS(4)

	mt, rho := NewParatelluriteMaterialTensor()
	mt = Tensor4Rotation(mt, NewChangMatrix(10 * math.Pi / 180))

	ap := 0.0017
	n := 500

	dst,_ := os.Create("scene.pov")
	ch := make(chan []byte)
	go func(){
		var chunk []byte
		for {
			chunk = <- ch
			dst.Write(chunk)
		}
	}()

	pool := NewPool(5)

	for p := 0; p < n; p++ {
		fmt.Println(p)
		for q := 0; q < n; q++ {
			sx := complex(-ap + 2.0 * ap * float64(p) / float64(n), 0)
			sy := complex(-ap + 2.0 * ap * float64(q) / float64(n), 0)
			
			pool.AddTask(&SpherePoolTask{ch, mt, rho, sx, sy})
		}
	}
	pool.WaitAll()
}

//calculate polarization by the given matrix with
//zero determinant
func CalcPol(mat Mat3) (pol Vec3, nm Mat3) {
	n := 0
	S := 0.0
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			n ++
			S += cmplx.Abs(mat[p][q])
		}
	}

	S /= float64(n)

	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			nm[p][q] = mat[p][q] / complex(S, 0)
		}
	}

	return calcPolImpl(nm), nm 
}

//calculation of polarization for matrix with normalized rows
func calcPolImpl(mat Mat3) (ret Vec3) {
	//fmt.Println("mat = ", mat)

	//fmt.Println("try to find p")
	for p:=0; p < 3; p++ {
		//fmt.Println("\tp = ", p)
		i1 := (p + 1) % 3
		i2 := (p + 2) % 3
		//fmt.Println("\t p i1 i2 = ", Abs(mat[p]),  Abs(mat[i1]), Abs(mat[i2]))

		if Abs(mat[p]) > 1e-3 && Abs(mat[i1]) < 1e-3 && Abs(mat[i2]) < 1e-3 {
			c1 := Vec3{-mat[p][1],  mat[p][0], 0}
			c2 := Vec3{         0, -mat[p][2], mat[p][1]}

			if Abs(c1) > 1e-8 {
				return Normalized(c1)
			}

			if Abs(c2) > 1e-8 {
				return Normalized(c2)
			}
		}
	}

	//fmt.Println("\t can't find p")
	//fmt.Println("calcPolImpl for matrix\n", mat, "\n")
	for p := 0; p < 3; p++  {
		q := (p + 1) % 3
		r := (p + 2) % 3
		a := Abs(mat[p])
		if math.Abs(a) < 1e-8 {
			//fmt.Println("answer found by rows ", q, " and ", r)
			return Normalized(CrossProduct(mat[q], mat[r]))
		}
	}

	//fmt.Println("before row normalization:\n", mat, "\n")

	for p := 0; p < 3; p++ {
		mat[p] = Normalized(mat[p])
	}

 	//fmt.Println("after row normalization:\n", mat, "\n")
	
	ret = CrossProduct(mat[0], mat[1])
	val := Abs(ret)
	for p := 0; p < 2; p++ {
		cand := CrossProduct(mat[p + 1], mat[(p + 2) % 3])
		cand_val := Abs(cand)
		if (cand_val > val) {
			ret = cand
			val = cand_val
		}
	}

	ret = Normalized(ret)
	//fmt.Println("pol = ", ret)
	//fmt.Println("z = ", MatrixVectorProduct(mat, ret))
	
	return 
}

//strain tensor from wavevector and polarization
func StrainByKQ(k, q Vec3) (strain Mat3) {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			strain[i][j] = 0.5 * (k[i] * q[j] + k[j] * q[i]) * complex(0, 1)
		}
	}

	return
}

//stress tensor by strain tensor and elasticity tensor
func StressByStrainElasticity(strain Mat3, mt Tensor4) (stress Mat3) {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					stress[i][j] += mt[i][j][k][l] * strain[k][l]
				}
			}
		}
	}

	return
}

//christoffel matrix with piezoeffect
func NewPiezoChristoffel(pt Tensor3, mt Tensor4, n Vec3, tensEps Mat3) (ret Mat3) {
	eps := complex(0, 0)

	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			eps += tensEps[p][q] * n[p] * n[q]
		}
	}

	for p := 0; p < 3; p ++ {
		for q := 0; q < 3; q ++ {
			var val, gamma_p, gamma_q complex128

			for r := 0; r < 3; r++ {
				for s := 0; s < 3; s++ {
					val += mt[p][r][q][s] * n[r] * n[s]
					gamma_p += pt[r][p][s] * n[r] * n[s]
					gamma_q += pt[r][q][s] * n[r] * n[s]
				}
			}
			
			ret[p][q] = val + gamma_p * gamma_q / eps
		}
	}

	return
}

//christoffel matrix
func NewChristoffel(mt Tensor4, n Vec3) (ret Mat3) {
	for p := 0; p < 3; p ++ {
		for q := 0; q < 3; q ++ {
			var val complex128

			for r := 0; r < 3; r++ {
				for s := 0; s < 3; s++ {
					val += mt[p][r][q][s] * n[r] * n[s]
				}
			}
			
			ret[p][q] = val
		}
	}

	return
}


//rotation of tensor4 by given mat3
func Tensor4Rotation(mt Tensor4, mat Mat3) (ret Tensor4) {
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			for r := 0; r < 3; r++ {
				for s := 0; s < 3; s++ {
					for i := 0; i < 3; i++ {
						for j := 0; j < 3; j++ {
							for k := 0; k < 3; k++ {
								for l := 0; l < 3; l++ {
									ret[p][q][r][s] += mt[i][j][k][l] *
										mat[p][i] *
										mat[q][j] *
										mat[r][k] *
										mat[s][l] 
								}
							}
						}
					}
				}
			}
		}
	}
	return
}

//rotation of tensor3 by given mat3
func PiezoTensorRotation(mt Tensor3, mat Mat3) (ret Tensor3) {
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			for r := 0; r < 3; r++ {
				for i := 0; i < 3; i++ {
					for j := 0; j < 3; j++ {
						for k := 0; k < 3; k++ {
							ret[p][q][r] += mt[i][j][k] *
								mat[p][i] *
								mat[q][j] *
								mat[r][k] 
						}
					}
				}
			}
		}
	}
	return
}

//rotation of mat3 by given mat3
func Mat3Rotation(mat, rot Mat3) (ret Mat3) {
	for k := 0; k < 3; k++ {
		for m := 0; m < 3; m++ {
			for i:=0; i < 3; i++ {
				for j:=0; j < 3; j++ {
					ret[k][m] += mat[i][j] * 
						rot[k][i] * 
						rot[m][j]
				}
			}
		}
	}
	return 
}

//generate matrix for chang filter with alpha cut
// ex = {-1/sq2,  1/sq2, 0}
// ey = {-s/sq2, -s/sq2, c}
// ez = { c/sq2,  c/sq2, s}
func NewChangMatrix(alpha float64) Mat3 {
	c := complex(math.Cos(alpha), 0)
	s := complex(math.Sin(alpha), 0)
	sq2 := complex(math.Sqrt(2), 0)
	return Mat3{
		{-1.0/sq2,  1.0/sq2, 0},
		{  -s/sq2,   -s/sq2, c},
		{   c/sq2,    c/sq2, s}}
}

//generate matrix for x rotation
// ex = {     1,   0,  0}
// ey = {     0,   c,  s}
// ez = {     0,  -s,  c}
func NewRotxMatrix(alpha float64) Mat3 {
	c := complex(math.Cos(alpha), 0)
	s := complex(math.Sin(alpha), 0)
	return Mat3{
		{ complex(1.0, 0),  complex(0, 0), complex(0, 0)},
		{ complex(0, 0)  ,              c,             s},
		{ complex(0, 0)  ,             -s,             c}}
}

//build data for plotting the slownes surface section
func SlownessSection(flnm string) {
	mt, rho := NewParatelluriteMaterialTensor()
	mt = Tensor4Rotation(mt, NewChangMatrix(10 * math.Pi / 180))
	//a := Normalized(Vec3{0.000222, 0.000036, 0})
	a := Normalized(Vec3{0, 1, 0})
	b := Normalized(Vec3{0, 0, 1})

	destFile, _ := os.Create(flnm)
	defer destFile.Close()

	for phi := 0.0; phi < 2 * math.Pi; phi += 0.01 {
		n := VectorSumm(
			VectorComplexProduct(a, complex(math.Cos(phi), 0)),
			VectorComplexProduct(b, complex(math.Sin(phi), 0)))
		gp, _ := NewChristoffelPoly(mt, n)
		ar := gp.AllRoots()
		for _, r := range ar {
			s := cmplx.Sqrt(rho / r)
			rn := VectorComplexProduct(n, s)
			fmt.Fprintln(destFile, real(DotProduct(a, rn)), real(DotProduct(b, rn))) 
		}
	}
}

//Construct christoffel polynome for given direction
func NewChristoffelPoly(mt Tensor4, n Vec3) (Poly, PolyMat) {
	var dat PolyMat
	for p := 0; p < 3; p ++ {
		for q := 0; q < 3; q++ {
			fv := 0.0 + 0.0i
			for r := 0; r < 3; r++ {
				for s := 0; s < 3; s++ {
					fv += mt[p][r][q][s] * n[r] * n[s]
				}
			}

			if p == q {
				dat[p][q] = Poly{fv, -1}
			} else {
				dat[p][q] = Poly{fv}
			}
		}
	}
	return PolyMatDet(dat), dat
}

func GenPolyMatrixForEpsilon(n1, n2 complex128, ep Mat3) (ret PolyMat)  {
	ret[0][0] = Poly{n2*n2 - ep[0][0], 0, 1}
	ret[0][1] = Poly{-n1 * n2 - ep[0][1]}
	ret[0][2] = Poly{-ep[0][2], -n1 }


	ret[1][0] = Poly{-n1*n2 - ep[1][0]}
	ret[1][1] = Poly{ n1*n1 - ep[1][1], 0, 1}
	ret[1][2] = Poly{- ep[1][2], -n2}

	ret[2][0] = Poly{- ep[2][0], -n1}
	ret[2][1] = Poly{- ep[2][1], -n2}
	ret[2][2] = Poly{n1*n1+n2*n2 - ep[2][2]}

	return
}

//print velocities for given direction
func PrintVelocitiesInDirection() {
	mt, rho := NewParatelluriteMaterialTensor()
	cp, _ := NewChristoffelPoly(mt, Normalized(Vec3{1, 1, 0}))
	ar := cp.AllRoots()
	for _, r := range ar {
		fmt.Println("v = ", math.Sqrt(real(r) / real(rho)))
	}
}

//show points for 4 on 2 slice
func PrintSixPoints() {
	dst, _ := os.Create("six_points.txt")
	defer dst.Close()

	x, y := 0.000222, 0.000036
	r := math.Sqrt(x*x + y*y)
	mt, rho := NewParatelluriteMaterialTensor()
	mt = Tensor4Rotation(mt, NewChangMatrix(10 * math.Pi / 180))
	ar := PolyMatDet(GenPolyMat(mt, complex(x, 0), complex(y, 0), rho)).AllRoots()
	for _, v := range(ar) {
		fmt.Fprintln(dst, r, real(v))
	}	
}

//test group velocity
func GroupVelocityTest() {
	omega := 1e8
	mt, rho := NewParatelluriteMaterialTensor()
	dest, _ := os.Create("xy_group.txt")
	defer dest.Close()

	mt = Tensor4Rotation(mt, NewChangMatrix(10.0 * math.Pi / 180))
	
	for phi := 0.0; phi < math.Pi/2; phi += 0.001 {
		n := Vec3{complex(math.Cos(phi), 0), complex(math.Sin(phi), 0), 0}
		cp, mat := NewChristoffelPoly(mt, n)
		ar := cp.AllRoots()
		for _, gamm := range(ar) {
			matVal := PolyMatEval(mat, gamm)
			q, _ := CalcPol(matVal)
			v := math.Sqrt(real(gamm) / real(rho))
			k := VectorComplexProduct(n, complex(omega / v, 0))
			S := StrainByKQ(k, q)
			T := StressByStrainElasticity(S, mt)
			vgr := MatrixVectorProduct(T, q)
			nvgr := Normalized(vgr)
			dpv := DotProduct(n, nvgr)
			cpv := CrossProduct(n, nvgr)[2]
			dp, cp := imag(dpv), imag(cpv)
			psi := math.Atan2(cp, dp)
			fmt.Fprintln(dest, phi / math.Pi * 180, psi / math.Pi * 180)
		}
	}	
}

//test dependence of walkoff angle on cut angle
func TestWalkoffDependence() {
	mt, rho := NewParatelluriteMaterialTensor()
	dest, _ := os.Create("walkoff.txt")
	destSection, _ := os.Create("section.txt")
	defer dest.Close()

	mt = Tensor4Rotation(mt, NewChangMatrix(0.0 * math.Pi / 180))
	
	n := 10000
	A := 0.0017

	for p := 0 ; p < n; p++ {
		x := A / float64(n) * float64(p)

		cv := NewCompositeWave(mt, rho, complex(x, 0), 0, 1e8, NewVec3(1, 0, 0))
			
		if cv.IsGoodForTest() {
			fmt.Fprintln(dest, Deg(cv.GetTestAlpha()), Deg(cv.GetTestWalkoff(destSection)))
		}
	}
}

//build transition matrix by new z value
func MatrixByZ(m Vec3) (ret Mat3) {
	//m3 := m
	//m1 := Normalized(Vec3{0, -m[2], m[1]})
	//m2 := CrossProduct(m3, m1)
	//return Mat3{m1, m2, m3}
	return Mat3{
		Normalized(Vec3{m[1], -m[0], 0}), 
		Normalized(Vec3{m[0]*m[2], m[1]*m[2], -(m[0]*m[0] + m[1]*m[1])}),			
		m}
}

type VelocityResult struct {
	Vel float64
	Q Vec3
}

func CalcPiezoVelocities(m_ac Vec3, rot Mat3) (ret []VelocityResult, rho complex128) {
	ptMainAxis := NewLinbo3PiezoTensor()
	mtMainAxis, rho := NewLinbo3MaterialTensor()
	epsMainAxis := Linbo3StaticEps()

	pt := PiezoTensorRotation(ptMainAxis, rot)
	mt := Tensor4Rotation(mtMainAxis, rot)
	staticEps := Mat3Rotation(epsMainAxis, rot)

	christ := NewPiezoChristoffel(pt, mt, m_ac, staticEps)

	var pm PolyMat
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			if p == q {
				pm[p][q] = Poly{christ[p][q], -1}
			} else {
				pm[p][q] = Poly{christ[p][q]}
			}
		}
	}

	pmd := PolyMatDet(pm)
	gammas := pmd.AllRoots()

	for _, gam := range(gammas) {
		mat := PolyMatEval(pm, gam)
		pol, _:= CalcPol(mat)
		ret = append(ret, VelocityResult{math.Sqrt(real(gam)/real(rho)), pol})
	}

	return
}

func CalcVelocities(m_ac Vec3, rot Mat3, mtMainAxis Tensor4, rho complex128) (ret []VelocityResult) {
	mt := Tensor4Rotation(mtMainAxis, rot)

	christ := NewChristoffel(mt, m_ac)

	var pm PolyMat
	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			if p == q {
				pm[p][q] = Poly{christ[p][q], -1}
			} else {
				pm[p][q] = Poly{christ[p][q]}
			}
		}
	}

	pmd := PolyMatDet(pm)
	gammas := pmd.AllRoots()

	for _, gam := range(gammas) {
		mat := PolyMatEval(pm, gam)
		pol, _:= CalcPol(mat)
		ret = append(ret, VelocityResult{math.Sqrt(real(gam)/real(rho)), pol})
	}

	return
}


//get theta and phi in spherical coordinate system by 3D vector
func CalcRealThetaPhi(xv Vec3, revMat Mat3) (theta, phi float64) {
	v := Normalized(Vec3Rotation(xv, revMat))
	theta = math.Acos(real(v[2]))
	phi = math.Atan2(real(v[1]), real(v[0]))
	//fmt.Println("\t\t\t\tcalcRealThetaPhi(", xv, ") = ", " phi =", phi, " theta =", theta)
	return
}

func NewVec3ByThetaPhi(theta, phi float64) Vec3 {
	return NewVec3(math.Sin(theta)*math.Cos(phi), math.Sin(theta)*math.Sin(phi), math.Cos(theta))
}

func DispatchCalcFigureOfMerit() {
	flagValue := flag.String("crystal", "teo2", "you can select a crystal type using this flag")
	flag.Parse()

	lambda := 0.633e-6+0i

	fmt.Println("flag = ", *flagValue)
	switch *flagValue {
	case "teo2":
		fmt.Println("calculate teo2")
		tens, rho := NewParatelluriteMaterialTensor()
		CalcFigureOfMerit(Teo2eps(), NewTeo2Photoelastic(), tens, rho, NewVec3ByThetaPhi(Rad(80.0), Rad(45.0)), lambda)

	case "linbo3":
		fmt.Println("calculate linbo3")
		tens, rho := NewLinbo3MaterialTensor()
		CalcFigureOfMerit(Linbo3eps(), NewLinbo3Photoelastic(), tens, rho, NewVec3ByThetaPhi(/*Rad(23.3), Rad(-90)*/Rad(13.0), Rad(90.0))/**/, lambda)

	case "te":
		fmt.Println("calculate te")
		tens, rho := NewTeMaterialTensor()
		CalcFigureOfMerit(NewTeEps(), NewTePhotoelastic(), tens, rho, NewVec3ByThetaPhi(Rad(90.0), Rad(0)), lambda)

	case "sr":
		fmt.Println("calculate sr")
		tens, rho := NewSrTetraborateMaterialTensor()
		CalcFigureOfMerit(NewSrTetraborateEps(), NewSrTetraboratePhotoelasticTensor(), tens, rho, NewVec3ByThetaPhi(Rad(89.9), Rad(0)), lambda)

	case "pb":
		fmt.Println("calculate pb")
		tens, rho := NewPbTetraborateMaterialTensor()
		CalcFigureOfMerit(NewPbTetraborateEps(), NewPbTetraboratePhotoelasticTensor(), tens, rho, NewVec3ByThetaPhi(Rad(89.9), Rad(0)), lambda)


	default:
		fmt.Println("can calculate only teo2, linbo3 or te")
	}
}

func CreateFile(flnm string) *os.File {
	dest, err := os.Create(flnm)
	
	if err != nil {
		log.Panic(err)
	}

	return dest
}

//calculate figure of acousto-optical merit
/*
New classification of diffraction:
0123 - four points vertical axis intersection
12 - points Within, so it is "W" points
00 - points outZide, it it is "Z" points

There are 12 pairs of indexes. Let us organize them in 6 groups: izz, iww, aszw, aswz, alzw, alwz
where i - isotropic diffraction, as - anisotropic short diffraction, al - anisotropic long diffraction
in terms of point number:
izz  = 03 + 30
iww  = 12 + 21
aszw = 01 + 32
aswz = 10 + 23
alzw = 02 + 31
alwz = 13 + 20

left parts of these equations serve as filenames. Right parts define indexes of start end points for FM calculations
*/
func CalcFigureOfMerit(epsMainAxis Mat3, peMainAxis, mtMainAxis Tensor4, rho complex128, dir_a Vec3, lambda complex128) {
	destAcousticModes, err := os.Create("modes_description.txt")	
	if err != nil {
		log.Panic(err)
	}
	defer destAcousticModes.Close()

	transMat := MatrixByZ(dir_a)
	revMat := Transpose(transMat)

	eps := Mat3Rotation(epsMainAxis, transMat)
	pe := Tensor4Rotation(peMainAxis, transMat)

	m_ac := Vec3{0, 0, 1}
	//acousticData := CalcVelocities(m_ac, transMat, mtMainAxis, rho) 
	acousticData, _ := CalcPiezoVelocities(m_ac, transMat)

	for m_n := 0; m_n < 3; m_n ++ {
		fmt.Println("selected acoustic mode: ", acousticData[m_n])
		fmt.Println("\tacoustic polarization in crystal system:", Vec3Rotation(acousticData[m_n].Q, revMat))
		fmt.Println()

		fmt.Fprintln(destAcousticModes, "acoustic mode number ", m_n, "\tvelocity = ", acousticData[m_n].Vel)

		mIndexToFrequency := complex(acousticData[m_n].Vel, 0) / lambda

		/*
		nameTest := fmt.Sprintf("%d%s", m_n, "_test_ww_0.txt")
		destTest, err := os.Create(nameTest)
		if err != nil {
			log.Panic(err)
		}
		defer destTest.Close()
*/

		//create file descriptors for acousto-optic figure of merit
		destMizz := CreateFile(fmt.Sprintf("%d%s", m_n, "_m_izz.txt"))
		defer destMizz.Close()

		destMiww := CreateFile(fmt.Sprintf("%d%s", m_n, "_m_iww.txt"))
		defer destMiww.Close()

		destMaszw := CreateFile(fmt.Sprintf("%d%s", m_n, "_m_aszw.txt"))
		defer destMaszw.Close()
 
		destMaswz := CreateFile(fmt.Sprintf("%d%s", m_n, "_m_aswz.txt"))
		defer destMaswz.Close()

		destMalzw := CreateFile(fmt.Sprintf("%d%s", m_n, "_m_alzw.txt"))
		defer destMalzw.Close()
 
		destMalwz := CreateFile(fmt.Sprintf("%d%s", m_n, "_m_alwz.txt"))
		defer destMalwz.Close()
 
		//create file descriptors for the frequencies 
		destFreqIzz := CreateFile(fmt.Sprintf("%d%s", m_n, "_freq_izz.txt"))
		defer destFreqIzz.Close()

		destFreqIww := CreateFile(fmt.Sprintf("%d%s", m_n, "_freq_iww.txt"))
		defer destFreqIww.Close()

		destFreqAszw := CreateFile(fmt.Sprintf("%d%s", m_n, "_freq_aszw.txt"))
		defer destFreqAszw.Close()

		destFreqAswz := CreateFile(fmt.Sprintf("%d%s", m_n, "_freq_aswz.txt"))
		defer destFreqAswz.Close()

		destFreqAlzw := CreateFile(fmt.Sprintf("%d%s", m_n, "_freq_alzw.txt"))
		defer destFreqAlzw.Close()

		destFreqAlwz := CreateFile(fmt.Sprintf("%d%s", m_n, "_freq_alwz.txt"))
		defer destFreqAlwz.Close()

		//create file descriptors for 3d sphere representation
		destX := CreateFile("coord_x.txt")
		defer destX.Close()

		destY := CreateFile("coord_y.txt")
		defer destY.Close()

		destZ := CreateFile("coord_z.txt")
		defer destZ.Close()

		destSphereMizz := CreateFile(fmt.Sprintf("sphere_%d%s", m_n, "_m_izz.txt"))
		defer destSphereMizz.Close()

		destSphereMiww := CreateFile(fmt.Sprintf("sphere_%d%s", m_n, "_m_iww.txt"))
		defer destSphereMiww.Close()

		destSphereMaszw := CreateFile(fmt.Sprintf("sphere_%d%s", m_n, "_m_aszw.txt"))
		defer destSphereMaszw.Close()
 
		destSphereMaswz := CreateFile(fmt.Sprintf("sphere_%d%s", m_n, "_m_aswz.txt"))
		defer destSphereMaswz.Close()

		destSphereMalzw := CreateFile(fmt.Sprintf("sphere_%d%s", m_n, "_m_alzw.txt"))
		defer destSphereMalzw.Close()
 
		destSphereMalwz := CreateFile(fmt.Sprintf("sphere_%d%s", m_n, "_m_alwz.txt"))
		defer destSphereMalwz.Close()
		


		n := 200
		for p := 0; p < n; p++ {
			for q := 0; q <= n; q++	{
				theta_o := (math.Pi) * float64(p) / float64(n-1) 
				phi_o := 2 * math.Pi * float64(q) / float64(n)
				dir_o := NewVec3ByThetaPhi(theta_o, phi_o)

				fmt.Fprint(destX, real(dir_o[0]))
				fmt.Fprint(destY, real(dir_o[1]))
				fmt.Fprint(destZ, real(dir_o[2]))

				if q != n {
					fmt.Fprint(destX, ",")
					fmt.Fprint(destY, ",")
					fmt.Fprint(destZ, ",")
				}
				
				directionIndexes := IndexByDirection(dir_o, eps)
				if cmplx.Abs(directionIndexes[0].N) > cmplx.Abs(directionIndexes[1].N) {
					directionIndexes[0], directionIndexes[1] = directionIndexes[1], directionIndexes[0]
				}
				
				val := directionIndexes[0] //inner cavity
				{
					nx, ny, _ := dir_o[0] * val.N, dir_o[1] * val.N, dir_o[2] * val.N
					indexes := IndexByProjection(nx, ny, eps)
					sort.Sort(indexes)
					
					rn := 0
					for _, val_p := range(indexes) {
						if math.Abs(imag(val_p.N3)) < 1e-5 {
							rn ++
						}
					}

					if rn==4 {
						//iww group 12+21
						{
							//1-2 iww in-in
							fm1, freq1, ir1 := FinalFmCalc(nx, ny, *indexes[1], *indexes[2], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
							incDir := Vec3{nx, ny, indexes[1].N3}
							theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())
							
							real_theta_inc1, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
							if ir1 {
								fmt.Fprintf(destMiww, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc1), real(fm1))
								fmt.Fprintf(destFreqIww, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc1), real(freq1))
							}

							//2-1 iww in-in
							fm2, freq2, ir2 := FinalFmCalc(nx, ny, *indexes[2], *indexes[1], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
							incDir = Vec3{nx, ny, indexes[2].N3}
							theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

							real_theta_inc2, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
							if ir2 {
								fmt.Fprintf(destMiww, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc2), real(fm2))
								fmt.Fprintf(destFreqIww, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc2), real(freq2))
							}

							if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
								fmt.Fprint(destSphereMiww, real(fm1))
							} else {
								fmt.Fprint(destSphereMiww, real(fm2))
							}
							if q!=n {
								fmt.Fprint(destSphereMiww, ", ")
							}
						}


						//aswz group 10+23
						{
							//1-0 aswz in-out
							fm1, freq1, ir1 := FinalFmCalc(nx, ny, *indexes[1], *indexes[0], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
							incDir := Vec3{nx, ny, indexes[1].N3}
							theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

							real_theta_inc1, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
							if ir1 {
								fmt.Fprintf(destMaswz, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc1), real(fm1))
								fmt.Fprintf(destFreqAswz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc1), real(freq1))
							}
							
							//2-3 aswz in-out
							fm2, freq2, ir2 := FinalFmCalc(nx, ny, *indexes[2], *indexes[3], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
							incDir = Vec3{nx, ny, indexes[2].N3}
							theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

							real_theta_inc2, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
							if ir2 {
								fmt.Fprintf(destMaswz, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc2), real(fm2))
								fmt.Fprintf(destFreqAswz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc2), real(freq2))
							}

							if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
								fmt.Fprint(destSphereMaswz, real(fm1))
							} else {
								fmt.Fprint(destSphereMaswz, real(fm2))
							}
							if q!=n {
								fmt.Fprint(destSphereMaswz, ", ")
							}
						}						

						//alwz group 13+20
						{
							//1-3 alwz in-out
							fm1, freq1, ir1 := FinalFmCalc(nx, ny, *indexes[1], *indexes[3], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
							incDir := Vec3{nx, ny, indexes[1].N3}
							theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

							real_theta_inc1, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
							if ir1 {
								fmt.Fprintf(destMalwz, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc1), real(fm1))
								fmt.Fprintf(destFreqAlwz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc1), real(freq1))
							}

							//2-0 alwz in-out
							fm2, freq2, ir2 := FinalFmCalc(nx, ny, *indexes[2], *indexes[0], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
							incDir = Vec3{nx, ny, indexes[2].N3}
							theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

							real_theta_inc2, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
							if ir2 {
								fmt.Fprintf(destMalwz, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc2), real(fm2))
								fmt.Fprintf(destFreqAlwz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc2), real(freq2))
							}

							if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
								fmt.Fprint(destSphereMalwz, real(fm1))
							} else {
								fmt.Fprint(destSphereMalwz, real(fm2))
							}
							if q!=n {
								fmt.Fprint(destSphereMalwz, ", ")
							}
						}

					}  else {
						//log.Panic("number of real roots = ", rn)
					}				


					val := directionIndexes[1] 
					{
						nx, ny, _ := dir_o[0] * val.N, dir_o[1] * val.N, dir_o[2] * val.N
						indexes := IndexByProjection(nx, ny, eps)
						sort.Sort(indexes)
						
						rn := 0
						for _, val_p := range(indexes) {
							if math.Abs(imag(val_p.N3)) < 1e-5 {
								rn ++
							}
						}


						if rn==4 {
							//aszw group 01+32
							{
								//0-1 aszw out-in
								fm1, freq1, ir1 := FinalFmCalc(nx, ny, *indexes[0], *indexes[1], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir := Vec3{nx, ny, indexes[0].N3}
								theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
								if ir1 {
									fmt.Fprintf(destMaszw, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc), real(fm1))
									fmt.Fprintf(destFreqAszw, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq1))
								}

								//3-2 aszw out-in
								fm2, freq2, ir2 := FinalFmCalc(nx, ny, *indexes[3], *indexes[2], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir = Vec3{nx, ny, indexes[3].N3}
								theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc = CalcRealThetaPhi(incDir, revMat)
								if ir2 {
									fmt.Fprintf(destMaszw, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc), real(fm2))
									fmt.Fprintf(destFreqAszw, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq2))
								}

								if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
									fmt.Fprint(destSphereMaszw, real(fm1))
								} else {
									fmt.Fprint(destSphereMaszw, real(fm2))
								}
								if q!=n {
									fmt.Fprint(destSphereMaszw, ", ")
								}
							}

							//izz group 03+30
							{
								//0-3 izz out-out
								fm1, freq1, ir1 := FinalFmCalc(nx, ny, *indexes[0], *indexes[3], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir := Vec3{nx, ny, indexes[0].N3}
								theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
								if ir1 {							
									fmt.Fprintf(destMizz, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc), real(fm1))
									fmt.Fprintf(destFreqIzz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq1))
								}

								//3-0 izz out-out
								fm2, freq2, ir2 := FinalFmCalc(nx, ny, *indexes[3], *indexes[0], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir = Vec3{nx, ny, indexes[3].N3}
								theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc = CalcRealThetaPhi(incDir, revMat)
								if ir2 {
									fmt.Fprintf(destMizz, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc), real(fm2))
									fmt.Fprintf(destFreqIzz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq2))
								}

								if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
									fmt.Fprint(destSphereMizz, real(fm1))
								} else {
									fmt.Fprint(destSphereMizz, real(fm2))
								}
								if q!=n {
									fmt.Fprint(destSphereMizz, ", ")
								}
							}

							//alzw group 02+31
							{
								//0-2 alzw out-in
								fm1, freq1, ir1 := FinalFmCalc(nx, ny, *indexes[0], *indexes[2], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir := Vec3{nx, ny, indexes[0].N3}
								theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
								if ir1 {
									fmt.Fprintf(destMalzw, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc), real(fm1))
									fmt.Fprintf(destFreqAlzw, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq1))
								}
								
								//3-1 alzw out-in
								fm2, freq2, ir2 := FinalFmCalc(nx, ny, *indexes[3], *indexes[1], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir = Vec3{nx, ny, indexes[3].N3}
								theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc = CalcRealThetaPhi(incDir, revMat)
								if ir2 {
									fmt.Fprintf(destMalzw, "%e, %e, %e\n",  Deg(real_phi_inc), Deg(real_theta_inc), real(fm2))
									fmt.Fprintf(destFreqAlzw, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq2))
								}

								if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
									fmt.Fprint(destSphereMalzw, real(fm1))
								} else {
									fmt.Fprint(destSphereMalzw, real(fm2))
								}
								if q!=n {
									fmt.Fprint(destSphereMalzw, ", ")
								}

							}
						} else if rn==2 {
							flagFirst := true
							var source,dest IndexRetProj
							for _, val_p := range(indexes) {
								if math.Abs(imag(val_p.N3)) < 1e-5 {
									if flagFirst {
										flagFirst = false
										source = *val_p
									} else {
										dest = *val_p
									}
								}
							}
							
							//fmt.Fprintf(destTest, "%e, %e, %e\n", real(nx), real(ny), real(source.N3))
							//fmt.Fprintf(destTest, "%e, %e, %e\n", real(nx), real(ny), real(dest.N3))
							
							//izz group
							{
								//izz 0-1 out-out
								fm1, freq1, ir1 := FinalFmCalc(nx, ny, source, dest, eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir := Vec3{nx, ny, source.N3}
								theta_inc1, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())

								real_theta_inc, real_phi_inc := CalcRealThetaPhi(incDir, revMat)
								if ir1  {
									fmt.Fprintf(destMizz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(fm1))
									fmt.Fprintf(destFreqIzz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq1))
								}
						

								//izz 1-0 out-out						
								fm2, freq2, ir2 := FinalFmCalc(nx, ny, dest, source, eps, pe, acousticData[m_n], m_ac, real(rho), transMat, mIndexToFrequency)
								incDir = Vec3{nx, ny, source.N3}
								theta_inc2, _ := CalcRealThetaPhi(incDir, NewUnityMatrix())
								
								real_theta_inc, real_phi_inc = CalcRealThetaPhi(incDir, revMat)
								if ir2 {
									fmt.Fprintf(destMizz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(fm2))
									fmt.Fprintf(destFreqIzz, "%e, %e, %e\n", Deg(real_phi_inc), Deg(real_theta_inc), real(freq2))
								}

								if math.Abs(theta_o - theta_inc1) < math.Abs(theta_o - theta_inc2) {
									fmt.Fprint(destSphereMizz, real(fm1))
								} else {
									fmt.Fprint(destSphereMizz, real(fm2))
								}
								if q!=n {
									fmt.Fprint(destSphereMizz, ", ")
								}

							}
						}
					}
				}
			}

			fmt.Fprintln(destX)
			fmt.Fprintln(destY)
			fmt.Fprintln(destZ)

			fmt.Fprintln(destSphereMiww)
			fmt.Fprintln(destSphereMaswz)
			fmt.Fprintln(destSphereMalwz)

			fmt.Fprintln(destSphereMizz)
			fmt.Fprintln(destSphereMaszw)
			fmt.Fprintln(destSphereMalzw)
		}
	}
}

/*

//calculate figure of acousto-optical merit
func CalcFigureOfMerit(epsMainAxis Mat3, peMainAxis, mtMainAxis Tensor4, rho complex128, dir_a Vec3) {
	destAcousticModes, err := os.Create("modes_description.txt")	
	if err != nil {
		log.Panic(err)
	}
	defer destAcousticModes.Close()

	destThetaTest, err := os.Create("test_theta.txt")
	if err != nil {
		log.Panic(err)
	}
	defer destThetaTest.Close()

	if err != nil {
		log.Panic(err)
	}

	//theta_a := Rad(89.999) 
	//phi_a := Rad(60.0)
	//dir_a := Vec3{complex(math.Cos(phi_a) * math.Sin(theta_a), 0), complex(math.Sin(phi_a) * math.Sin(theta_a), 0), complex(math.Cos(theta_a), 0)}
	fmt.Println("dir_a = ", dir_a)
	
	transMat := MatrixByZ(dir_a)
	revMat := Transpose(transMat)

	fmt.Println("transMat = ")
	MatrixPPrint(transMat, "")
	fmt.Println("revMat = ", revMat)

	fmt.Println("epsMainAxis = ", epsMainAxis)

	eps := Mat3Rotation(epsMainAxis, transMat)
	pe := Tensor4Rotation(peMainAxis, transMat)

	fmt.Println("eps = ", eps)

	// peTest := Tensor4Rotation(pe, Transpose(transMat))
	// fmt.Println()
	// for t1 := 0; t1 < 3; t1 ++ {
	// 	for t2 := 0; t2 < 3; t2 ++ {
	// 		for t3 := 0; t3 < 3; t3 ++ {
	// 			for t4 := 0; t4 < 3; t4 ++ {
	// 				fmt.Println(t1+1, t2+1, t3+1, t4+1, real(peTest[t1][t2][t3][t4]))
	// 			}
	// 		}
	// 	}
	// }					
	// fmt.Println()


	testEps := Mat3Rotation(eps, Transpose(transMat))
	fmt.Println("test epsilon = ", testEps, "\n")
	fmt.Println("mat = ", transMat)

	fmt.Println("epsilon main axis", epsMainAxis, "\n")
	fmt.Println("epsilon rotated ", eps, "\n")

	//IndexByProjection(0.1, 0.1, eps)


	//mtMainAxis, rho := NewLinbo3MaterialTensor() 

	m_ac := Vec3{0, 0, 1}
	acousticData := CalcVelocities(m_ac, transMat, mtMainAxis, rho) //CalcPiezoVelocities(m_ac, transMat)

	for m_n := 0; m_n < 3; m_n ++ {
		fmt.Println("selected acoustic mode: ", acousticData[m_n])
		fmt.Println("\tacoustic polarization in crystal system:", Vec3Rotation(acousticData[m_n].Q, revMat), "\n")

		fmt.Fprintln(destAcousticModes, "acoustic mode number ", m_n, "\tvelocity = ", acousticData[m_n].Vel)

		var nameExtraordinary bytes.Buffer
		fmt.Fprint(&nameExtraordinary, m_n, "_extraordinary.txt")
		destExtraordinary, err := os.Create(nameExtraordinary.String())
		if err != nil {
			log.Panic(err)
		}
		defer destExtraordinary.Close()

		var nameOrdinary bytes.Buffer
		fmt.Fprint(&nameOrdinary, m_n, "_ordinary.txt")
		destOrdinary, err := os.Create(nameOrdinary.String())
		if err != nil {
			log.Panic(err)
		}
		defer destOrdinary.Close()

		debug := false

		n := 180
		for p := 0; p < n; p++ {
			for q := 0; q <= n; q++	{
				//theta_o := Rad(80)
				//phi_o := Rad(40.0)
				theta_o := math.Pi * float64(p) / float64(n) 
				phi_o := 2 * math.Pi * float64(q) / float64(n) 
		                //dir_o := Vec3Rotation(NewVec3ByThetaPhi(theta_o, phi_o), transMat)
				dir_o := NewVec3ByThetaPhi(theta_o, phi_o)
				if debug {
					fmt.Println("\t\t\tdir_o = ", dir_o)
					fmt.Println("\t\t\ttheta, phi = ", Deg(theta_o), Deg(phi_o))
					fmt.Println()
				}

				for iter_num, val := range(IndexByDirection(dir_o, eps)) {
					if debug {
						fmt.Println("\t\titernum = ", iter_num)
						fmt.Println("\t\t\tBy direction val = ", val)
					}

					nx, ny, nz := dir_o[0] * val.N, dir_o[1] * val.N, dir_o[2] * val.N

					if debug {
						fmt.Println("\t\t\tnx, ny, nz = ", nx, ny, nz)
					}

					indexes := IndexByProjection(nx, ny, eps)
					if debug {
						fmt.Println("\t\t\tindex by projection = ")
					}

					sort.Sort(indexes)
					rn := 0
					for _, val_p := range(indexes) {
						if math.Abs(imag(val_p.N3)) < 1e-5 {
							rn ++
						}
						//fmt.Println("\t\t\t\t", *val_p)
					}
					//fmt.Println()

					if rn==4 {
						//fmt.Println("\tfour real indexes")

						fmt.Fprintln(destThetaTest, theta_o, indexes[0].QD, indexes[3].QD)

						//fmt.Println("\t\t\tour debug nx, ny = ", nx, ny)

						fm := FinalFmCalc(nx, ny, *indexes[1], *indexes[2], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, debug)						
						real_theta_o, real_phi_o := CalcRealThetaPhi(Vec3{nx, ny, indexes[1].N3}, revMat)
						real_theta_o_d, real_phi_o_d := CalcRealThetaPhi(Vec3{nx, ny, indexes[2].N3}, revMat)

						fmt.Fprintf(destExtraordinary, "%e %e %e\n",  Deg(real_phi_o), Deg(real_theta_o), real(fm))
						if debug {
							fmt.Println("\t\t\t\tincident theta, phi, fm = ",Deg(real_theta_o), Deg(real_phi_o), real(fm))
							fmt.Println("\t\t\t\tdiffracted theta, phi, fm = ",Deg(real_theta_o_d), Deg(real_phi_o_d), real(fm))
							fmt.Println("\t\t\t\textra case theta_o =", Deg(real_theta_o), " phi_o =  ", Deg(real_phi_o), " fm = ", real(fm), " n1 = ", Abs(Vec3{nx, ny, indexes[1].N3}), " n2 = ", Abs(Vec3{nx, ny, indexes[2].N3}) )
							fmt.Println("\t\t\t\tpolarization from = ", indexes[1].QD)
							fmt.Println("\t\t\t\tpolarization to = ", indexes[2].QD)
							fmt.Println("\t\t\tfm = ", fm)
							fmt.Println()
						}

						fm = FinalFmCalc(nx, ny, *indexes[2], *indexes[1], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, debug)
						real_theta_o, real_phi_o = CalcRealThetaPhi(Vec3{nx, ny, indexes[2].N3}, revMat)
						real_theta_o_d, real_phi_o_d = CalcRealThetaPhi(Vec3{nx, ny, indexes[1].N3}, revMat)
						fmt.Fprintf(destExtraordinary, "%e %e %e\n", Deg(real_phi_o), Deg(real_theta_o), real(fm))
						if debug {
							fmt.Println("\t\t\t\tincident theta, phi, fm = ", Deg(real_theta_o), Deg(real_phi_o), real(fm))
							fmt.Println("\t\t\t\tdiffracted theta, phi, fm = ",Deg(real_theta_o_d), Deg(real_phi_o_d), real(fm))
							fmt.Println("\t\t\t\textra case  theta_o = ", Deg(real_theta_o), "phi_o = ", Deg(real_phi_o), " fm = ", real(fm), " n1 = ", Abs(Vec3{nx, ny, indexes[2].N3}), " n2 = ", Abs(Vec3{nx, ny, indexes[1].N3}))
							fmt.Println("\t\t\t\tpolarization from = ", indexes[2].QD)
							fmt.Println("\t\t\t\tpolarization to = ", indexes[1].QD)
							fmt.Println("\t\t\tfm = ", fm)
							fmt.Println()
						}


						fm = FinalFmCalc(nx, ny, *indexes[0], *indexes[3], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, debug)
						real_theta_o, real_phi_o = CalcRealThetaPhi(Vec3{nx, ny, indexes[0].N3}, revMat)
						real_theta_o_d, real_phi_o_d = CalcRealThetaPhi(Vec3{nx, ny, indexes[3].N3}, revMat)

						fmt.Fprintf(destOrdinary,  "%e %e %e\n", Deg(real_phi_o), Deg(real_theta_o), real(fm))
						if debug {
							fmt.Println("\t\t\t\tincident theta, phi, fm = ",Deg(real_theta_o), Deg(real_phi_o), real(fm))
							fmt.Println("\t\t\t\tdiffracted theta, phi, fm = ",Deg(real_theta_o_d), Deg(real_phi_o_d), real(fm))
							fmt.Println("\t\t\tord case theta_o = ", Deg(real_theta_o), " phi_o = ", Deg(real_phi_o), " fm = ", real(fm), " n1 = ", Abs(Vec3{nx, ny, indexes[0].N3}), " n2 = ", Abs(Vec3{nx, ny, indexes[3].N3}))
							fmt.Fprintln(destThetaTest, Deg(real_theta_o), real(fm))
							fmt.Println("\t\t\t\tpolarization from = ", indexes[0].QD)
							fmt.Println("\t\t\t\tpolarization to = ", indexes[3].QD,)
							fmt.Println("\t\t\tfm = ", fm)
							fmt.Println()
						}

						fm = FinalFmCalc(nx, ny, *indexes[3], *indexes[0], eps, pe, acousticData[m_n], m_ac, real(rho), transMat, debug)
						real_theta_o, real_phi_o = CalcRealThetaPhi(Vec3{nx, ny, indexes[3].N3}, revMat)
						real_theta_o_d, real_phi_o_d = CalcRealThetaPhi(Vec3{nx, ny, indexes[0].N3}, revMat)

						fmt.Fprintf(destOrdinary,  "%e %e %e\n", Deg(real_phi_o), Deg(real_theta_o), real(fm))
						if debug {
							fmt.Println("\t\t\t\tincident theta, phi, fm = ",Deg(real_theta_o), Deg(real_phi_o), real(fm))
							fmt.Println("\t\t\t\tdiffracted theta, phi, fm = ",Deg(real_theta_o_d), Deg(real_phi_o_d), real(fm))
							fmt.Println("\t\t\t\tord case theta_o = ", Deg(real_theta_o), " phi_o = ", Deg(real_phi_o), " fm = ", real(fm), " n1 = ", Abs(Vec3{nx, ny, indexes[3].N3}), " n2 = ", Abs(Vec3{nx, ny, indexes[0].N3}))
							fmt.Fprintln(destThetaTest, Deg(real_theta_o), real(fm))
							fmt.Println("\t\t\t\tpolarization from = ", indexes[3], revMat)
							fmt.Println("\t\t\t\tpolarization to = ", indexes[0], revMat)
							fmt.Println("\t\t\tfm = ", fm)
							fmt.Println()
						}


					} else if rn==2 {
						//fmt.Println("two real roots")
						flagFirst := true
						var source,dest IndexRetProj
						for _, val_p := range(indexes) {
							if math.Abs(imag(val_p.N3)) < 1e-5 {
								if flagFirst {
									flagFirst = false
									source = *val_p
								} else {
									dest = *val_p
								}
							}
						}
						fm := FinalFmCalc(nx, ny, source, dest, eps, pe, acousticData[m_n], m_ac, real(rho), transMat, debug)
						real_theta_o, real_phi_o := CalcRealThetaPhi(Vec3{nx, ny, source.N3}, revMat)
						real_theta_o_d, real_phi_o_d := CalcRealThetaPhi(Vec3{nx, ny, dest.N3}, revMat)

						fmt.Fprintf(destOrdinary, "%e %e %e\n", Deg(real_phi_o), Deg(real_theta_o), real(fm))
						if debug {
							fmt.Println("\t\t\t\tincident theta, phi, fm = ",Deg(real_theta_o), Deg(real_phi_o), real(fm))
							fmt.Println("\t\t\t\tdiffracted theta, phi, fm = ",Deg(real_theta_o_d), Deg(real_phi_o_d), real(fm))
							fmt.Println("\t\t\t\tord case theta_o = ", Deg(real_theta_o), " phi_o = ", Deg(real_phi_o), " fm = ", real(fm), " n1 = ", Abs(Vec3{nx, ny, source.N3}), " n2 = ", Abs(Vec3{nx, ny, dest.N3}))
							fmt.Println("\t\t\t\tpolarization from = ", source.QD)
							fmt.Println("\t\t\t\tpolarization to = ", dest.QD)
							fmt.Println("\t\t\tfm = ", fm)
							fmt.Println()
						}

						fm = FinalFmCalc(nx, ny, dest, source, eps, pe, acousticData[m_n], m_ac, real(rho), transMat, debug)
						real_theta_o, real_phi_o = CalcRealThetaPhi(Vec3{nx, ny, dest.N3}, revMat)
						real_theta_o_d, real_phi_o_d = CalcRealThetaPhi(Vec3{nx, ny, source.N3}, revMat)

						fmt.Fprintf(destOrdinary, "%e %e %e\n", Deg(real_phi_o), Deg(real_theta_o), real(fm))
						if debug {
							fmt.Println("\t\t\t\tincident theta, phi, fm = ",Deg(real_theta_o), Deg(real_phi_o), real(fm))
							fmt.Println("\t\t\t\tdiffracted theta, phi, fm = ",Deg(real_theta_o_d), Deg(real_phi_o_d), real(fm))
							//fmt.Println("\t\t\t\tord case theta_o = ", Deg(real_theta_o), " phi_o = ", Deg(real_phi_o), " fm = ", real(fm), " n1 = ", Abs(Vec3{nx, ny, source.N3}), " n2 = ", Abs(Vec3{nx, ny, dest.N3}))
							fmt.Println("\t\t\t\tpolarization from = ", dest.QD)
							fmt.Println("\t\t\t\tpolarization to = ", source.QD)
							fmt.Println("\t\t\tfm = ", fm)
							fmt.Println()
						}


					} else {
						log.Panic("number of real roots = ", rn)
					}
					//fmt.Println()
					//fmt.Println()
				}
			}
			//fmt.Fprintln(destOrdinary)
			//fmt.Fprintln(destExtraordinary)
		}
	}
}

*/

//radians to degrees
func Deg(x float64) float64 {
	return x * 180.0 / math.Pi;
}

//degrees to radians
func Rad(x float64) float64 {
	return x * math.Pi / 180
}

func testRotateBack(sourceRotPol, destRotPol Vec3, epsRot Mat3, peRot Tensor4, acoustRotPol Vec3, m_acRot Vec3, transMat Mat3) {
	
	fmt.Println("testRotateBack")

	revMat := Transpose(transMat)

	sourcePol := Vec3Rotation(sourceRotPol, revMat)
	destPol := Vec3Rotation(destRotPol, revMat)
	eps := Mat3Rotation(epsRot, revMat)
	pe := Tensor4Rotation(peRot, revMat)
	acoustPol := Vec3Rotation(acoustRotPol, revMat)
	m_ac := Vec3Rotation(m_acRot, revMat)
	
	fmt.Println()
	fmt.Println("\t\t\ttestRotateBack")
	fmt.Println("\t\t\tsourcePol = ", sourcePol)
	fmt.Println("\t\t\tdestPol = ", destPol)
	fmt.Println("\t\t\teps = ")
	MatrixPPrint(eps, "\t\t\t\t")
	//fmt.Println("\t\t\tpe = ", pe)
	//tensor4PPrint(os.Stdout, pe, "\t\t\t")
	fmt.Println("\t\t\tacoustPol = ", acoustPol)
	fmt.Println("\t\t\tm_ac = ", m_ac)
	fmt.Println()

	var leftPart, rightPart Mat3
	for j := 0; j < 3; j++ {
		for k := 0; k < 3; k++ {
			for i := 0; i < 3; i++ {
				for l := 0; l < 3; l++ {
					leftPart[j][k] += sourcePol[i]*eps[i][j]*eps[k][l]*destPol[l]
				}
			}
		}
	}

	for j := 0; j < 3; j++ {
		for k := 0; k < 3; k++ {
			for m := 0; m < 3; m++ {
				for n := 0; n < 3; n++ {
					rightPart[j][k] += pe[j][k][m][n] *
								(acoustPol[m] * m_ac[n] +
								acoustPol[n] * m_ac[m])
				}
			}
		}
	}

	fmt.Println("\t\t\tleftPart = ")
	MatrixPPrint(leftPart, "\t\t\t\t")
	fmt.Println("\t\t\trightPart = ")
	MatrixPPrint(rightPart, "\t\t\t\t")

	var testSumm complex128

	for p := 0; p < 3; p++ {
		for q := 0; q < 3; q++ {
			testSumm += leftPart[p][q] * rightPart[p][q]
		}
	}

	fmt.Println("\t\t\ttestSumm = ", testSumm)
	fmt.Println()

	fm := complex(0, 0)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					for m := 0; m < 3; m ++ {
						for n := 0; n < 3; n ++ {
							val := sourcePol[i]*eps[i][j]*eps[k][l]*
								destPol[l] * pe[j][k][m][n] *
								(acoustPol[m] * m_ac[n] +
								acoustPol[n] * m_ac[m])
							
							if (cmplx.Abs(val) > 1e-5) {
								fmt.Println("\t\t\t\t", i+1, j+1, k+1, l+1, m+1, n+1, val)
							}

							fm = fm + val 
						}
					}
				}
			}
		}
	}

	fmt.Println("\t\t\tsumm in RotateBack = ", fm)
	fmt.Println()
}

func FinalFmCalc(nx, ny complex128, source, dest IndexRetProj, eps Mat3, pe Tensor4, acousticData VelocityResult, m_ac Vec3, rho float64, transMat Mat3, mIndexToFrequency complex128) (fm, frequency complex128, ir bool) {

	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					for m := 0; m < 3; m ++ {
						for n := 0; n < 3; n ++ {

							val := source.QD[i]*eps[i][j]*eps[k][l]*
								dest.QD[l] * pe[j][k][m][n] *
								(acousticData.Q[m] * m_ac[n] +
								acousticData.Q[n] * m_ac[m])
							
							fm = fm + val 
						}
					}
				}
			}
		}
	}

	V := acousticData.Vel

	fm = fm * fm

	fm = fm / complex(4 * source.N * dest.N * rho * V * V * V, 0)
	
	frequency = complex(cmplx.Abs((source.N3 - dest.N3) * mIndexToFrequency), 0)

	if false && (cmplx.Abs(frequency) <10e6 || cmplx.Abs(frequency)>500e6) {
		frequency = 0
		ir = false
	} else {
		ir = true
	}

	frequency /= 1e6

	return 
}

func TestFresnel() {
	epsMainAxis := Linbo3eps()

	transMat := NewRotxMatrix(Rad(20))
	eps := Mat3Rotation(epsMainAxis, transMat)
	
	for theta := 0.0; theta < 95.0; theta += 1.0 {
		dir_o := NewVec3ByThetaPhi(Rad(theta), Rad(20))
		for _, val := range(IndexByDirection(dir_o, eps)) {
			if math.Abs(real(val.N) - 2.32) < 1e-5 {
				fmt.Println(theta, "val = ", val.N, val.QD, Vec3Rotation(val.QD, Transpose(transMat)))
			}
		}
	}
}

func MixedCase() {

	destOrdinary, err := os.Create("ordinary.txt")
	if err != nil {
		log.Panic(err)
	}
	defer destOrdinary.Close()

	destExtra, err := os.Create("extra.txt")
	if err != nil {
		log.Panic(err)
	}
	defer destExtra.Close()

	//for theta_grad := 0.0; theta_grad < 180.0; theta_grad ++ 
	{
		//for phi_grad := 0.0; phi_grad < 90.0; phi_grad += 1.0 
		{
			theta_a := Rad(13)
			phi_a := Rad(90)
			dir_a := Vec3{complex(math.Cos(phi_a) * math.Sin(theta_a), 0), complex(math.Sin(phi_a) * math.Sin(theta_a), 0), complex(math.Cos(theta_a), 0)}
			
			transMat := MatrixByZ(dir_a)
			revMat := Transpose(transMat)
			epsMainAxis := Linbo3eps()
			peMainAxis := NewLinbo3Photoelastic()

			eps := Mat3Rotation(epsMainAxis, transMat)
			acousticData, rho := CalcPiezoVelocities(NewVec3(0,0,1), transMat)
			for ind := range(acousticData) {
				acousticData[ind].Q = Vec3Rotation(acousticData[ind].Q, revMat)
			}

			n_ac := 0
			modeData := acousticData[n_ac]

			theta_o, phi_o := Rad(103), Rad(90)
			dir_o := Vec3Rotation(NewVec3ByThetaPhi(theta_o, phi_o), transMat)
			for _, val := range(IndexByDirection(dir_o, eps)) {
				nx, ny, nz := dir_o[0] * val.N, dir_o[1] * val.N, dir_o[2] * val.N
				trueSource := Vec3Rotation(Vec3{nx, ny, nz}, revMat)
				indexes := IndexByProjection(nx, ny, eps)
				sort.Sort(indexes)
				rn := 0
				for _, val_p := range(indexes) {
					if math.Abs(imag(val_p.N3)) < 1e-5 {
						rn ++
					} 
				}

				if rn == 4 {
					ProcessFm(destOrdinary, trueSource, real(nx), real(ny), *indexes[0], *indexes[3], revMat, true , epsMainAxis, modeData, real(rho), peMainAxis, dir_a)
					ProcessFm(destExtra, trueSource, real(nx), real(ny), *indexes[1], *indexes[2], revMat, false, epsMainAxis, modeData, real(rho), peMainAxis, dir_a)
				}

				if rn==2 {
					flagFirst := true
					var source,dest IndexRetProj
					for _, val_p := range(indexes) {
						if math.Abs(imag(val_p.N3)) < 1e-5 {
							if flagFirst {
								flagFirst = false
								source = *val_p
							} else {
								dest = *val_p
							}
						}
					}
					ProcessFm(destOrdinary, trueSource, real(nx), real(ny), source, dest, revMat, true, epsMainAxis, modeData, real(rho), peMainAxis, dir_a)
				}
			}
		}
	}
}

func ProcessFm(destFile *os.File, trueSource Vec3, nx, ny float64, source IndexRetProj, dest IndexRetProj, revMat Mat3, isOrd bool, 
	eps Mat3, modeData VelocityResult, rho float64, pe Tensor4, dir_a Vec3) {
	sourceVec := Vec3Rotation(NewVec3(nx, ny, real(source.N3)), revMat)
	destVec := Vec3Rotation(NewVec3(nx, ny, real(dest.N3)), revMat)

	if Abs(VectorDiff(trueSource, destVec)) < 1e-6 {
		sourceVec, destVec = destVec, sourceVec
	} 

	if Abs(VectorDiff(trueSource, sourceVec)) > 1e-6 {
		return
	}
	fmt.Println()

	var sourcePol, destPol Vec3
	ok := false
	for _, val := range(IndexByDirection(Normalized(sourceVec), eps)) {
		if math.Abs(real(val.N) - Abs(sourceVec)) < 1e-6 {
			ok = true
			sourcePol = val.QD
		}
	}

	if !ok {
		log.Panic("not ok in ProcessFm")
	}

	ok = false
	for _, val := range(IndexByDirection(Normalized(destVec), eps)) {
		if math.Abs(real(val.N) - Abs(destVec)) < 1e-8 {
			ok = true
			destPol = val.QD
		}
	}

	if !ok {
		log.Panic("not ok in ProcessFm")
	}

/*
	if isOrd {
		fmt.Println("ordinary")
	} else {
		fmt.Println("extraordinary")
	}*/

	//fmt.Println("source polariz =", sourcePol)
	//fmt.Println("dest polariz = ", destPol)
	//fmt.Println()

	V := modeData.Vel
	fmMult := 1.0 / (4 * Abs(sourceVec) * Abs(destVec) * rho * V * V * V)

	fmt.Printf("sourcePol := %#v\n", sourcePol)
	fmt.Printf("destPol := %#v\n", destPol)
	fmt.Printf("eps := %#v\n", eps)
	fmt.Printf("pe := %#v\n", pe)
	fmt.Printf("dir_a := %#v\n", dir_a)
	fmt.Printf("pol_a := %#v\n", modeData.Q)

	summ := 0.0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					for m := 0; m < 3; m ++ {
						for n := 0; n < 3; n ++ {
							summ += real(sourcePol[i]) * real(eps[i][j]) * real(eps[k][l]) * real(destPol[l]) * 
								real(pe[j][k][m][n]) * real(dir_a[n] * modeData.Q[m] + dir_a[n] * modeData.Q[m])
						}
					}
				}
			}
		}
	}
	
	fmt.Println("summ = ", summ)
	fm := fmMult * summ * summ

	sourceTheta, sourcePhi := CalcRealThetaPhi(sourceVec, NewUnityMatrix())
	fmt.Fprintln(destFile, Deg(sourceTheta), Deg(sourcePhi), fm)
	
}

func CalcSumm(sourcePol, destPol, dir_a, pol_a Vec3, eps Mat3, pe Tensor4) (summ float64){
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					for m := 0; m < 3; m ++ {
						for n := 0; n < 3; n ++ {
							summ += real(sourcePol[i]) * real(eps[i][j]) * real(eps[k][l]) * real(destPol[l]) * 
								real(pe[j][k][m][n]) * real(dir_a[n] * pol_a[m] + dir_a[n] * pol_a[m])
						}
					}
				}
			}
		}
	}
	return
}

/*
	difference := VectorDiff(destVec, sourceVec)

	diffTheta, diffPhi := CalcRealThetaPhi(difference)


	fmt.Println("difference theta phi = ", Deg(diffTheta), Deg(diffPhi))
	fmt.Println("abs(source) = ", Abs(sourceVec))
	fmt.Println("abs(dest) = ", Abs(destVec))
	fmt.Println()

*/

func TestMeritRotation() {
	fmt.Println("TestMeritRotation")

	/*sourcePol := Normalized(NewVec3(1, 2, 3))
	destPol := Normalized(NewVec3(4, 5, 6))
	dir_a := Normalized(NewVec3(1, 0, 0))
	pol_a := Normalized(NewVec3(0, 1, 0))
	eps := Mat3{{11, 12, 13}, {21, 22, 23}, {31, 32, 33}}
	pe := NewTetragonalPhotoelasticTensor(1, 2, 3, 4, 5, 6, 7, 8)*/

	sourcePol := Vec3{(0+0i), (0.22495105434387122+0i), (0.9743700647852337+0i)}
	destPol := Vec3{(0+0i), (0.25981033005129545+0i), (0.9656596669627643+0i)}
	eps := Mat3{[3]complex128{(5.3824+0i), (0+0i), (0+0i)}, [3]complex128{(0+0i), (5.3824+0i), (0+0i)}, [3]complex128{(0+0i), (0+0i), (4.9729+0i)}}
	pe := Tensor4{[3][3][3]complex128{[3][3]complex128{[3]complex128{(-0.021+0i), (0+0i), (0+0i)}, [3]complex128{(0+0i), (0.06+0i), (-0.052+0i)}, [3]complex128{(0+0i), (-0.052+0i), (0.172+0i)}}, [3][3]complex128{[3]complex128{(0+0i), (-0.0405+0i), (-0.052+0i)}, [3]complex128{(-0.0405+0i), (0+0i), (0+0i)}, [3]complex128{(-0.052+0i), (0+0i), (0+0i)}}, [3][3]complex128{[3]complex128{(0+0i), (0.109+0i), (0.121+0i)}, [3]complex128{(0.109+0i), (0+0i), (0+0i)}, [3]complex128{(0.121+0i), (0+0i), (0+0i)}}}, [3][3][3]complex128{[3][3]complex128{[3]complex128{(0+0i), (-0.0405+0i), (-0.052+0i)}, [3]complex128{(-0.0405+0i), (0+0i), (0+0i)}, [3]complex128{(-0.052+0i), (0+0i), (0+0i)}}, [3][3]complex128{[3]complex128{(0.06+0i), (0+0i), (0+0i)}, [3]complex128{(0+0i), (-0.021+0i), (0.052-0i)}, [3]complex128{(0+0i), (0.052-0i), (0.172+0i)}}, [3][3]complex128{[3]complex128{(0.109+0i), (0+0i), (0+0i)}, [3]complex128{(0+0i), (-0.109-0i), (0.121+0i)}, [3]complex128{(0+0i), (0.121+0i), (0+0i)}}}, [3][3][3]complex128{[3][3]complex128{[3]complex128{(0+0i), (0.109+0i), (0.121+0i)}, [3]complex128{(0.109+0i), (0+0i), (0+0i)}, [3]complex128{(0.121+0i), (0+0i), (0+0i)}}, [3][3]complex128{[3]complex128{(0.109+0i), (0+0i), (0+0i)}, [3]complex128{(0+0i), (-0.109-0i), (0.121+0i)}, [3]complex128{(0+0i), (0.121+0i), (0+0i)}}, [3][3]complex128{[3]complex128{(0.141+0i), (0+0i), (0+0i)}, [3]complex128{(0+0i), (0.141+0i), (0+0i)}, [3]complex128{(0+0i), (0+0i), (0.118+0i)}}}}
	dir_a := Vec3{(1.377427943335181e-17+0i), (0.224951054343865+0i), (0.9743700647852352+0i)}
	pol_a := Vec3{(1.1835572873289225e-17+0i), (0.19348525243984363+0i), (0.981103183711219+0i)}

	fmt.Println("CalcSumm =", CalcSumm(sourcePol, destPol, dir_a, pol_a, eps, pe));	

	rot := NewRotxMatrix(Rad(-13));
	sourcePolr := Vec3Rotation(sourcePol, rot)
	destPolr := Vec3Rotation(destPol, rot)
	dir_ar := Vec3Rotation(dir_a, rot)
	pol_ar := Vec3Rotation(pol_a, rot)
	epsr := Mat3Rotation(eps, rot)
	per := Tensor4Rotation(pe, rot)
	fmt.Println("CalcSumm =", CalcSumm(sourcePolr, destPolr, dir_ar, pol_ar, epsr, per));	

}

func TestAngles() {
	theta_a := Rad(13.0) 
	phi_a := Rad(90.0)
	dir_a := Vec3{complex(math.Cos(phi_a) * math.Sin(theta_a), 0), complex(math.Sin(phi_a) * math.Sin(theta_a), 0), complex(math.Cos(theta_a), 0)}
	fmt.Println("dir_a = ", dir_a)
	
	transMat := MatrixByZ(dir_a)
	revMat := Transpose(transMat)
	epsMainAxis := Linbo3eps()
	peMainAxis := NewLinbo3Photoelastic()

	eps := Mat3Rotation(epsMainAxis, transMat)
	pe := Tensor4Rotation(peMainAxis, transMat)

	testEps := Mat3Rotation(eps, Transpose(transMat))

	m_ac := Vec3{0, 0, 1}
	acousticData, rho := CalcPiezoVelocities(m_ac, transMat)

	//fmt.Print("acoustic = ", acousticData)

	theta_o_c := Rad(45)
	phi_o_c := Rad(45)

	dir_o_c := NewVec3ByThetaPhi(theta_o_c, phi_o_c)
	dir_o_i_l := Vec3Rotation(dir_o_c, transMat)

	fmt.Println("crystal system: vec =", dir_o_c)
	fmt.Println("laboratory system: vec =", dir_o_i_l)

	for _, val := range(IndexByDirection(dir_o_i_l, eps)) {
		nx, ny, _ := dir_o_i_l[0] * val.N, dir_o_i_l[1] * val.N, dir_o_i_l[2] * val.N
		fmt.Println("\tnx, ny = ", nx, " ", ny)
		fmt.Println("\tN = ", val.N)
		indexes := IndexByProjection(nx, ny, eps)
		sort.Sort(indexes)

		for _, val_d := range(indexes) {
			dir_o_d_l := Vec3{nx, ny, val_d.N3}
			fmt.Println("\t\tN = ", Abs(dir_o_d_l), "\n")
			fmt.Println("\t\tdifracted in laboratory = ", dir_o_d_l)
			dir_o_d_c := Vec3Rotation(dir_o_d_l, revMat)
			fmt.Println("\t\tdifracted in crystal = ", dir_o_d_c)
			theta_d_c, phi_d_c := CalcRealThetaPhi(dir_o_d_c, NewUnityMatrix())
			fmt.Println("\t\ttheta, phi = ", Deg(theta_d_c), Deg(phi_d_c))
			pol_o_d_l := val_d.QD
			pol_o_d_c := Vec3Rotation(pol_o_d_l, revMat)
			fmt.Println("\t\tpolarization in crystal = ", pol_o_d_c)
			fmt.Println()
		}
	}

	_ = revMat
	_ = pe
	_ = testEps
	_ = rho
	_ = acousticData
}

func TestLoadPicture() {
	mat := LoadMatrix("a.png")

	for p := 0; p < len(mat); p++ {
		for q := 0; q<len(mat[0]); q ++ {
			if cmplx.Abs(mat[p][q]) < 1e-5 {
				fmt.Print(" ")
			} else {
				fmt.Print("*")
			}
		}
		fmt.Println()
	}
}

func TestSaveMatrix() {
	nx := 500
	ny := 600

	mat := fftw.Alloc2d(ny, nx)

	for p := 0; p < ny; p++ {
		for q := 0; q<nx; q ++ {
			mat[p][q] = complex(math.Sin(math.Sqrt(float64(p*p + q*q)) / 4), 0)
		}
	}

	SaveMatrix(mat, "sin_dreams.png", false)
}

func CalcParatellurite() {
	tens := NewTeo2Photoelastic()
	Tensor4SixPrint(os.Stdout, tens)
}

func TestTeO2() {
	dir := Normalized(Vec3{0,1,0})
	tens, rho := NewTeMaterialTensor()
	cm := Christoffel(dir, tens)
	fmt.Println("christoffel matrix is ", cm)
	fmt.Println("rho = ", rho)

	gp, _ := NewChristoffelPoly(tens, dir)
	ar := gp.AllRoots()
	for _, r := range ar {
		v := 1/cmplx.Sqrt(rho / r)
		fmt.Println("v = ", real(v)) 
	}
}

func TestCopy() {
	a := make([]int, 5)
	b := make([]int, 5)

	a[1] = 6
	a[2] = 7
	a[3] = 8

	fmt.Println("a = ", a)
	fmt.Println("b = ", b)

	copy(b,a)
	
	fmt.Println()
	fmt.Println("a = ", a)
	fmt.Println("b = ", b)

}

func NonMoveCopy(b, a [][]complex128) {
	for p := 0; p < len(a); p++ {
		for q := 0; q < len(a[0]); q++ {
			b[p][q] = a[p][q]
		}
	}
}

func TestConvolution() {
	fmt.Println("test convolution")
	apic := LoadMatrix("a.png")
	bpic := LoadMatrix("b.png")

	fmt.Println("apic module = ", MatrixAbs(apic))
	fmt.Println("bpic module = ", MatrixAbs(bpic))
	
	nx := len(apic[0])
	ny := len(apic)

	fftwSource := fftw.Alloc2d(ny, nx)
	fftwDest := fftw.Alloc2d(ny, nx)

	planFWD := fftw.PlanDft2d( fftwSource, fftwDest, fftw.Forward, fftw.Estimate)
	planBACK := fftw.PlanDft2d(fftwSource, fftwDest, fftw.Backward, fftw.Estimate)

	af := fftw.Alloc2d(ny, nx)
	bf := fftw.Alloc2d(ny, nx)
	cf := fftw.Alloc2d(ny, nx)
	cpic := fftw.Alloc2d(ny, nx)

	NonMoveCopy(fftwSource, apic)
	planFWD.Execute()
	NonMoveCopy(af, fftwDest)

	NonMoveCopy(fftwSource, bpic)
	planFWD.Execute()
	NonMoveCopy(bf, fftwDest)
	
	for p := 0; p < ny; p++ {
		for q := 0; q < nx; q++ {
			cf[p][q] = af[p][q] * bf[p][q]
		}
	}

	NonMoveCopy(fftwSource, cf)
	planBACK.Execute()
	NonMoveCopy(cpic, fftwDest)

	fmt.Println("cpic module = ", MatrixAbs(cpic))

	SaveMatrix(cpic, "c.png", false)
}

func CalcFieldStructure(f, Ax complex128, nz int, tens Tensor4, rho complex128, forceProfileFlnm string, forceDir Vec3, flnmTmpl string) {
	forceProfile := LoadMatrix(forceProfileFlnm) // n.b. forceProfile is slice of slices
	nx := len(forceProfile[0])
	ny := len(forceProfile)
	dz := Ax / complex(float64(nx), 0)
	Ay := dz * complex(float64(ny), 0)
	fmt.Println("force profile loaded", "nx=", nx, "ny=", ny)


	wm := NewWaveMatrix(ny, nx, tens, rho, Ax, Ay, f, forceDir) 
	fmt.Println("wave matrix calculated")

	fftwSource := fftw.Alloc2d(ny, nx)
	fftwDest := fftw.Alloc2d(ny, nx)
	planFWD := fftw.PlanDft2d(fftwSource, fftwDest, fftw.Forward, fftw.Estimate)
	planBACK := fftw.PlanDft2d(fftwSource, fftwDest, fftw.Backward, fftw.Estimate)

	fpFourier := fftw.Alloc2d(ny, nx)
	
	NonMoveCopy(fftwSource, forceProfile);
	planFWD.Execute()
	NonMoveCopy(fpFourier, fftwDest)

	wm.LoadProfileImage(fpFourier)

	ows := NewStorage(PlaneWave{}.GetNumberOfComponents(), ny, nx)

	volumeData := NewStorage(nz + 1, ny, nx)

	for s := 0; s <= nz; s++ {
		fmt.Println("0 < ", s, "<", nz)
		
		wm.FillStorage(ows)
		wm.ShiftZ(dz);

		ows.LayeredTransform(fftwSource, fftwDest, planBACK)

		//plane wave matrix
		pwm := NewPlaneWaveMatrix(ows)

		ds := pwm.GetDisplacementSlice()
		volumeData.SetZLayer(s, ds)

		//var bufFlnm bytes.Buffer
		//fmt.Fprint(&bufFlnm, "wd", s, ".png")
		//SaveMatrix(ds, bufFlnm.String(), false)
	}

	
	for s := 0; s < ny; s++ {
		fmt.Println(s, MatrixAbs(volumeData.SliceHW(s)))
	}

	SaveMatrix(volumeData.SliceHW(ny/2), flnmTmpl + "hw.png", false)
	SaveMatrix(volumeData.SliceHD(nx/2), flnmTmpl + "hd.png", false)
}

func TestLayeredTransform() {
	nx, ny := 200, 100
	h := 10

	fftwSource := fftw.Alloc2d(ny, nx)
	fftwDest := fftw.Alloc2d(ny, nx)
	planFWD := fftw.PlanDft2d(fftwSource, fftwDest, fftw.Forward, fftw.Estimate)
	//planBACK := fftw.PlanDft2d(fftwSource, fftwDest, fftw.Backward, fftw.Estimate)


	tstStor := NewStorage(h, ny, nx)

	for s := 0; s < h; s++ {
		for p := 0; p < ny; p++ {
			for q := 0; q < nx; q++ {
				dx := q - nx / 2
				dy := p - ny / 2
				if (dx*dx + dy*dy) < (10*s)*(10*s) {
					tstStor.Set(s,p,q,1)
				}
			}
		}
	}

	SaveMatrix(tstStor.GetMatrix(0),  "a00.png", false)
	SaveMatrix(tstStor.GetMatrix(1),  "a01.png", false)
	SaveMatrix(tstStor.GetMatrix(2),  "a02.png", false)
	SaveMatrix(tstStor.GetMatrix(3),  "a03.png", false)
	SaveMatrix(tstStor.GetMatrix(4),  "a04.png", false)
	SaveMatrix(tstStor.GetMatrix(5),  "a05.png", false)
	SaveMatrix(tstStor.GetMatrix(6),  "a06.png", false)
	SaveMatrix(tstStor.GetMatrix(7),  "a07.png", false)
	SaveMatrix(tstStor.GetMatrix(8),  "a08.png", false)
	SaveMatrix(tstStor.GetMatrix(9),  "a09.png", false)

	tstStor.LayeredTransform(fftwSource, fftwDest, planFWD)

	SaveMatrix(tstStor.GetMatrix(0),  "b00.png", false)
	SaveMatrix(tstStor.GetMatrix(1),  "b01.png", false)
	SaveMatrix(tstStor.GetMatrix(2),  "b02.png", false)
	SaveMatrix(tstStor.GetMatrix(3),  "b03.png", false)
	SaveMatrix(tstStor.GetMatrix(4),  "b04.png", false)
	SaveMatrix(tstStor.GetMatrix(5),  "b05.png", false)
	SaveMatrix(tstStor.GetMatrix(6),  "b06.png", false)
	SaveMatrix(tstStor.GetMatrix(7),  "b07.png", false)
	SaveMatrix(tstStor.GetMatrix(8),  "b08.png", false)
	SaveMatrix(tstStor.GetMatrix(9),  "b09.png", false)
}

func TestLinbPolarization() {
	mt, rho := NewLinbo3MaterialTensor()
	phi := Rad(13)
	dir := NewVec3(0, math.Sin(phi), math.Cos(phi))
	chr, chrMat := NewChristoffelPoly(mt, dir)
	allRoots := chr.AllRoots()
	
	for _, r := range(allRoots) {
		fmt.Println("r = ", r)
		fmt.Println("\tv = ", cmplx.Sqrt(r/rho))
		pme := PolyMatEval(chrMat, r)
		pol, _ := CalcPol(pme)
		fmt.Println("\tpol = ", pol)
	}

	_ = chrMat
	
}

type SortableArray []float64

func (sa SortableArray) Len() int {
	return len(sa)
}

func (sa SortableArray) Less(i,j int) bool {
	return sa[i] < sa[j]
}

func (sa SortableArray) Swap(i, j int) {
	sa[i], sa[j] = sa[j], sa[i]
}

func (sb SurfaceBuilder) NewDiscretePosition(p, q int) Vec3 {
	theta := math.Pi / float64(sb.n) * float64(p)
	phi := math.Pi * 2.0 / float64(sb.n) * float64(q)

	dir := Vec3{
		complex(math.Cos(phi)*math.Sin(theta), 0), 
		complex(math.Sin(phi)*math.Sin(theta), 0), 
		complex(math.Cos(theta), 0)}

	chr, _ := NewChristoffelPoly(*sb.mt, dir)

	ar := chr.AllRoots()
	arr := make(SortableArray, len(ar))
	for ind, val := range(ar) {
		arr[ind] = real(val)
	}
	sort.Sort(arr)
	
	val := complex(math.Sqrt(real(sb.rho)/arr[sb.ind]), 0)

	return VectorComplexProduct(dir, 1000*val) 
}

func (sb SurfaceBuilder) formTriangle(dest *os.File, p1, q1, p2, q2, p3, q3 int) {
	v1 := sb.NewDiscretePosition(p1, q1)
	v2 := sb.NewDiscretePosition(p2, q2)
	v3 := sb.NewDiscretePosition(p3, q3)

	cv := 1.0

	fmt.Fprintln(dest, "triangle {")
	fmt.Fprintln(dest, "\t<", real(v1[0]), ",", real(v1[2]), ",", real(v1[1]), ">")
	fmt.Fprintln(dest, "\t<", real(v2[0]), ",", real(v2[2]), ",", real(v2[1]), ">")
	fmt.Fprintln(dest, "\t<", real(v3[0]), ",", real(v3[2]), ",", real(v3[1]), ">")
	fmt.Fprintln(dest, "\tpigment{ color<",cv,",",cv,",",cv,">}")
	fmt.Fprintln(dest, "}")
	fmt.Fprintln(dest)
}

type SurfaceBuilder struct {
	rho complex128
	mt *Tensor4
	n int
	ind int
}

func NewSurfaceBuilder(mt *Tensor4, rho complex128, n, ind int) *SurfaceBuilder {
	return &SurfaceBuilder{rho, mt, n, ind}
}

func BwSurface() {
	mt, rho := NewLinbo3MaterialTensor()

	mtt := Tensor4Rotation(mt, NewChangMatrix(Rad(10)))

	n := 500
	sb := NewSurfaceBuilder(&mtt, rho, n, 0)
	dest, err := os.Create("bw_surface.pov")

	if err != nil {
		log.Panic(err)
	}
	
	for q:=0; q<n; q++ {
		sb.formTriangle(dest, 0, 0, 1, q, 1, q+1)
		sb.formTriangle(dest, n-1, q, n, 0, n-1, q+1)
	}

	for p:=0; p<n; p++ {
		for q:=0; q<n; q++ {
			sb.formTriangle(dest, p, q, p+1, q, p+1, q+1)
			sb.formTriangle(dest, p, q, p+1, q+1, p, q+1)
		}
	}
}
