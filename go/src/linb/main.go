package main

import (
	"fmt"
	"bsl"
	"log"
	"math"
	//"math/cmplx"
	"os"
)

func drawYZCrossSection(tens bsl.Mat3, flnm string) {
	n := 500.0

	dest, err := os.Create(flnm)
	defer dest.Close()
	if err != nil {
		log.Panic(err)
	}

	for phi := 0.0; phi < 2.0 * math.Pi; phi += 2.0 * math.Pi / n {
		n := bsl.Vec3{
			0, 
			complex(math.Cos(phi), 0), 
			complex(math.Sin(phi), 0)}
		indexes := bsl.IndexByDirection(n, tens)
		for _, val := range(indexes) {
			fmt.Fprintln(dest, real(val.N) * math.Cos(phi), real(val.N) * math.Sin(phi), real(val.N), phi * 180 / math.Pi)
		}
	}
}

func tstEpsRotation() {
	epss := bsl.Linbo3eps()

	epssRot := bsl.Mat3Rotation(epss, bsl.NewRotxMatrix(45.0/180*math.Pi))

	drawYZCrossSection(epss, "rotate_no.txt")
	drawYZCrossSection(epssRot, "rotate_10.txt")
}

func tstProjectionSolution() {
	phi := 13.0 / 180 * math.Pi
	dir := bsl.Normalized(bsl.Vec3{0.0, complex(math.Sin(phi), 0), complex(math.Cos(phi), 0)})
	epss := bsl.Linbo3eps()
	for _, val := range(bsl.IndexByDirection(dir, epss)) {
		indexes := bsl.IndexByProjection(val.N * dir[0], val.N * dir[1], epss)
		for _, pdr := range(indexes) {
			fmt.Println("\tn3 = ", pdr.N3)
			fmt.Println("\tqe = ", pdr.QE)
			fmt.Println("\tqd = ", pdr.QD)
			fmt.Println()
		}
	}

}

func tst3d() {
	dest, err := os.Create("3d_index.inc")
	defer dest.Close()

	if err != nil {
		log.Panic(err)
	}

	epss := bsl.Linbo3eps()
	for x := -3.0; x < 3.0; x += 0.01 {
		for y := -3.0; y < 3.0; y += 0.01 {
			indexes := bsl.IndexByProjection(complex(x, 0), complex(y, 0), epss)
			for _, pdr := range(indexes) {
				z := pdr.N3
				if math.Abs(imag(z)) < 1e-10 {
					fmt.Fprintln(dest, "sphere{< ", x, ", ", real(z), ", ", y, " > 0.01 pigment{color<1, 1, 1>}}")
				}
			}
		}
	}
}

//povDest, err := os.Create("indeces.inc")
//defer povDest.Close()

//if err != nil {
//	log.Panic(err)
//}

/*for m_n := 0; m_n < 3; m_n++*/ 
	//fmt.Println("acoustic result:\n", acousticData)
	
	//theta_o := 0.0 / 180 * math.Pi
	//phi_o := 90.0 / 180 * math.Pi




/*
		n := 60.0
		for z := -1.0; z <= 1; z += 2.0/float64(n) {
			r := math.Sqrt(1 - z * z)
			for phi := 0.0; phi < math.Pi * 2.0; phi += 2.0 * math.Pi / n {
				o_dir := bsl.NewVec3(r * math.Cos(phi), r * math.Sin(phi), z)
				for _, val := range(bsl.IndexByDirection(o_dir, eps)) {
					if real(o_dir[0]) < 0 {
						fmt.Fprintf(povDest, "sphere{<%g, %g, %g> 0.002 pigment{color<0,0,1>}}\n", real(val.N * o_dir[0]), real(val.N * o_dir[2]), real(val.N * o_dir[1]))
					}
				}
			}
		}
*/

func criteria(x, y, z float64) bool {
	return x < 0
}

func tstRotationToZ() {
	theta_a := 13.0 / 180 * math.Pi
	dir_a := bsl.Vec3{complex(0, 0), complex(-math.Sin(theta_a), 0), complex(-math.Cos(theta_a), 0)}
	//transMat := bsl.MatrixByZ(dir_a)
	_ = dir_a
	transMat := bsl.NewRotxMatrix(theta_a)
	epsMainAxis := bsl.Linbo3eps()
	eps := bsl.Mat3Rotation(epsMainAxis, transMat)
	//_ = transMat
	//eps := epsMainAxis

	dest, err := os.Create("index.pov")
	defer dest.Close()

	if (err != nil) {
		log.Panic(err)
	}

	for theta := 0.0; theta < math.Pi; theta += 0.01 {
		for phi := 0.0; phi < math.Pi * 2; phi += 0.01 {
			dir_o := bsl.Vec3{complex(math.Cos(phi) * math.Sin(theta), 0), complex(math.Sin(phi) * math.Sin(theta), 0), complex(math.Cos(theta), 0)}
			for _, val := range(bsl.IndexByDirection(dir_o, eps)) {
				x, y, z := real(dir_o[0] * val.N), real(dir_o[1] * val.N), real(dir_o[2] * val.N)
				if criteria(x, y, z) && math.Abs(imag(val.N)) < 1e-5 {
					fmt.Fprintln(dest, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{color<0, 1, 0>}}")
				}
			}
		}
	}


	if (err != nil) {
		log.Panic(err)
	}

	for x := -3.0; x < 3.0; x += 0.01 {
		for y := -3.0; y < 3.0; y += 0.01 {
			indexes := bsl.IndexByProjection(complex(x, 0), complex(y, 0), eps)
			for _, pdr := range(indexes) {
				z := real(pdr.N3)
				if criteria(x, y, z) && (math.Abs(imag(pdr.N3)) < 1e-10) {
					fmt.Fprintln(dest, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{color<1, 0, 0>}}")
				}
			}
		}
	}
}

					//fmt.Println("source = ", source)
					//fmt.Println("dest = ", dest)

					//e1 := bsl.Vec3Rotation(source.QD, bsl.Transpose(transMat))
					//e2 := bsl.Vec3Rotation(dest.QD, bsl.Transpose(transMat))

					//fmt.Println("e1 = ", e1, 180 * math.Atan2(real(e1[2]), real(e1[1])) / math.Pi)
					//fmt.Println("e2 = ", e1, 180 * math.Atan2(real(e2[2]), real(e2[1])) / math.Pi)
					
					//sx, sz, sy := real(nx), real(source.N3), real(ny)
					/*			fmt.Fprintf(povDest, "sphere{<%g, %g, %g> 0.015 pigment{color<0,1,0>}}\n", sx, sz, sy)
			fmt.Fprintf(povDest, "cylinder{<%g, %g, %g> <%g, %g, %g> 0.005 pigment{color<1,0,1>}}\n", 
				sx - real( source.QD[0]) * 0.125, 
				sz - real( source.QD[2]) * 0.125, 
				sy - real( source.QD[1]) * 0.125,
				sx + real( source.QD[0]) * 0.125, 
				sz + real( source.QD[2]) * 0.125, 
				sy + real( source.QD[1]) * 0.125)


			dx, dz, dy := real(nx), real(dest.N3), real(ny)
			fmt.Fprintf(povDest, "sphere{<%g, %g, %g> 0.015 pigment{color<0,0,1>}}\n", dx, dz, dy)
			fmt.Fprintf(povDest, "cylinder{<%g, %g, %g> <%g, %g, %g> 0.005 pigment{color<1,0,1>}}\n", 
				dx - real( dest.QD[0]) * 0.125, 
				dz - real( dest.QD[2]) * 0.125, 
				dy - real( dest.QD[1]) * 0.125,
				dx + real( dest.QD[0]) * 0.125, 
				dz + real( dest.QD[2]) * 0.125, 
				dy + real( dest.QD[1]) * 0.125)


			fmt.Println("source angle = ", math.Atan2(real(source.N3), math.Sqrt(sx * sx + sy * sy)) / math.Pi * 180)
			fmt.Println("dest angle = ", math.Atan2(real(dest.N3), math.Sqrt(sx * sx + sy * sy)) / math.Pi * 180)
*/			

	//destIndProj, err := os.Create("projection_index.pov")
	//defer destIndProj.Close()



func testMatrix() {
	theta_a := 30.0 / 180 * math.Pi
	dir_a := bsl.Vec3{complex(0, 0), complex(math.Sin(theta_a), 0), complex(math.Cos(theta_a), 0)}
	
	fmt.Println("matrix:\n", bsl.NewRotxMatrix(-theta_a))

	fmt.Println("matrix by z:\n", bsl.MatrixByZ(dir_a))

}

func testAcousticProblem() {
	theta_a := -13.0 * math.Pi / 180
	//theta_a := 0.0
	/*fmt.Println(
		bsl.CalcPiezoVelocities(bsl.NewVec3(0, math.Sin(theta_a), math.Cos(theta_a)), 
		bsl.NewRotxMatrix(0)))
	 */
	_ = theta_a

	dest, err := os.Create("linb_cross.pov")
	defer dest.Close()

	if err != nil {
		log.Panic(err)
	}

	m := 5.0
	newCoordSys := bsl.NewRotxMatrix(theta_a)
	fmt.Println("transformation matrix:")
	bsl.MatrixPPrint(newCoordSys, "\t")

	for phi := 0.0; phi < 2 * math.Pi; phi += 0.01 {
		dir := bsl.NewVec3(0, math.Cos(phi), math.Sin(phi))
		vels, _ := bsl.CalcPiezoVelocities(dir, newCoordSys)
		for _, val := range(vels) {
			s := 1000.0 / val.Vel * m
			x, y, z := real(dir[0]) * s , real(dir[1]) * s, real(dir[2]) * s
			fmt.Fprintln(dest, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{color<1, 0, 0>}}")
		}
	}

	
}


func testRotation() {
	theta_a := 13.0 / 180 * math.Pi
	fmt.Println("transformation into our system:")
	bsl.MatrixPPrint(bsl.NewRotxMatrix(-theta_a), "\t")

	/*
	dir_a := bsl.Vec3{complex(0, 0), complex(math.Sin(theta_a), 0), complex(math.Cos(theta_a), 0)}
	transMat := bsl.MatrixByZ(dir_a)

	a := bsl.NewVec3(0, 1, 0)
	fmt.Println("a = ", a)
	b := bsl.Vec3Rotation(a, transMat)
	fmt.Println("b = ", b) 
*/
}

func fitSlownessData() {
	theta_a := -14.0 * math.Pi / 180
	//theta_a := 0.0
	/*fmt.Println(
		bsl.CalcPiezoVelocities(bsl.NewVec3(0, math.Sin(theta_a), math.Cos(theta_a)), 
		bsl.NewRotxMatrix(0)))
	 */
	_ = theta_a

	dest, err := os.Create("linb_cross_e.pov")
	defer dest.Close()

	if err != nil {
		log.Panic(err)
	}

	destAngle, err := os.Create("linb_e_angle.txt")
	defer dest.Close()

	if err != nil {
		log.Panic(err)
	}	

	m := 10.0
	newCoordSys := bsl.NewRotxMatrix(theta_a)

	for phi := 0.0; phi < 2 * math.Pi; phi += 0.003 {
		dir := bsl.NewVec3(0, math.Cos(phi), math.Sin(phi))
		vels, _ := bsl.CalcPiezoVelocities(dir, newCoordSys)
		for _, val := range(vels) {
			s := 1000.0 / val.Vel * m
			x, y, z := real(dir[0]) * s , real(dir[1]) * s, real(dir[2]) * s
			fmt.Fprintln(dest, "sphere{< ", x, ", ", z, ", ", y, " > 0.01 pigment{color<1, 0, 0>}}")
		}
	}

	source, err := os.Open("extra.txt")
	if err != nil {
		log.Panic(err)
	}

	var x0, y0, Mx, My, x, y, intensity float64
	fmt.Fscanf(source, "%g %g %g %g", &x0, &y0, &Mx, &My)
	fmt.Fscanf(source, "\n")

	for ;; {
		_, err := fmt.Fscanf(source, "%g %g %g", &x, &y, &intensity)
		if err != nil {
			fmt.Print("error of file reading: ", err)
			break
		}

		L := 1.72e3
		lambda := 532e-9
		f := 88e6
		scale := L * lambda * f
		
		x3d :=  (x-x0) / (Mx * scale)
		y3d := -(y-y0) / (My * scale)

		fmt.Fprintln(dest, "sphere{< ", 0, ", ", y3d * 1000 * m, ", ", x3d * 1000 * m, " > ", 0.02 * math.Log10(intensity) , " pigment{color<0, 1, 0>}}")
	
		angle := (math.Atan2(y3d, x3d) + theta_a) * 180 / math.Pi
		fmt.Fprintln(destAngle, angle, intensity)
	}
}

func main() {
	//fitSlownessData()
	//testRotation()
	//testAcousticProblem()
	//bsl.TestLinbPolarization()
	//testMatrix()
	//tstRotationToZ()
	//tstEpsRotation()
	//tstProjectionSolution()
	//tst3d()
	//bsl.TestLaguerre()
	//bsl.TestFresnel()
	//bsl.TestMeritRotation()
	bsl.DispatchCalcFigureOfMerit()
	//tstAngles()
	//bsl.MixedCase()
	//bsl.CalcParatellurite()
	//bsl.TestTeO2()
}
