package main

import (
	"log"
	"fmt"
	"os"
	"bsl"
	"math"
	"math/cmplx"
)

func drawGraph() {
	//a := 0.0000111
	//b := 0.0011446
	a,b := 0.0, 0.002

	alpha := 45.0 / 180 * 3.14159265358

	n := 10000.0
	d := (b - a) / n

	dest, err := os.Create("teo2_5.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer dest.Close()

	mt, rho := bsl.NewParatelluriteMaterialTensor()

	for t := a; t <= b; t += d {
		cw := bsl.NewCompositeWave(mt, rho, complex(t * math.Cos(alpha), 0), complex(t * math.Sin(alpha), 0), 1e6, bsl.NewVec3(1, 1, 1))
		fmt.Fprintln(dest, t, cmplx.Abs(cw.GetFarnDet()))
	}
}

func debugOnePoint() {
	mt, rho := bsl.NewParatelluriteMaterialTensor()

	t := 0.0004/math.Sqrt(2)
	cw := bsl.NewCompositeWave(mt, rho, complex(t, 0), complex(t, 0), 1, bsl.NewVec3(1, 1, 1))
	
	cw.Print(os.Stdout)
}

func main() {
	debugOnePoint()
	//drawGraph()
}
