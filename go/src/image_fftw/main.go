package main

import ("image/color"
	"fmt"
	"github.com/runningwild/go-fftw"
	"image/png"
	"log"
	"math/cmplx"
	"os"
	)

func loadPictureSaveSpectrum() {
	source, err := os.Open("circle_100_100.png")
	if err != nil {
		log.Fatal(err)
	}
	defer source.Close()

	smile, err := png.Decode(source)
	if err != nil {
		log.Fatal(err)
	}

	rc := smile.Bounds()
	minx, miny := rc.Min.X, rc.Min.Y
	maxx, maxy := rc.Max.X, rc.Max.Y

	for p := miny; p < maxy; p++ {
		for q:= minx; q < maxx; q ++ {
			fmt.Print(color.GrayModel.Convert(smile.At(q, p)).(color.Gray).Y/255)
		}
		fmt.Println()
	}
	
	nx, ny := (maxx - minx), (maxy - miny)
	datSource := fftw.Alloc2d(ny, nx)
	datDest := fftw.Alloc2d(ny, nx)

	plan := fftw.PlanDft2d(datSource, datDest, fftw.Forward, fftw.Estimate)

	for p := miny; p < maxy; p++ {
		for q:= minx; q < maxx; q ++ {
			datSource[p][q] = complex(float64(color.GrayModel.Convert(smile.At(q, p)).(color.Gray).Y)/255.0, 0)
		}
	}

	plan.Execute()


	destFour, err := os.Create("four.txt")
	if err != nil {
		log.Panic(err)
	}

	for p := miny; p < maxy; p++ {
		for q:= minx; q < maxx; q ++ {
			fmt.Fprintln(destFour, p, q, cmplx.Abs(datDest[p][q]))
		}
	}
}



func main() {
     	
}
