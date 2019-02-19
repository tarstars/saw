package main

import (
	"bsl"
	//"fmt"
)


func BeamStructCalculation() {
	f := complex(30e6, 0)
	Ax := complex(1e-2, 0)
	forceProfile := "fp_500x300_center.png"
	nz := 300
	forceDir := bsl.Vec3{1, 0, 0}
	tens, rho := bsl.NewParatelluriteMaterialTensor()
	tens = bsl.Tensor4Rotation(tens, bsl.NewChangMatrix(bsl.Rad(10)))

	bsl.CalcFieldStructure(f, Ax, nz, tens, rho, forceProfile, forceDir, "teo2_x_30mhz")	
}

func main() {
	//bsl.TestLoadPicture()
	BeamStructCalculation()
	//bsl.TestSaveMatrix()
	//bsl.TestConvolution()
	//bsl.TestCopy()
	//bsl.TestLayeredTransform()
}
