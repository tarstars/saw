package main

import (
	"bsl"
	"fmt"
)

func main() {
	mt, rho := bsl.NewPbTetraborateMaterialTensor()
	n := bsl.NewVec3(1, 0, 0)
	vels := bsl.CalcVelocities(n, bsl.NewUnityMatrix(), mt, rho)
	
	fmt.Println("for the direction ", n)
	for _, val := range(vels) {
		fmt.Println("v=", val.Vel, "q=", val.Q)
	}
}
