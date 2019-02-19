package bsl


//Make material matrix for paratellurite crystal
func NewParatelluriteMaterialTensor() (Tensor4, complex128) {
	return NewTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10), 5.96e3
}

//Make material matrix for linbo3 crystal
func NewLinbo3MaterialTensor() (Tensor4, complex128) {
	return NewTrigonalMaterialTensor(19.886e10, 5.467e10, 6.799e10, 0.783e10, 23.418e10, 5.985e10, 7.209e10), /*4642.8*/4700
}

//Make material matrix for tellurium
func NewTeMaterialTensor() (Tensor4, complex128) {
	//from 1979_Royer
	c11:=complex(3.257e10,0)
	c12:=complex(0.845e10,0)
	c13:=complex(2.57e10,0)
	c14:=complex(1.238e10,0)
	c33:=complex(7.17e10,0)
	c44:=complex(3.094e10,0)
	c66:=complex(1.206e10,0)
	return NewTrigonalMaterialTensor(c11, c12, c13, c14, c33, c44, c66), 6250
}

//epsilon for tellurium crystal
func NewTeEps() Mat3 {
	return Uniaxial(4.8*4.8, 6.2*6.2)
}

//photoelastic tensor for lithium niobate
func NewTePhotoelastic() Tensor4 {
	return NewTrigonalPhotoelasticTensor(0.164, 0.138, 0.146, -0.086, 0.038, -0.04, 0.28, 0.14)
}

//piezo tensor for lithium niobate 
func NewLinbo3PiezoTensor() Tensor3 {
	return NewTrigonalPiezoTensor(3.655, 2.407, 0.328, 1.894)
}

//photoelastic tensor for lithium niobate
func NewLinbo3Photoelastic() Tensor4 {
	return NewTrigonalPhotoelasticTensor(-0.021, 0.06, 0.172, 0.141, 0.118, -0.052, 0.109, 0.121)
}

//photoelastic tensor for paratellurite
func NewTeo2Photoelastic() Tensor4 {
	return NewTetragonalPhotoelasticTensor(0.0074, 0.187, 0.34, 0.0905, 0.24, -0.17, -0.463)
}

//static epsilon constants for lithium niobate
func Linbo3StaticEps() Mat3 {
	eps0 := complex(8.8542e-12, 0)
	exx, ezz := eps0 * 44.9, eps0 * 26.7
	return Uniaxial(exx, ezz)
}

//tensor epsilon for LiNbO3 crystals
func Linbo3eps() Mat3 {
	return Uniaxial(/*5.3973, 4.9890*/ 5.3824, 4.9729)
}

//tensor epsilon for TeO2 crystals
func Teo2eps() Mat3 {
	/*no:=complex(2.1895, 0)
	ne:=complex(2.3309, 0)*/

	return Uniaxial(5.1984, 6.0025)
}

//tensor epsilon for uniaxial crystals
func Uniaxial(exx, ezz complex128) Mat3 {
	return Mat3{
		{exx, 0, 0}, 
		{0, exx, 0}, 
		{0, 0, ezz}}
}

//tensor epsilon for uniaxial crystals
func Biaxial(exx, eyy, ezz complex128) Mat3 {
	return Mat3{
		{exx, 0, 0}, 
		{0, eyy, 0}, 
		{0, 0, ezz}}
}

//lead tetraborate data according to I. Martynyuk-Lototska et al.
//Optical Materials 31 (2009) 660-667

func NewSrTetraborateMaterialTensor() (Tensor4, complex128) {
	c11:=complex(304e9, 0)
	c12:=complex(70e9, 0)
	c13:=complex(49e9, 0)
	c22:=complex(268e9, 0)
	c23:=complex(55e9, 0)
	c33:=complex(378e9, 0)
	c44:=complex(139e9, 0)
	c55:=complex(120e9, 0)
	c66:=complex(133e9, 0)

	return NewOrthorombicMaterialTensor(c11, c12, c13, c22, c23, c33, c44, c55, c66), 4011
}

func NewSrTetraboratePhotoelasticTensor() Tensor4 {
	p11 := complex(-0.9, 0)
	p12 := complex(-0.6, 0)
	p13 := complex(0.21, 0)

	p21 := complex(-0.53, 0)
	p22 := complex(-0.4, 0)
	p23 := complex(0.14, 0)

	p31 := complex(-0.9, 0)
	p32 := complex(-0.6, 0)
	p33 := complex( 0.4, 0)

	p44 := complex(0.5, 0)
	p55 := complex(0.5, 0)
	p66 := complex(0.5, 0)

	return NewOrthorombicPhotoelasticTensor(p11, p12, p13,    p21, p22, p23,    p31, p32, p33,   p44, p55, p66)
}

func NewSrTetraborateEps() Mat3 {
	nxx := complex(1.7333, 0)
	nyy := complex(1.7323, 0)
	nzz := complex(1.7356, 0)
	return Biaxial(nxx*nxx, nyy*nyy, nzz*nzz)
}


//pb constants
func NewPbTetraborateMaterialTensor() (Tensor4, complex128) {
	c11:=complex(342e9, 0)
	c12:=complex(9.6e9, 0)
	c13:=complex(67e9, 0)
	c22:=complex(304e9, 0)
	c23:=complex(64e9, 0)
	c33:=complex(367e9, 0)
	c44:=complex(135.4e9, 0)
	c55:=complex(117.5e9, 0)
	c66:=complex(129.8e9, 0)

	return NewOrthorombicMaterialTensor(c11, c12, c13, c22, c23, c33, c44, c55, c66), 5852
}

func NewPbTetraboratePhotoelasticTensor() Tensor4 {
	p11 := complex(-0.06, 0)
	p12 := complex(0.4, 0)
	p13 := complex(0.8, 0)

	p21 := complex(-0.1, 0)
	p22 := complex(-0.3, 0)
	p23 := complex(0.9, 0)

	p31 := complex(0.4, 0)
	p32 := complex(-0.6, 0)
	p33 := complex(0.8, 0)

	p44 := complex(-0.4, 0)
	p55 := complex(0.17, 0)
	p66 := complex(0.5, 0)

	return NewOrthorombicPhotoelasticTensor(p11, p12, p13,    p21, p22, p23,    p31, p32, p33,  p44, p55, p66)
}

func NewPbTetraborateEps() Mat3 {
	nxx := complex(1.9326, 0)
	nyy := complex(1.9275, 0)
	nzz := complex(1.9309, 0)
	return Biaxial(nxx*nxx, nyy*nyy, nzz*nzz)
}



