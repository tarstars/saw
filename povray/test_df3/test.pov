#version 3.7;

global_settings {assumed_gamma 1.0}

#declare csSize = 3;
#declare cameraDist = 13;

camera {
  location <cameraDist * cos(2 * pi * clock), 2, cameraDist * sin(2 * pi * clock)>
  look_at <0, 1, 0>
  
  angle 30
}

cylinder{
  <-0.3, 0, 0> * csSize
  < 1, 0, 0> * csSize
  0.05
  pigment{
    color<1,0,0>
  }
}

cylinder{
  <0, 0, -0.3> * csSize
  <0, 0,  1.0> * csSize
  0.05
  pigment{
    color<0,1,0>
  }
}

cylinder{
  < 0, -0.3, 0> * csSize
  < 0, 1, 0> * csSize
  0.05
  pigment{
    color<0,0,1>
  }
}

light_source {
  <1, 0.5, 20>
  color <1.1, 1.1, 1.1>
}

#declare box_interior = interior {
  media {
    //emission <1,1,1> * 1
    absorption <1,1,1> * 20
    //scattering { 1, <1.0, 1.0, 1.0> * 0.0001 }
    density {
      density_file df3 "sphere.df3"
      color_map {
	[0.00 rgb <0.0, 0.0, 0.0>]
	[0.25 rgb <0.0, 0.0, 0.0>]
	[0.25 rgb <1.0, 0.0, 0.0> * 0.1]
	[0.50 rgb <1.0, 0.0, 0.0> * 0.1]
	[0.50 rgb <0.0, 1.0, 0.0> * 0.5]
	[0.75 rgb <0.0, 1.0, 0.0> * 0.5]
	[0.75 rgb <0.0, 0.0, 1.0> * 1.0]
	[1.00 rgb <0.0, 0.0, 1.0> * 1.0]
//	[0.4 rgb <0,1.0,0>]
//	[0.6 rgb <0,0,1.0>]
	// [0.0 rgb <1,1,1> * 0.0]
	// [0.1 rgb <1,1,1> * 0.1]
	// [0.2 rgb <1,1,1> * 0.2]
	// [0.3 rgb <1,1,1> * 0.3]
	// [0.4 rgb <1,1,1> * 0.4]
	// [0.5 rgb <1,1,1> * 0.5]
	// [0.6 rgb <1,1,1> * 0.6]
	// [0.7 rgb <1,1,1> * 0.7]
	// [0.8 rgb <1,1,1> * 0.8]
	// [0.9 rgb <1,1,1> * 0.9]
	// [1.0 rgb <1,1,1> * 1.0]
	
	//[0.2 rgb <.01,0,0>]
	//[1.0 rgb <0,0,1>]
	//[0.6 rgb <0,1,0>]
	//[1.0 rgb <1,0,0>]
	//[1.0 rgb <1,1,1>]
      }
    }
  }
}

box {
  <0.0, 0.0, 0.0>
  <1.0, 1.0, 1.0>
  pigment {rgbf 1.0}
  interior {box_interior}
  hollow
  translate<1.0, 0, 1.0>
}

box {
  <0.0, 0.0, 0.0>
  <1.0, 1.0, 1.0>
  pigment {rgbf 1.0}
  interior {box_interior}
  hollow
  translate <-1.0, 0, 1.0>
}

plane {
  y, -1
  pigment {
    color <1,1,1>
  }
}

plane {
  z, -3
  pigment {
    color <1,1,1>
  }
}
