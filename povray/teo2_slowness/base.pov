#version 3.7;

global_settings {assumed_gamma 1.0}

#include "colors.inc"

#declare viewAngle = clock * 2 * pi;
#declare cameraDist = 10;     
#declare customColor = <1, 1, 0>;
#declare csSize = 2;

camera{
  location < cos(viewAngle) * cameraDist, 
             2.5, 
             sin(viewAngle) * cameraDist >

  angle 40  
  look_at  <0, 0, 0> 
}             

plane{
   y, -3
   pigment{White}
}     

light_source {
    <cos(viewAngle) * cameraDist, 10, sin(viewAngle) * cameraDist>
    color<1, 1, 1>
}

#include "coordinate_system.inc"
#include "scene.pov"
