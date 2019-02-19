#version 3.7;

global_settings {assumed_gamma 1.0}

#include "colors.inc"

#declare viewAngle = clock * 2 * pi;
#declare cameraDist = 10;     
#declare customColor = <1, 1, 0>;
#declare csSize = 3;

#switch(clock)
  #range (0,1)
    camera{
      location < cos(viewAngle) * cameraDist, 
      2.5, 
      sin(viewAngle) * cameraDist >
      
      angle 50  
      look_at  <0, 1, 0> 
    }             
  #break
  
  #range(1,2)
    camera{
      location < cos(viewAngle) * cameraDist, 
      2.5, 
      sin(viewAngle) * cameraDist >
      
      angle 50  
      look_at  <0, 1, 0> 
    }             
  #break
#end

plane{
   y, -3
   pigment{White}
}     

light_source {
    <cos(viewAngle - pi / 4) * cameraDist, 10, sin(viewAngle - pi / 4) * cameraDist>
    color<1, 1, 1>
}

light_source {
    <cos(viewAngle + pi / 4) * cameraDist, 10, sin(viewAngle + pi / 4) * cameraDist>
    color<1, 1, 1>
}

#include "coordinate_system.inc"

intersection {
  union {
    #include "3d_index.inc"
  }
  plane {
    <1, 1, 2>,
    1.0
  }
}
