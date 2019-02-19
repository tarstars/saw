#version 3.7;

global_settings {
  assumed_gamma  1.0
}

#declare csSize = 3;
#declare planeCameraDist = 15;
#declare rotationClock = 0.03;

camera{
  location <1 + planeCameraDist * cos(2 * pi * rotationClock), csSize, 1 + planeCameraDist * sin(2 * pi * rotationClock)>
  look_at <1,1,1>
  angle 30
}

light_source{
  <5,5,10>
  color<1,1,1>
}

cylinder{
  <-0.3, 0, 0> * csSize
  < 1, 0, 0> * csSize
  0.05
  pigment{
    color<1,0,0>
  }
}

text {
  ttf "timrom.ttf" "X" 0.05, 0
  scale 0.7            
  rotate <0,50,0>
  translate <csSize, 0.1, 0>
  pigment {
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

text {
  ttf "timrom.ttf" "Y" 0.05, 0
  scale 0.7
  rotate <0, 90, 0>
  translate <0, 0.1, csSize>
  pigment {
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

text {
  ttf "timrom.ttf" "Z" 0.05, 0
  scale 0.7
  rotate <0, -150, 0>
  translate <0.1, csSize, 0>
  pigment {
    color<0, 0, 1>
  }
}


#declare initialPosition = <0.1, 0.1, 0.1>;

#declare shiftx = <0.1, 0.0, 0.0>;
#declare shifty = <0.0, 0.0, 0.1>;
#declare shiftz = <0.0, 0.1, 0.0>;

#declare mv = 0.5;

#declare sxx = 0;
#declare syy = 0;
#declare szz = 0;
#declare syx = 0;
#declare sxy = 0;
#declare szx = 0;
#declare sxz = 0; 
#declare szy = 0;
#declare syz = 0;

#switch(clock)
  #range (0, 1)
    #declare sxx = sin(2 * pi * (clock - 0)) * mv;
    text {
      ttf "timrom.ttf" "Sxx" 0.05, 0
      rotate <0, -90, 0>
      translate <csSize / 2, csSize * 0.8, csSize / 3>
      pigment {
        color<1, 0, 0>
      }
    }
  #break
  #range (1, 2)
    #declare syy = sin(2 * pi * (clock - 1)) * mv;
    text {
      ttf "timrom.ttf" "Syy" 0.05, 0
      rotate <0, -90, 0>
      translate <csSize / 2, csSize * 0.8, csSize / 3>
      pigment {
        color<0, 1, 0>
      }
    }
  #break
  #range (2, 3)
    #declare szz = sin(2 * pi * (clock - 2)) * mv;
    text {
      ttf "timrom.ttf" "Szz" 0.05, 0
      rotate <0, -90, 0>
      translate <csSize / 2, csSize * 0.8, csSize / 3>
      pigment {
        color<0, 0, 1>
      }
    }
  #break
  #range (3, 4)
    #declare syz = 0 * sin(2 * pi * (clock - 3)) * mv;
    #declare szy = sin(2 * pi * (clock - 3)) * mv;
    text {
      ttf "timrom.ttf" "Syz" 0.05, 0
      rotate <0, -90, 0>
      translate <csSize / 2, csSize * 0.8, csSize / 3>
      pigment {
        color<0, 1, 1>
      }
    }
  #break
  #range (4, 5)
    #declare sxz = 0 * sin(2 * pi * (clock - 4)) * mv;
    #declare szx = sin(2 * pi * (clock - 4)) * mv;
   text {
      ttf "timrom.ttf" "Sxz" 0.05, 0
      rotate <0, -90, 0>
      translate <csSize / 2, csSize * 0.8, csSize / 3>
      pigment {
        color<1, 0, 1>
      }
   }
  #break
  #range (5, 6)
    #declare sxy = 1 * sin(2 * pi * (clock - 5)) * mv;
    #declare syx = 0 * sin(2 * pi * (clock - 5)) * mv;
    text {
      ttf "timrom.ttf" "Sxy" 0.05, 0
      rotate <0, -90, 0>
      translate <csSize / 2, csSize * 0.8, csSize / 3>
      pigment {
        color<1, 1, 0>
      }
    }
  #break
#end


#for(indx, 0, 20) 
  #for(indy, 0, 20)
    #for(indz, 0, 6)
      #declare position = initialPosition + shiftx * indx + shifty * indy + shiftz * indz;
      #declare position = position + (indx * sxx + indy * syx + indz * szx) * shiftx;
      #declare position = position + (indx * sxy + indy * syy + indz * szy) * shifty;
      #declare position = position + (indx * sxz + indy * syz + indz * szz) * shiftz;

      sphere{
        position,
        0.05
        pigment {
          color<1,1,0>
        }
      }
    #end
  #end
#end

