cylinder {
  <-0.3, 0, 0> * csSize
  < 1.0, 0, 0> * csSize
  0.01 * csSize
  pigment{color <1, 0, 0> }
}                        

cone {
  < 1.0, 0, 0> * csSize, 0.05 * csSize, 
  < 1.1, 0, 0> * csSize, 0  
  pigment{color <1, 0, 0> }
}

cylinder {
  <0, 0, -0.3> * csSize
  <0, 0,  1.0> * csSize
  0.01 * csSize
  pigment{color <0, 1, 0>}
} 

cone {
  < 0.0, 0, 1.0> * csSize, 0.05 * csSize, 
  < 0.0, 0, 1.1> * csSize, 0  
  pigment{color <0, 1, 0> }
}


cylinder {
  <0, -0.3, 0> * csSize
  <0,  1.0, 0> * csSize
  0.01 * csSize
  pigment{color <0, 0, 1>}
}                   

cone {
  < 0.0, 1.0, 0> * csSize, 0.05 * csSize, 
  < 0.0, 1.1, 0> * csSize, 0  
  pigment{color <0, 0, 1> }
}


text{
        ttf "timrom.ttf" "X"
        0.1 0
        pigment{ color <1, 0, 0> }
        
        scale <1, 1, 1> * csSize / 4
        translate <1, 0.1, 0> * csSize        
}

text{
        ttf "timrom.ttf" "Y"
        0.1 0
        pigment{ color <0, 1, 0> }
        
        scale <1, 1, 1> * csSize / 4
        translate <0, 0.1, 1> * csSize        
}

text{
        ttf "timrom.ttf" "Z"
        0.1 0
        pigment{ color <0, 0, 1> }

        rotate <0, 180, 0>        
        scale <1, 1, 1> * csSize / 4
        translate <0.25, 1, 0> * csSize        
}
