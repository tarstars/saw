#ifndef _COEFF_
#define _COEFF_

class Coeff
{
 public:
  double c11, c12, c13, c33, c44, c66, c16;

 Coeff(double d11, double d12, double d13, double d33, double d44, double d66, double d16):
  c11(d11), c12(d12), c13(d13), c33(d33), c44(d44), c66(d66), c16(d16){}

  static Coeff 
    paratellurite()
  {
    return Coeff(5.6, 5.1, 2.2, 10.6, 2.65, 6.6, 0);
  }
};

#endif
