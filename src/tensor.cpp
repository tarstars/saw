#include "tensor.h"
#include "coeff.h"
#include "mat3.h"
#include "vec3.h"
#include "poly.h"
#include "mat3poly.h"

using namespace std;

Tensor::Tensor()
{
  int p;
  int q;
  int m;
  int n;

  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      for(m = 0; m < 3; ++m)
	for(n = 0; n < 3; ++n)
	  dat[p][q][m][n] = 0;
}

Tensor::Tensor(const Coeff& r)
{
  int p;
  int q;
  int m;
  int n;

  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      for(m = 0; m < 3; ++m)
	for(n = 0; n < 3; ++n)
	  dat[p][q][m][n] = 0;

  dat[0][0][0][0] = dat[1][1][1][1] = r.c11;

  dat[0][0][1][1] = dat[1][1][0][0] = r.c12;

  dat[0][0][2][2] = dat[1][1][2][2] = dat[2][2][0][0] = dat[2][2][1][1] = r.c13;

  dat[2][2][2][2] = r.c33;

  dat[1][2][1][2] = dat[1][2][2][1] = dat[2][1][1][2] = dat[2][1][2][1] = 
    dat[0][2][0][2] = dat[0][2][2][0] = dat[2][0][0][2] = dat[2][0][2][0] = r.c44;

  dat[0][1][1][0] = dat[0][1][0][1] = dat[1][0][0][1] = dat[1][0][1][0] = r.c66;
}

Mat3
Tensor::crist(const Vec3& v)const
{
  Mat3 ret;

  for(int p = 0; p < 3; ++p)
    for(int q = 0; q < 3; ++q)
      for(int r = 0; r < 3; ++r)		
	for(int s = 0; s < 3; ++s)
	  ret.dat[p][q] += dat[p][r][q][s] * v[r] * v[s];

  return ret;
}


std::ostream& 
operator<<(std::ostream& os, const Tensor& ten)
{
  int p, q, r, s;

  for(p = 0; p < 3; ++p)
    {
      for(q = 0; q < 3; ++q)
	{
	  for(r = 0; r < 3; ++r)
	    {
	      for(s = 0; s < 3; ++s)
		{
		  os << ten.dat[p][r][q][s] << " ";
		}
	      os << " | ";
	    }
	  os << endl;
	}
      os << endl;
    }
  return os;
}

Tensor
Tensor::rotate(const Mat3& rm)
{
  Tensor ret;
  int i, j, k, l;
  int p, q, r, s;

  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      for(k = 0; k < 3; ++k)
	for(l = 0; l < 3; ++l)
	  {
	    ret.dat[i][j][k][l] = 0;
	    for(p = 0; p < 3; ++p)
	      for(q = 0; q < 3; ++q)
		for(r = 0; r < 3; ++r)
		  for(s = 0; s < 3; ++s)
		    ret.dat[i][j][k][l] += 
		      dat[p][q][r][s] * 
		      rm.dat[i][p] * 
		      rm.dat[j][q] * 
		      rm.dat[k][r] * 
		      rm.dat[l][s];
	  }
  return ret;
}

Tensor
Tensor::rotate_bidir(const Mat3& rm, const Mat3& im)
{
  Tensor ret;
  int i, j, k, l;
  int p, q, r, s;

  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      for(k = 0; k < 3; ++k)
	for(l = 0; l < 3; ++l)
	  ret.dat[i][j][k][l] = 0;

  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      for(k = 0; k < 3; ++k)
	for(l = 0; l < 3; ++l)
	  for(p = 0; p < 3; ++p)
	    for(q = 0; q < 3; ++q)
	      for(r = 0; r < 3; ++r)
		for(s = 0; s < 3; ++s)
		  ret.dat[i][j][k][l] += dat[p][q][r][s] * rm.dat[p][i] * im.dat[j][q] * rm.dat[k][r] * im.dat[l][s];

  return ret;
}


void
Tensor::test(ostream& os)
{
  int i, j, k, l;
  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      for(k = 0; k < 3; ++k)
	for(l = 0; l < 3; ++l)
	  {
	    if(dat[i][j][k][l] != dat[j][i][k][l])
	      os << i << j << k << l << "!=" << j << i << k << l << endl;
	    if(dat[i][j][k][l] != dat[k][l][i][j])
	      os << i << j << k << l << "!=" << j << i << k << l << endl;
	  }
}

int delta(int p, int q)
{
  return (p == q) ? 1 : 0;
}


Poly 
Tensor::poly_by_n3(double n1, double n2, double g)const
{
  //cerr << "poly by n3" << endl;

  Mat3poly polymatrix;

  CD select[3] = {n1, n2, 1};

  //cerr << "building the polynomial matrix" << endl;

  int p;
  int q;
  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      {
	//cerr << "p q = " << p << " " << q << endl;
	CD coeff[3];
	int i;
	int j;

	for(i = 0; i < 3; ++i)
	  for(j = 0; j < 3; ++j)
	    {
	      //cerr << "\ti j = " << i << " " << j << endl;
	      int ind = delta(i, 2) + delta(j, 2);
	      //cerr << "\t\tind = " << ind << endl;
	      CD ninj = select[i] * select[j];
	      CD tn = dat[p][i][q][j];
	      coeff[ind] += ninj * tn;
	    }

	//cerr << "\t coeff are " << coeff[0] << " " << coeff[1] << " " << coeff[2] << endl;
	int deg = 0;
	int t;

	for(t = 0; t < 3; ++t)
	  if(abs(coeff[t]) > 1e-10)
	    deg = t;
	
	//cerr << "\t deg = " << deg << endl;

	//cerr << "polymatrix " << p << " " << q << " = " << polymatrix[p][q] << endl;
	switch(deg)
	  {
	  case 0: polymatrix[p][q] = Poly(coeff[0]); /*cerr << "\t\tcase 0 " << endl;*/ break;
	  case 1: polymatrix[p][q] = Poly(coeff[0], coeff[1]); /*cerr << "\t\tcase 1 " << endl;*/ break;
	  case 2: polymatrix[p][q] = Poly(coeff[0], coeff[1], coeff[2]); /*cerr << "\t\tcase 2 " << endl;*/ break;
	  }
	//cerr << "polymatrix " << p << " " << q << " = " << polymatrix[p][q] << endl;
      }

  for(int p = 0; p < 3; ++p)
    polymatrix[p][p] = polymatrix[p][p] - Poly(g);

  /*
  cerr << "polynoms are:" << endl;
  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      cerr << p << " " << q << polymatrix[p][q] << endl;
  */    

  Poly ret;
  int prog[6][7] = {{0,0,1,1,2,2,1},{0,1,1,2,2,0,1},{0,2,1,0,2,1,1},{2,0,1,1,0,2,-1},{2,1,1,2,0,0,-1},{2,2,1,0,0,1,-1}};
  for(int t = 0; t < 6; ++t)
    {
      Poly summary(prog[t][6]);
      for(int p = 0; p < 3; ++p)
	{
	  summary = summary * polymatrix[prog[t][2 * p]][prog[t][2 * p + 1]];
	}
      ret = ret + summary;
    }

  //cerr << "end poly by n3" << endl;
  return ret;
}

void 
Tensor::mupad_style_output(ostream &os)const
{
  os << "array(1..3" << endl 
     << "1..3" << endl 
     << "1..3" << endl
     << "1..3" << endl;

  int p, q, r, s;
  for(p = 0; p < 3; ++p)
    for(q = 0; q < 3; ++q)
      for(r = 0; r < 3; ++r)
	for(s = 0; s < 3; ++s)
	  os << "(" << p + 1
	     << "," << q + 1
	     << "," << r + 1
	     << "," << s + 1
	     << ") = " << real(dat[p][q][r][s]) << endl;

  os << "):";
}
