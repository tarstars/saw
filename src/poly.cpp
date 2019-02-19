#include <iostream>
#include <cmath>
#include "poly.h"

using namespace std;

ostream&
operator<<(ostream &os, const Poly &r)
{
  for(int t = 0; t < (int)r.dat.size(); ++t)
    os << r.dat[t] << " ";
  return os;
}

template<int s>
Poly ar2pol(double ar[s])
{
  Poly ret;
  ret.dat.resize(s);
  for(int t = 0; t < s; ++t)
    ret.dat[t] = ar[t];
  return ret;
}

Poly
Poly::operator/(const Poly& r)const
{
  poly_type a = dat, b = r.dat;
  int n = a.size() - 1;
  int m = b.size() - 1;

  Poly ret;
  ret.dat = poly_type(n - m + 1);


  for(int t = n - m; t >= 0; --t)
    {
      poly_type::value_type k = a[t + m] / b[m];
      ret.dat[t] = k;
      for(int p = 0; p < m; ++p)
	a[t + p] -= k * b[p];
    }

  return ret;

}

Poly
Poly::operator*(const Poly& r)const
{
  Poly ret;
  const poly_type &a = dat;
  const poly_type &b = r.dat;

  int n = a.size() - 1;
  int m = b.size() - 1;

  ret.dat = poly_type(n + m + 1);

  for(int t = 0; t <= m; ++t)
    for(int p = 0; p <= n; ++p)
	ret.dat[p + t] += a[p] * b[t];

  return ret;
}

Poly
Poly::operator+(const Poly& r) const
{
  //cerr << "operator + " << endl;
  Poly ret;
  int n = dat.size();
  int m = r.dat.size();

  ret.dat = poly_type(max(m, n));
  for(int t = 0; t < min(m,n); ++t)
    ret.dat[t] = dat[t] + r.dat[t];

  for(int t = min(m,n); t < max(m,n); ++t)
    if(n > m)
      ret.dat[t] = dat[t];
    else
      ret.dat[t] = r.dat[t];

  return ret;
  
  }

Poly
Poly::operator-(const Poly& r) const
{
  Poly ret;
  int n = dat.size();
  int m = r.dat.size();

  ret.dat = poly_type(max(m, n));
  for(int t = 0; t < min(m,n); ++t)
    ret.dat[t] = dat[t] - r.dat[t];

  for(int t = min(m,n); t < max(m,n); ++t)
    if(n > m)
      ret.dat[t] = dat[t];
    else
      ret.dat[t] = - r.dat[t];

  return ret;
  
}

template<typename T>
void
Poly::fst(T x, T *pv, T *pdv, T *pddv)const
{
  T f = 0;
  T s = 0;
  T u = 0;

  for(int t = dat.size() - 1; t >= 0; --t)
    {
      u = u * x + s;
      s = s * x + f;
      f = f * x + dat[t];
    }

  if(pv) *pv = f;
  if(pdv) *pdv = s;
  if(pddv) *pddv = u;
}

poly_type::value_type 
Poly::operator()(poly_type::value_type x)const
{
  poly_type::value_type ret = 0;
  for(int t = dat.size() - 1; t >= 0; --t)
    ret = ret * x + dat[t];
  return ret;
}

//#define _TEST1_

template<typename T>
T laguerr(const Poly& r)
{
#ifdef _TEST1_
  cerr << "begin laguerr " << endl;
#endif

  T x = 1e-19;
  T a = 1; 

  for(int k = 0; k < 1000 && !isnan(abs(a)) && abs(a) > 1e-14; ++k)
    {
#ifdef _TEST1_
      cerr << "x = " << x << endl;
#endif

      T v;
      T dv;
      T ddv;

      r.fst(x, &v, &dv, &ddv);
      
      int n = r.dat.size() - 1;
      T g = dv / v;
      T h = g * g - ddv / v;
      T ddenum = sqrt(complex<long double>(n - 1.0) * (h * complex<long double>(n) - g * g));
      T denum1 = g + ddenum;
      T denum2 = g - ddenum;
      T denum = 0;
      if(abs(denum1) > abs(denum2))
	denum = denum1;
      else
	denum = denum2;

      a = complex<long double>(n) / denum;

#ifdef _TEST1_
      cerr << "a = " << a << " " << isnan(abs(a)) << endl;
#endif

      if(!isnan(abs(a)))
	x -= a;
    }

#ifdef _TEST1_
  cerr << "end laguerr" << endl;
#endif

  return x;
}

poly_type
all_roots(const Poly& r)
{
#ifdef _TEST_
  cerr << "all roots" << endl;
  cerr << "start from " << r << endl;
#endif

  poly_type ret;
  Poly lt(r);
  while(lt.deg() > 0)
    {
      complex<long double> root = laguerr<complex<long double> >(lt);
#ifdef _TEST_ 
      cerr << "\troot is " << root << endl;
#endif

      ret.push_back(root);
      lt = lt / Poly(-root,1);
#ifdef _TEST_
      cerr << "\tpoly is " << lt << endl;
#endif

    }
  return ret;
}

#ifdef _TEST_
int main()
{
  double a[]={1,3,3,1}, b[] = {1,2,1}, c[] = {1,3,3,1};
  double d[] = {1,4,6,4,1};
  cout << ar2pol<4>(a) << endl;
  cout << ar2pol<3>(b) << endl;
  cout << ar2pol<5>(d) / ar2pol<3>(b) << endl;
  cout << ar2pol<4>(a) / ar2pol<3>(b) << endl;
  cout << ar2pol<4>(a) / ar2pol<4>(c) << endl;

  double e[] = {1,1}, f[] = {1};
  Poly pa = ar2pol<2>(e);
  Poly pb = ar2pol<1>(f);

  cout << endl;

  for(int t = 0; t < 10; ++t)
    {
      cout << pb << endl;
      pb = pb * pa;
    }

  cout << endl;

  cout << ar2pol<4>(a) << endl;
  cout << ar2pol<3>(b) << endl;
  cout << ar2pol<4>(a) * ar2pol<3>(b) << endl;

  double tlar[] = {1, 0, 0, 0, 0, 1};
  Poly lt = ar2pol<6>(tlar);
  cerr << "laguerr test = " << laguerr<complex<double> >(lt) << endl;

  cerr << "full laguerr: " << endl;
  poly_type roots = all_roots(lt);
  for(int t = 0; t < roots.size(); ++t)
    cerr << roots[t] << " ";
  cerr << endl;

}
#endif
