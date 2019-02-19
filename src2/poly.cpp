/// @file poly.cpp
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#include "poly.h"

#include <gsl/gsl_poly.h>

#include <ostream>
#include <cmath>

#include <boost/assign/list_of.hpp>

using namespace boost;
using namespace boost::assign;
using namespace std;

namespace farn {

    ostream& operator<<(ostream &os, Poly const& p) {
        for (Poly::CoeffVec::const_iterator i = p.c.begin(); i != p.c.end(); ++i) {
            os << *i << " ";
        }
        return os;
    }

    ostream& operator<<(ostream& os, Poly::RootVec const& p) {
        for (Poly::RootVec::const_iterator i = p.begin(); i != p.end(); ++i) {
            os << *i << " ";
        }
        return os;
    }

    Poly Poly::operator/(const Poly& r) const {
        CoeffVec a = c;
        CoeffVec const& b = r.c;
        int n = a.size() - 1;
        int m = b.size() - 1;

        Poly ret;
        ret.c = CoeffVec(n - m + 1);

        for(int t = n - m; t >= 0; --t)
        {
            CoeffType k = a[t + m] / b[m];
            ret.c[t] = k;
            for(int p = 0; p < m; ++p)
                a[t + p] -= k * b[p];
        }

        return ret;
    }

    Poly Poly::operator*(const Poly& r) const {
        Poly ret;
        const CoeffVec &a = c;
        const CoeffVec &b = r.c;

        int n = a.size() - 1;
        int m = b.size() - 1;

        ret.c = CoeffVec(n + m + 1);

        for(int t = 0; t <= m; ++t)
            for(int p = 0; p <= n; ++p)
                ret.c[p + t] += a[p] * b[t];

        return ret;
    }

    Poly Poly::operator+(const Poly& r) const {
        //cerr << "operator + " << endl;
        Poly ret;
        int n = c.size();
        int m = r.c.size();

        ret.c = CoeffVec(max(m, n));
        for(int t = 0; t < min(m,n); ++t)
            ret.c[t] = c[t] + r.c[t];

        for(int t = min(m,n); t < max(m,n); ++t)
            if(n > m)
                ret.c[t] = c[t];
            else
                ret.c[t] = r.c[t];

        return ret;
    }

    Poly Poly::operator-(const Poly& r) const {
        Poly ret;
        int n = c.size();
        int m = r.c.size();

        ret.c = CoeffVec(max(m, n));
        for(int t = 0; t < min(m,n); ++t)
            ret.c[t] = c[t] - r.c[t];

        for(int t = min(m,n); t < max(m,n); ++t)
            if(n > m)
                ret.c[t] = c[t];
            else
                ret.c[t] = - r.c[t];

        return ret;  
    }

    template<typename T>
    void Poly::fst(T x, T *pv, T *pdv, T *pddv, double *perr) const {
        T f = 0;
        T s = 0;
        T u = 0;
	double err = abs(c.back());
	double abx = abs(x);

        for(int t = c.size() - 1; t >= 0; --t)
        {
            u = u * x + s;
            s = s * x + f;
            f = f * x + c[t];
	    err = abs(c[t]) + abx * err;
        }

        if(pv) *pv = f;
        if(pdv) *pdv = s;
        if(pddv) *pddv = 2. * u;
	if(perr) *perr = err;
    }

	template<typename T>
	bool isnan(T x) {
		return !(x == x);
	}

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

            r.fst(x, &v, &dv, &ddv, 0);
      
            int n = r.deg();//.size() - 1;
            T g = dv / v;
            T h = g * g - ddv / v;
            T ddenum = sqrt(Poly::RootType(n - 1.0) * (h * Poly::RootType(n) - g * g));
            T denum1 = g + ddenum;
            T denum2 = g - ddenum;
            T denum = 0;
            if(abs(denum1) > abs(denum2))
                denum = denum1;
            else
                denum = denum2;

            a = Poly::RootType(n) / denum;

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

    template<typename T>
    T laguerr_from_book(const Poly& r)
    {
      const int MR=8, MT=10, MAXIT = MR * MT;
      const double EPS=numeric_limits<double>::epsilon();
      static const double frac[MR + 1] = 
	{0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};


      Poly::RootType x, dx, x1, b, d, f, g, h, sq, gp, gm, g2, its;

      x = 0;

      int m = r.deg();
      for(int iter = 1; iter <= MAXIT; ++iter){
	its = iter;
	b = r.c[m];
	double err;
	d = f = 0;
	r.fst(x, &b, &d, &f, &err);
	f *= 2;
	err *= EPS;
	if (abs(b) <= err)
	  return x;
	g = d / b;
	g2 = g * g;
	h = g2 - 2.0 * f / b;
	sq = sqrt(double(m - 1) * (double(m) * h - g2));
	gp = g + sq; 
	gm = g - sq;
	double abp = abs(gp);
	double abm = abs(gm);
	if (abp < abm) 
	  gp = gm;
	dx = max(abp, abm) > 0.0 ? double(m) / gp : polar(1 + abs(x), double(iter));
	x1 = x - dx;
	if (x == x1) return x;
	if (iter % MT != 0) 
	  x = x1;
	else
	  x -= frac[iter/MT]*dx;
      }

      return 0;
    }

     
    Poly::RootVec
    Poly::all_roots() const
    {
#ifdef _TEST_
        cerr << "all roots" << endl;
        cerr << "start from " << (*this) << endl;
#endif

        RootVec ret;
        Poly lt(*this);
        while(lt.deg() > 0)
        {
            RootType root = laguerr<RootType>(lt);
#ifdef _TEST_ 
            cerr << "\troot is " << root << endl;
#endif

            ret.push_back(root);

            lt = lt / Poly((CoeffVec const&)(list_of(-root)(1)));
            

#ifdef _TEST_
            cerr << "\tpoly is " << lt << endl;
#endif

        }
        return ret;
    }
  

  /*
  Poly::RootVec
  Poly::all_roots()const{
    Poly::RootVec ret(deg());
    vector<double> coeff(c.size());
    for(int t = 0; t < int(c.size()); ++t)
      coeff[t] = c[t];

    vector<double> res(2 * deg());
    
    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(c.size());

    gsl_poly_complex_solve(&coeff[0], c.size(), w, &res[0]);

    gsl_poly_complex_workspace_free(w);

    for(int t = 0; t < int(deg()); ++t){
      ret[t] = complex<double>(res[2 * t], res[2 * t + 1]);
    }

    return ret;
  }
  */

#ifdef _TEST_
    /*
    template<int s>
    Poly ar2pol(CoeffType ar[s])
    {
        Poly ret;
        ret.dat.resize(s);
        for(int t = 0; t < s; ++t)
            ret.dat[t] = ar[t];
        return ret;
    }

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
        cerr << "laguerr test = " << laguerr<RootType>(lt) << endl;

        cerr << "full laguerr: " << endl;
        poly_type roots = all_roots(lt);
        for(int t = 0; t < roots.size(); ++t)
            cerr << roots[t] << " ";
        cerr << endl;

    }
    */
#endif

}  // namespace farn
