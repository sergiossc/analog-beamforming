//
//  Program: testsugiyama.cc

#include <math.h>
#include "ModAr.h"
#include "polynomialT.cc"
#include "gcdpoly.cc"

// create instantiations of the polynomial class of type ModAr
template class polynomialT<ModAr>;
template class polytemp<ModAr>;
// create instantiations of the polynomial class of type double
template class polynomialT<double>;
template class polytemp<double>;

// Create an instantiation of the gcd function, for polynomials of type ModAr
template <class ModAr> void
gcd(const polynomialT<ModAr> &a, const polynomialT<ModAr> &b, 
	polynomialT<ModAr> &g,
	polynomialT<ModAr> &s, polynomialT<ModAr> &t, int sdeg);

// template <double> void
// gcd(const polynomialT<double> &a, const polynomialT<double> &b, 
// 	polynomialT<double> &g,
// 	polynomialT<double> &s, polynomialT<double> &t, int sdeg);


main()
{
   int i,j;
   ModAr::setdefaultm(5);
   ModAr::setshowmod(0);

   ModAr sum = 0;
   int p = 3;

   POLYC(ModAr,b,{2,3,4,2,2,3});
   POLYC(ModAr,t1,{1,3,4,2});

   cout << "b=" << b << endl;
   cout << "t=" << t1 << endl;
   cout << "b*t1=" <<  b*t1 << endl;


   polynomialT<ModAr> a;
   polynomialT<ModAr> s, t, g;
   a.setc(2*p,1);  // x^{2p}

   cout << "a=" << a << endl;
   cout << "b=" << b << endl;
   gcd(a,b,g,s,t,p);
   cout << "g=" << g << endl; 
   cout << "s=" << s << endl;
   cout << "t=" << t << endl;
   cout << "b*t=" << b*t << endl;

   POLYC(ModAr,p1,{4,1,0,0,4,3,1});
   POLYC(ModAr,p2,{4,2,2,2,3});
   cout << "p1=" << p1 << "  p2=" << p2 << endl;
   gcd(p1,p2,g,s,t);
   cout << "g=" << g << endl; 
   cout << "s=" << s << endl;
   cout << "t=" << t << endl;
   cout << "b*t=" << b*t << endl;

   POLYC(ModAr,b1,{1,4,2,2,4,1});
   p = 3;  						// 2p = 6 coefficients
   a = ModAr(0);
   a.setc(2*p,1);
   gcd(a,b1,g,s,t,p);
   cout << "g=" << g << endl; 
   cout << "s=" << s << endl;
   cout << "t=" << t << endl;
   cout << "b*t=" << b1*t << endl;

}

/*
Local Variables:
compile-command:"g++ -o testsugiyama -g testsugiyama.cc ModAr.cc"
End:
*/


