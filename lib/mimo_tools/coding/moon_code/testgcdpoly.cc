//  Program: testpolygcd.cc

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "ModArnew.h"
#include "polynomialT.cc"
#include "gcdpoly.cc"

// make some typedefs for shorthand
typedef ModAr<5> Z5;
typedef polynomialT<double> polyR;
typedef polynomialT<Z5> polyZ5;

// cause template classes to instantiate with the desired types
template class polynomialT<double>;
template class polynomialT<Z5>;


main()
{

   Z5 d1[] = {1,0,0,1,3,0,4,3}; // 1+x^3 + 3x^4 + 4x^6 + 3x^7
   Z5 d2[] = {0,1,0,1,4};       // x+x^3+4x^4

   polyZ5 a(7,d1);      // build polynomial with d1 coefficients
   polyZ5 b(4,d2);      // build polynomial with d2 coefficients
   polyZ5 s, t, g;
   a.setprintdir(1);

   cout << "a=" << a << endl;
   cout << "b=" << b << endl;
   gcd(a,b,g,s,t);                  // compute the gcd of a and b
   cout << "g=" << g << endl; 
   cout << "s=" << s << endl;
   cout << "t=" << t << endl;

   double d1d[] = {2,8,10,4};      // 2+8x+10x^2+4x^3
   double d2d[] = {1,7,14,8};      // 1+7x+14x^2+8x^3
   polynomialT<double> ad(3,d1d);   // build polynomial with d1d coefficients
   polynomialT<double> bd(3,d2d);   // build polynomial with d2d coefficients
   polynomialT<double> sd,td,gd;
   cout << "a=" << ad << endl;
   cout << "b=" << bd << endl;
   gcd(ad,bd,gd,sd,td);             // compute gcd
                                    // (this will automatically pick
                                    // out the 'double' instantiation
   cout << "g=" << gd << endl; 
   cout << "s=" << sd << endl;
   cout << "t=" << td << endl;
}

/*
Local Variables:
compile-command:"g++ testgcdpoly.cc -g -o testgcdpoly"
End:
*/
