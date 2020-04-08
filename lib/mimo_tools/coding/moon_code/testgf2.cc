//  Program: testgf2.cc
//  test basic operations on GF2 class
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "GF2.h"

main()
{
   GF2 a,b,c;
   int b1 = 1;
   char b2 = 1;
   unsigned char b3=1;
   a = b1;
   b = b2;
   c = b3;
   cout << a << endl;
   cout << b << endl;
   cout << "a+b=" << a+b << endl;
   cout << "a+=b=" << (a+= b) << endl;
   a = b1;
   cout << "a+1=" << a+1 << endl;
   cout << "1+a=" << 1+a << endl;
   cout << "a=" << a << "  b=" << b << endl;
   cout << "a*b=" << a*b << endl;
   cout << "a*=b=" << (a*=b) << endl;
   a = b1;
   cout << "a*1=" << a*1 << endl;
   cout << "1*a=" << 1*a << endl;
   cout << "a^3=" << (a^3) << endl;
   cout << "a/b=" << a/b << endl;
   cout << "a/=b=" << (a/=b) << endl;
   a = b1;
   cout << "1/a=" <<  (1/a) << endl;
   cout << "a==b" << (a==b) << endl;
   cout << "a!=b" << (a!=b) << endl;
   cout << "-a=" << -a << endl;
}

/*
Local Variables:
compile-command: "g++ -g -o testgf2 testgf2.cc"
End:
*/


