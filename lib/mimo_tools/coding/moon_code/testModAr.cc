// ******************************************************************
// testModAR.cc -- test the modularized arithmetic
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "ModAr.h"

#include <iostream>
using namespace std;

main()
{
   ModAr::setdefaultm(5);		// set the default modulo
   ModAr::setshowmod(0);		// set so modulo is not displayed

   ModAr a(3);
   ModAr b(4);
   ModAr c;
   int aint;
   c = a+b;
   cout << "a=" << a << "  b=" << b << "  c=a+b=" << c << "\n";
   a = 2;
   cout << "a=" << a << "\n";
   c = a+1;
   cout << "a+1=" << a+1 << "\n";
   cout << "c=" << c << "\n";
   aint = a.getv();
   ModAr f = -a;
   cout << "f=-a= " << f << "\n";
   cout << "f^5= " << (f^5) << "\n";
   cout << "a=" << a << "  b=" << b << "\n";
   cout << "a*b= " << a*b << "\n";
   cout << "a/b= " << a/b << "\n";
   a /= b;
   cout << "a /= b: " << a << "\n";
}

/*
Local Variables:
compile-command:"g++ -o testModAr -g testModAr.cc ModAr.cc"
End:
*/
