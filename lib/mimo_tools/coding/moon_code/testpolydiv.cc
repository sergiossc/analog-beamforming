// testpolydiv.cc -- test division operations
// Copyright 2003 Todd K. Moon, Utah State University

#include "ModAr.h"
#include "polynomialT.cc"

template class polynomialT<ModAr>;
template class polytemp<ModAr>;


#include <iostream.h>

main()
{


   // Now check polynomials in Z_5[x]

   ModAr::setdefaultm(5);
   ModAr::setshowmod(0);

   ModAr d1m[] = {3,4,3}; // 3 + 4x + 3x^2
   ModAr d2m[] = {1,3};   // 1+3x
   polynomialT<ModAr> p1m(2,d1m);
   polynomialT<ModAr> p2m(1,d2m);
//    ModAr d1m[] = {3,3}; // 3+3x
//    ModAr d2m[] = {3};   // 3;
//    polynomialT<ModAr> p1m(1,d1m);
//    polynomialT<ModAr> p2m(0,d2m);


   cout << "\n\n\nPolynomial mod 5\n";
   cout << "p1m=" << p1m << endl;
   cout << "p2m=" << p2m << endl;
   polynomialT<ModAr> qm = p1m/p2m;
   polynomialT<ModAr> rm = polynomialT<ModAr>::getlastremainder();
   cout << "quotient (p1/p2): " << qm << endl;
   cout << "remainder (p1/p2): " << rm << endl;

   ModAr::setdefaultm(2);

   
}

/*
Local Variables:
compile-command:"g++ -o testpolydiv -g -fno-implicit-templates testpolydiv.cc instantiate.cc ModAr.cc"
End:
*/
