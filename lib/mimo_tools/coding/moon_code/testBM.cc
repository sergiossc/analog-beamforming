//  Program: testBM.cc --- test the Berlekamp Massey algorithm
//  Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "GF2.h"                // provide binary arithmetic
#include "GFNUM2m.h"            // provide GF(2^m) arithmetic
#include "ModAr.h"              // Provide arithmetic mod 5
#include "polynomialT.cc"       // provide polynomial operations
#include "berlmass.cc"          // Berlekamp-Massey algorithm
#include "berlmass2.cc"         // Berlekamp-Massey, w/o polynomials

// create instantiations of the polynomial class of type GF2
template class polynomialT<GF2>;
template class polytemp<GF2>;

// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;

// create instantiations of the polynomial class of type ModAr
template class polynomialT<ModAr>;
template class polytemp<ModAr>;

main()
{
   int j;
   int l;
   // GF2 test
   GF2 d1[] = {1,1,1,0,1,0,0};           // array of data
   GF2 c1a[6];
   polynomialT<GF2> c = berlmass(d1,7);  // polynomial form
   cout << "c=" << c << endl;
   berlmass2(d1,7,c1a,l);                // array form
   for(j = 0; j <= l; j++) cout << c1a[j] << " "; cout << endl; cout << "\n\n";

   // ModAr test
   ModAr::setdefaultm(5);               // set up the field
   ModAr d2[] = {2,3,4,2,2,3};          // array of data
   ModAr c2a[6];
   polynomialT<ModAr> c2 = berlmass(d2,6); // polynomial form
   cout << "c=" << c2 << endl;
   berlmass2(d2,6,c2a,l);               // array form
   for(j = 0; j <= l; j++) cout << c2a[j] << " "; cout << endl; cout << "\n\n";


   GFNUM2m::initgf(4,0x13);             // Initialize field: 1 011 = 1+x+x^4
      
   GFNUM2m d4[] = {1,1,A^10,1,A^10,A^5};// array of data
   GFNUM2m c4a[6];
   polynomialT<GFNUM2m> c4 = berlmass(d4,6);// polynomial form
   cout << "c=" << c4 << endl;
   berlmass2(d4,6,c4a,l);               // array form
   for(j = 0; j <= l; j++) cout << c4a[j] << " "; cout << endl; cout << "\n\n";

}

/*
Local Variables:
compile-command: "g++ -o testBM -g testBM.cc ModAr.cc GFNUM2m.cc"
End:
*/
