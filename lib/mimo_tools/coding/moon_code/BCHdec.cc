// BCHdec.cc -- a general BCH decoder
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "BCHdec.h"
//#include "berlmass2.cc"
#include "berlmassBCH2.cc"

//#include "berlmassBCH.cc"
//#include "polynomialT.cc"

//#include "gcdpoly.cc"

BCHdec::BCHdec(int t_in, int n_in, GFNUM2m A1_in) : // constructor
   Searcher(t_in)
{
   t = t_in;
   t2 = 2*t_in;
   n = n_in;
   A1 = A1_in;
   s = new GFNUM2m[2*t];
   Lambda = new GFNUM2m[t+1];
}

void
BCHdec::decode(GF2 *r, GF2 *dec)
{
   // Fill in the blanks

   // Step 1: evaluate the syndromes

   // Step 2: Determine the error locator polynomial

   // Step 3: Find the roots of the error locator polynomial

   // Step 4: Find error values: Not necessary for binary BCH codes

   // Step 5: Correct errors corresponding to inverse roots
}


/*
Local Variables:
compile-command: "g++ -c -g BCHdec.cc"
End:
*/
	  
