//  Program: testBCH.cc --- test the BCH decoder
//  Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "BCHdec.h"

// instantiate the berlekamp-massey algorithm for GFNUM coefficients
template <class GFNUM2m> void
berlmass2(const GFNUM2m* s, int n, GFNUM2m* c, int& L);

template <class GFNUM2m> void
berlmass2(const GFNUM2m* s, int n, GFNUM2m* c, int& L);


main()
{
   int i,j,l;
   int j1,j2,j3;
   GFNUM2m::initgf(4,0x13);  // 1 0011 = 1+d+d^4
   int t = 3;
   int n = 15;
   BCHdec decoder(t,n);
   GF2 r[n];    // the received vector
   GF2 dec[n];    // the decoded vector
   for(j = 0; j < n; j++) r[j] = 0;  // clear out previous
   for(j1 = 0; j1 < n; j1++) {
	  cout << "j1=" << j1 << endl;
	  for(j2 = 0; j2 < n; j2++) {
		 for(j3 = 0; j3 < n; j3++) {
			r[j1] = r[j2] = r[j3] = 1;
			decoder.decode(r,dec);
			for(i = 0; i < n; i++) {
			   if(dec[i] != 0) {
				  cout << "Uncorrected error: (" << j1 <<"," << j2 << ","
					   << j3 << ")" << endl;
				  break;
			   }
			}
			r[j1] = r[j2] = r[j3] = 0;  // clear out previous
		 }
	  }
   }
}

/*
Local Variables:
compile-command: "g++ -o testBCH -g testBCH.cc BCHdec.cc GFNUM2m.cc ChienSearch.cc"
End:
*/


