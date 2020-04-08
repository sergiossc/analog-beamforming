// BCHdec.h -- a general BCH decoder
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BCHdec_H
#define BCHdec_h

#include "GF2.h"
#include "GFNUM2m.h"
#include "polynomialT.h"
#include "ChienSearch.h"

class BCHdec {
   int t; 						// error correcting capability
   int t2;						// 2*t
   int n;						// block length
   int nu;						// degree of error locator polynomial
   GFNUM2m A1;					// primitive element used
   ChienSearch Searcher;		// Chien search object
   GFNUM2m *s;					// syndrome storage
   GFNUM2m *Lambda;				// error locator polynomial
public:
   BCHdec(int t_in, int n_in, GFNUM2m A1_in = ALPHA);
   ~BCHdec() {  delete[] s; delete[] Lambda;};
   void decode(GF2 *r, GF2 *dec);
};



#endif

/*
Local Variables:
compile-command: "g++ -c -g BCHdec.cc"
End:
*/

