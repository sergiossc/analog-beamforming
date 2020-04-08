// Turbodec.h -- a Turbo decoder
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef TURBODEC_H
#define TURBODEC_H

#include "BinConvIIR.h"   // the systematic convolutional encoder
#include "interleave.h"   // the random interleaver
#include "Turboenc.h"			// the turbo encoder
#include "BCJR.h"				// the BCJR object

class Turbodec {
   Turboenc enc;				// encoder we are working on
   BCJR bcjr;					// BCJR object
   int blocklen;				// length of data block
   double **r1;					// data for encoder 1
   double **r2;					// data for encoder 2
   double **prior1, **prior2;	// locations for prior probabilities
   int numbranch;				// number of branches = 2^k
   unsigned char **P;			// puncture matrix
   int puncturelen;				// width of puncture matrix
   int puncturecycle;			// column counter for puncture matrix
   double R;					// code rate, default=1/3 with no puncturing
public:
   Turbodec(int deg, unsigned int h_in,unsigned int g_in,
			int in_blocklen,
			double in_simga2,
			unsigned int interleaveseed = 0,
			unsigned char **punctureP=0, int puncturelen=0);
   ~Turbodec() { 
	  FREEMATRIX(prior1); FREEMATRIX(prior2);
      FREEMATRIX(r1);  FREEMATRIX(r2);
   };
   double **decode(const double *r,int numit, 
			 unsigned int finalstate1=0, unsigned int finalstate2=-1);
   // turbo decoder function.  Returns a pointer to the array
   // of posterior probabilities

   // make this public so the channel variance can be set
   double sigma2;				// noise variance
};



#endif
/*
Local Variables:
compile-command: "g++ -c Turbodec.cc"
End:
*/

