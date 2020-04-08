//  Program: testturbodec.cc
//
//  Todd K. Moon
//
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <iostream>
using namespace std;
#include <math.h>

#include "BinConvIIR.h"
// #include "matalloc.h"
#include <iostream>
using namespace std;
#include "BPSKmodvec.h"
#include "Turboenc.h"
#include "Turbodec.h"
double gran(void);

int main()
{
   int i,j;
   int decout;

   int interleaveseed = 1;
   int k = 1;
   int n = 2;
//    unsigned int gnum = 0x11;   // 1 0001
//    unsigned int gden = 0x1F;     // 1 1111
//    int p = 4;					// degree of denominator
   unsigned int gnum = 4;		// 100 = 1
   unsigned int gden = 5;		// 101 = 1+D^2
   int p = 2;

//   unsigned char d[] = {1,0,0,0,0,0,0,0,0,0};  // for the impulse response
   unsigned char d[] = {1,1,0,0,1,0,1,0,1,1};
   int N = sizeof(d)/k;			// number of time steps in test array
   unsigned char **P=0;			// puncture matrix
//    CALLOCMATRIX(P,unsigned char,2,2);
//    P[0][0] = 1;  P[0][1] = 0;
//    P[1][0] = 0;  P[1][1] = 1;
   Turboenc encoder(p,gnum,gden,N,interleaveseed,P,2);
   unsigned char *out;

   unsigned char *encout;
   encout = encoder.encode(d);
   unsigned int finalstate1 = encoder.getstate1();
   unsigned int finalstate2 = encoder.getstate2();
   cout << "Final state1=" << finalstate1 << endl;
   cout << "Final state2=" << finalstate2 << endl;
   cout << "Encoded data:" << endl;
   if(P) {
	  for(i = 0; i < 2*N; i+=2) {
		 cout << int(encout[i]) << int(encout[i+1]) << " ";
	  }
   }
   else {
	  for(i = 0; i < 3*N; i+=3) {
		 cout << int(encout[i]) << int(encout[i+1]) << int(encout[i+2]) << " ";
	  }
   }
   cout << endl;
   int cblocklen = int(N/encoder.R);
   BPSKmodvec modulator(cblocklen);
   double *allouts;
   allouts = modulator.mod(encout);
   if(P) {
	  for(i = 0; i < cblocklen; i+=2) {
		 cout << allouts[i] << allouts[i+1] << " ";
	  }
   }
   else {
	  for(i = 0; i < 3*N; i+=3) {
		 cout << allouts[i] << allouts[i+1] << allouts[i+2] << " ";
	  }
   }
   cout << endl;
   // add the noise
   double sigma2 = 0.45;
   for(i = 0; i < cblocklen; i++) {
	  allouts[i] += sqrt(sigma2)*gran();
   }
//   sigma2 = 1;
   cout << "Outputs with noise: " << endl;
   if(P) {
	  for(i = 0; i < cblocklen; i+=2) {
		 cout << allouts[i] <<" " << allouts[i+1] << " | ";
	  }
   }
   else {
	  for(i = 0; i < 3*N; i+=3) {
		 cout<<allouts[i] << " " << allouts[i+1] <<" "<<allouts[i+2] << " | ";
	  }
   }
   cout << endl;
   Turbodec decoder(p,gnum,gden,N,sigma2,interleaveseed,P,2);
   decoder.decode(allouts,1,finalstate1,finalstate2);
}


/*
Local Variables:
compile-command: "g++ -o testturbodec -g testturbodec.cc BCJR.cc BinConvIIR.cc interleave.cc Turboenc.cc Turbodec.cc gran.cc"
End:
*/


