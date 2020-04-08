// BinConvdecBPSK.h -- Convolutional decoder for BPSK modulated data
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BinConvdecBPSK_H
#define BinConvdecBPSK_H
#include "Convdec.h"
#include "matalloc.h"

class BinConvdecBPSK : public Convdec {
   double ***outputmat; // [numstates][numbranches][n] -- outputs
   double *data;            // pointer to input data
public:
   BinConvdecBPSK(BinConv & encoder, int in_pathmem) 
      : Convdec(encoder,in_pathmem)
      {  
         buildoutputmat(encoder);       // set up the outputmat
      };
   ~BinConvdecBPSK() { FREETENSOR(outputmat,numstates);  delete[] data;   }

   virtual double metric(unsigned int state, int branch) {
      double t,sum=0;
      for(int i = 0; i < n; i++) {
         t = data[i]-outputmat[state][branch][i];
         sum += t*t;
      }
      return sum;
   };
   virtual double metric_dfree(unsigned int state, int branch) {
	  return 0.0;  // this does nothing here.
   }

   int decode(double *outs) {
      data = outs;
      return viterbi();
   };

   void buildoutputmat(BinConv& encoder); // to be completed
};

#endif
/*
Local Variables:
compile-command: "g++ -c BinConvdecBPSK.cc"
End:
*/
