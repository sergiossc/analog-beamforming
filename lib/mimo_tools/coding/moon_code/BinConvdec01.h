// BinConvdec01.h -- Convolutional decoder for binary (0,1) data
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BinConvdec01_H
#define BinConvdec01_H
#include "Convdec.h"
#include "matalloc.h"

class BinConvdec01 : public Convdec {
   unsigned char ***outputmat;  // [numstates][numbranches][n] -- outputs
   unsigned char *data;         // pointer to input data
public:
   BinConvdec01(BinConv & encoder, int in_pathmem) 
      : Convdec(encoder,in_pathmem)
      {  
         buildoutputmat(encoder);       // set up the outputmat
      };
   ~BinConvdec01() { FREETENSOR(outputmat,numstates);  delete[] data;   }

   // Compute the metric using the Hamming distance
   virtual double metric(unsigned int state, int branch) {
      double sum=0;
      for(int i = 0; i < n; i++) {
         sum += (data[i] != outputmat[state][branch][i]);
      }
      return sum;
   };

   virtual double metric_dfree(unsigned int state, int branch)
   {
	  return 0.0;  // this does nothing here.
   }

   int decode(unsigned char *outs) {
      data = outs;
      return viterbi();
   };
   void buildoutputmat(BinConv& encoder);
                              // Fill in the blanks
};
#endif
/*
Local Variables:
compile-command: "g++ -c BinConvdec01.cc"
End:
*/

