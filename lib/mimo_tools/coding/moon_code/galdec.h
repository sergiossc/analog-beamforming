// galdec.h --- Galager decoder class declarations
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef GALDEC_H
#define GALDEC_H
#include "matalloc.h"

class GALDEC {
public:  // for debugging purposes
   int N;						// block length
   int K;						// message length
   int M;						// N-K -- redundancy
   int *Nmlen;					// [M] lengths of row weight vectors
   int **Nm;					// [M][rowwt]
                                //  set of bits that participate in check m
   int *Mnlen;					// [N] lengths of column weight vectors
   int **Mn;					// [N][colwt]
                                // set of checks in which bit n participates
   int maxcolwt;				// maximum weight of columns
   int maxrowwt;				// maximum weight of rows

   int *na;						// [N] # of elements above this in column
   double *deltaq;				// [maxrowwt] temporary row info holding deltaq
   double **r1;					// [maxcolwt][N]  r1[m][n]
   double **r0;					// [maxcolwt][N]  r0[m][n]
   double **q0;					// [maxcolwt][N]  q0[m][n]
   double **q1;					// [maxcolwt][N]  q1[m][n]
   double **deltar;				// [maxcolwt][N]  deltar[m][n]
   double *q0p;					// [N] -- pseudopriors
   double *q1p;					// [N] -- pseudopriors
   unsigned char *c;			// decoded bit values

   void allocdecodedat(void);	// allocate the decoder memory
   void freedecodedat(void);	// free the decoder memory
public:
   GALDEC(char *fname, int offset=0); // constructor --- read from file
   //  Some definitions use base index=1, while some uses base index=0.
   // Use offset=1 when reading a file with base index=1
   // Use offset=0 when reading a file with base index=0
   ~GALDEC() { freedecodedat(); };
   int decode(double *pn, int printstuff, int maxnumloop, int & numloops);
   void printsparse(char *name, double **sparsemat);
   void printsparseA(void);
   void meastoprobAWGN(double *y, double *p, double a, double sigma2);
   // convert measured data to a probability.  In this case, for an AWGN
   // with variance sigma2 and channel modulation amplitude a
};

#endif
/*
Local Variables:
compile-command: "g++ -c -g galdec.cc"
End:
*/
