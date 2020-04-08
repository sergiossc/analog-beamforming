// galdec.cc -- Galager decoder class definitions
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only


#include "galdec.h"
#include <fstream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>

// const double TINYDIV = 1.e-10;
// const double EPS = 0;
// const double CLIPONE = 1-EPS;
// #define CLIPCHECK(q0,q1,one) if(q0>one){q0=one;q1=1-one;} \
//                              else if(q1>one){q1=one;q0=1-one;}

GALDEC::GALDEC(char *fname, int offset)
{
   int n,m,k,d;
   const int maxline = 1024;	// maximum length of line expected 
   char line[maxline];

   ifstream infile(fname);
   if(!infile) {
	  cerr << "Error: unable to open input file" << fname << endl;
	  exit(-1);
   }
   infile >> N;
   infile >> M;
   K = N-M;
   infile >> maxcolwt;
   infile >> maxrowwt;
   CALLOCMATRIX(Mn,int,N,maxcolwt);
   CALLOCMATRIX(Nm,int,M,maxrowwt);
   // read in the column weights
   Mnlen = new int[N];
   for(n = 0; n < N; n++) {
	  infile >> Mnlen[n];
   }
   // read in the row weights
   Nmlen = new int[M];
   for(m = 0; m < M; m++) {
	  infile >> Nmlen[m];
   }
   // read in the Mn data (checks for each bit)
   for(n = 0; n < N; n++) {
	  for(m = 0; m < maxcolwt; m++) {
		 infile >> d;
		 Mn[n][m] = d - offset;
	  }
   }

   // read in the Nm data (bits for each check)
   for(m = 0; m < M; m++) {
	  for(n = 0; n < maxrowwt; n++) {
		 infile >> d;
		 Nm[m][n] = d-offset;
	  }
   }
   allocdecodedat();
}

void
GALDEC::allocdecodedat(void)
// allocate memory used in the decoder
{
   na = new int[N];
   deltaq = new double[maxrowwt];
   CALLOCMATRIX(r1,double,maxcolwt,N);
   CALLOCMATRIX(r0,double,maxcolwt,N);
   CALLOCMATRIX(q1,double,maxcolwt,N);
   CALLOCMATRIX(q0,double,maxcolwt,N);
   CALLOCMATRIX(deltar,double,maxcolwt,N);
   q0p = new double[N];
   q1p = new double[N];
   c = new unsigned char[N];
}

void
GALDEC::freedecodedat(void)
{
   delete [] c;
   delete[] q1p;
   delete[] q0p;
   FREEMATRIX(deltar);
   FREEMATRIX(q0);
   FREEMATRIX(q1);
   FREEMATRIX(r0);
   FREEMATRIX(r1);
   delete[] deltaq;
   delete[] na;

   delete[] Mnlen;
   delete[] Nmlen;
}

   
int
GALDEC::decode(double *pn, int printstuff, int maxnumloop, int &numloops)
// pn = channel posterior probabilities (one for each N)
// printstuff - set to print intermediate results
// maxnumloop = maximum number of decoding iterations
// numloops = number of decoding iterations actually used
// returns: 1 if decoding succeeds; 0 for a decoding failure
{
   int i,l,k,m,n,row;
   double prod, prod0, prod1,alpha;
   int paritycheck = 0;
   char z;
   int loopcount = 0;
   double sum;
   int idx;

   // initialize the q data
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 q1[m][n] = pn[n];
	  }
   }
   if(printstuff) printsparse("q1",q1);


   do {  // loop until completion of decoding, or loop count
	  for(i = 0; i < N; i++) {
		 na[i] = 0;				// clear the "number above" offset
	  }
	  loopcount++;
	  if(printstuff) cout << "Iteration: " << loopcount << endl;
	  // Horizontal step

	  // Fill in the blanks ...



	  if(printstuff) printsparse("r1",r1);
	  
	  // Vertical step

	  // Fill in the blanks ...

	  // Decode using the pseudoposterior probabilities

	  // Fill in the blanks ...

	  if(printstuff) printsparse("q1",q1);
	  if(printstuff) VECDUMP(q1p,N);
	  if(printstuff) VECDUMP2(c,N,int);


	  // Check the parity condition 
	  paritycheck = 1;
	  for(m = 0; m < M; m++) {
		 z = 0;
		 for(l = 0; l < Nmlen[m]; l++) {
			z += c[Nm[m][l]];
		 }
		 z %= 2;
		 if(z) {
			// Parity check fails.  We could bail out of check at this 
			//   point.
			paritycheck = 0;
			break;	// break out of parity check loop
		 }
	  }
	  if(paritycheck) break;	// break out of decode loop
   } while(loopcount < maxnumloop);
   numloops = loopcount;
   return(paritycheck);
}


void GALDEC::printsparse(char *name,double **sparsemat)
{
   int l,n, m,lastn,n1;
   int *na;

   int sp = cout.precision();
   cout.precision(3);
   na = new int[N];
   for(l = 0; l < N; l++) na[l] = 0;
   cout << name << endl;
   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << "0\t";
			// printf("%.2f\t",0.0);
		 }
		 // print the new number
		 cout << sparsemat[na[n]++][n] << "\t";
		 // printf("%.2f\t",sparsemat[na[n]++][n]);
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << "0.0\t";
		 // printf("%.2f\t",0.0);
	  }
	  cout << endl;
	  // printf("\n");
   }
   delete[] na;
   cout.precision(sp);
}

void
GALDEC::printsparseA(void)  // print the sparse binary A matrix
{
   int l,n, m,lastn,n1;
   int lastm,m1;

   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << "0 ";
			// printf("0 ");
		 }
		 // print the new number
		 cout << "1 "; 
		 // printf("1 ");
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << "0 ";
		 // printf("0 ");
	  }
	  cout << endl;
	  // printf("\n");
   }
}

void
GALDEC::meastoprobAWGN(double *y,double *p,double a, double sigma2)
// convert measured data to a probability.  In this case, for an AWGN
// with variance sigma2 and channel modulation amplitude a
{
   for(int i = 0; i < N; i++) {
	  p[i] = 1/(1+exp(-2*a*y[i]/sigma2));
   }
}

/*
Local Variables:
compile-command: "g++ -c -g galdec.cc"
End:
*/
