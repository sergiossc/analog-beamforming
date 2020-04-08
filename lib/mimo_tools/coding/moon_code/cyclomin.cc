//
//
//  Program: cyclomin.cc --- compute cyclotomic sets and 
//   minimal polynomials for GF(2^m)
//
//  Todd K. Moon
//  Utah State University
//

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <iostream>
using namespace std;
#include "GFNUM2m.h"
#include "polynomialT.cc"

// template class polynomialT<GFNUM2m>;

main(int argc, char **argv)
{

   if(argc==1) {
	  cerr << "Prints cyclotomic cosets modulo 2^m-1 and minimal polynomials";
	  cerr << endl;
	  cerr << "Usage: " << argv[0]  << " m \n";
	  return(-1);
   }
   int m = atoi(argv[1]);
   GFNUM2m::initgf(m);			// initialize the field
   int N = (1<<m);
   int N1 = N-1;
   polynomialT<GFNUM2m> p;		// minimal polynomial 
   polynomialT<GFNUM2m> l1; { GFNUM2m l1a[] = {0,1}; l1.setc(l1a,1); }

   int *values = new int[N];
   int i;
   for(i = 0; i < N; i++) {
	  values[i] = 0;
   }
   int q = 2;
   int a;
   for(i = 0; i < N1; i++) {
	  if(values[i]) continue;
	  a = i;
	  cout << "{" << a;
	  l1[0] = A^a;
	  p = l1;
	  values[i] = 1;
	  while(1) {
		 a = (a*q) % N1;
		 if(a != i) {
			values[a] = 1;
			cout << " " << a;
			l1[0] = A^a;
			p = p * l1;
		 }
		 else {
			break;
		 }
	  }
	  cout << "}" << endl;
	  cout << "M(x)=" << p << endl;
   }
}



// This stuff is placed here, where hopefully it won't be too 
// readily discovered, since it provides a solution to a
// lab exercise

// Define static variables in GFNUM2m
int *GFNUM2m::p2v = 0;			// convert exponent to vector
int *GFNUM2m::v2p = 0;			// a list of elements to convert from
								// vector to exponential notation
int GFNUM2m::gfm = 0;			// vector size of field element
int GFNUM2m::gfN = 0;			// number of nonzero-elements in field
outformat GFNUM2m::outtype = power;
								// default to exponential output
// ALPHA elements from the field
GFNUM2m ALPHA;                  // define the element alpha
GFNUM2m& A= ALPHA;              // and a shorthand reference to it


void GFNUM2m::initgf(int m)
{
   //do the initialization by using only the size of the field
   // m in GF(2^m).
   // A fixed set of primitive polynomials is used.
   // Octal:
   unsigned int g[] = {1,1,7,013, 023, 045, 0103, 0211, 0435, 01021, 02011,
					   04005, 010123, 020033, 042103, 0100003};
   if(m>sizeof(g)/sizeof(unsigned int)) {
	  cerr << "Error: must specify connection polynomial for m" << endl;
   }
   GFNUM2m::initgf(m,g[m]);
}

// initgf: (1) Build the v2p and p2v tables
//         (2) set the global variable ALPHA
//         (3) set the static member variables gfm and gfN
void GFNUM2m::initgf(int m,unsigned int g)
// m -- GF(2^m)
// g -- generator polynomial, bits represent the coefficients:
// e.g. g = 0x13 = 1 0011 = D^4 + D + 1
{
   int i,j;

   if(m > sizeof(unsigned int)*8) {	// too many bits!
	  cerr << "Error: Degree too large in GFNUM2m" << endl;
   }
	  
   ALPHA.v = 2;					// set up alpha element
   gfm = m;
   gfN = (1<<m)-1;				// gfN = number of nonzero field elements,
								// gfN = 2^n -1

   if(v2p) delete[] v2p;		// delete any prior stuff
   if(p2v) delete[] p2v;
   
   v2p = new int[gfN+1];  		// table to convert vector to power form
   p2v = new int[gfN+1];		// tabel to convert power to vector form

   // Fill in the blanks:
   // Fill in the tables for v2p and p2v, using the LFSR ...
}


ostream& operator<<(ostream& s,const GFNUM2m& arg)
{
   if(arg.getouttype() == power) {
	  if(arg.v < 2) return s << arg.v;
	  else {
		 int e = arg.v2p[arg.v];
		 if(e == 1) return s << "A";
		 else return s << "A^" << e;
	  }
   }
   else {						// vector (numeric) output
	  return s << arg.v;
   }
}


  
/*
Local Variables:
compile-command: "g++ -o cyclomin -g cyclomin.cc"
End:
*/


