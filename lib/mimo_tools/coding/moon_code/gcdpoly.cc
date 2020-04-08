// gcd.cc --- Implement a gcd algorithm for a general
// polynomial type.  Also, implement a gcd for polynomials with real
// coefficients, truncating coefficients as necessary to avoid
// roundoff problems.
//
//  Also, implement the Sugiyama algorithm 

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <math.h>
#include "polynomialT.h"

// declare the gcd function, setting up a default value for the rdeg
template <class T> void
gcd(const polynomialT<T> &a, const polynomialT<T> &b, 
	polynomialT<T> &g,
	polynomialT<T> &s, polynomialT<T> &t, int rdeg=0);

// Function definition: gcd
template <class T> void
gcd(const polynomialT<T> &a, const polynomialT<T> &b, polynomialT<T> &g,
	polynomialT<T> &s, polynomialT<T> &t, int rdeg)
{
   // Fill in the blanks ...
}


static void chop(polynomialT<double> *p, double eps);

// This is a specialization for doubles, since it has to handle
// the roundoff more carefully
template <> void
gcd(const polynomialT<double> &a, const polynomialT<double> &b, 
	polynomialT<double> &g,	polynomialT<double> &s, polynomialT<double> &t,
	int rdeg)
{
   // Fill in the blanks.
}

void chop(polynomialT<double> *p, double eps)
{
   int i;
   int newdegree=0;
   int done = 0;
   for(i = p->getdegree(); i>= 0; i--) {
	  if(fabs((*p)[i]) < eps) {
		 (*p)[i] = 0;
		 if(i== p->getdegree()) {
			newdegree = i;  // set that new degree is necessary
		 }
	  }
   }
   if(newdegree) {
	  for(i = newdegree-1; i > 0; i--) {
		 if((*p)[i] != 0) {
			p->resizecopy(i);
			done = 1;
			break;
		 }
	  }
	  if(!done) p->resizecopy(0);
   }
}
