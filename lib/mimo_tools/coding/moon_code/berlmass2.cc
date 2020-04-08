//  Program: berlmass2.cc --- implement the Berlekamp-Massey algorithm
//  using arrays, not polynomials
//  Todd K. Moon

#include "polynomialT.h"

#define MAX2(a,b) (a> b? a:b)

template <class T> void
berlmass2(const T* s, int n, T* c, int& L)
// s = input coefficients s[0],s[1],... s[n-1]
// c = connection polynomial coefficients.  Must be allocated prior to calling
// L = degree of connection polynomial
{
   // fill in the blanks ..
}

#undef MAX2
/*
Local Variables:
compile-command: "g++ -g berlmass.cc"
End:
*/


