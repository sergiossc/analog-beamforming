// BinLFSR.h -- declarations for BinLFSR class
// Todd K. Moon
#ifndef BinPolyDiv_H
#define BinPolyDiv_H

#ifndef MIN
#define MIN(a,b) (a<b? a:b)
#endif

class BinPolyDiv {
   int g;  // divisor representation
   int p;  // degree of divisor
   int mask; // used to grab desired bits
   int state; // state of shift register
public:
   BinPolyDiv(void) { g=p=state=mask=0;}// default constructor
   BinPolyDiv(int g, int p);// constructor
   // g = g[0] + g[1]x + ... + g[p]x^p,
   // where g is represented by a bit field, 
   // g = g[p]|g[p-1]| ... | g[1]|g[0]
   // e.g., g = 0x23 = 10 0011 => D^5 + D + 1
   ~BinPolyDiv() {};                        // destructor
   void clearstate(void) { state = 0; };
   int div(unsigned char *d, int n, unsigned char *q, int &quotientdegree,
 	   int &remainderdegree);
   // divide d=d[0]+d[1]x + ... + d[n]x^n by g, and return quotient in q
   // and remainder as return value
   // d has d0 first, so it starts at the 
   // upper end of d and works back
   int remainder(unsigned char *d, int n, int &remainderdegree);
   // compute only the remainder, returning it as return value

   int rem1(unsigned char *d, int n);
   // rem1 does not bother to compute the remainder degree.
   // be computed.  It computes only the remainder
   int rem2(unsigned char *d, int n, int &remainderdegree);
   // rem2 skips the initial shift in, and may be used for 
   // incremental operations preceded either by remainder or rem1
};
#endif

/*
Local Variables:
compile-command: "g++ -c -g BinPolyDiv.cc"
End:
*/
