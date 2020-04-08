// BPSKmod.h --- a simple BPSK modulator to modulate a scalar of data
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BPSKMOD_H
#define  BPSKMOD_H

class BPSKmod 
{
public:
   double a;					// signal amplitude
   BPSKmod(double a_in=1) {
	  a = a_in;
   }
   ~BPSKmod() {}
   double mod(unsigned char bitin) {
	  double v = a*(2*bitin-1);
	  return v;
   }
};


#endif

/*
Local Variables:
compile-command: "g++ -c BPSKmod.h"
End:
*/

