// galtest.cc -- test the Gallager (low-density parity-check) code
// decoder
// This program tests using the simple example in the writeup
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "galdec.h"
#include <iostream>
using namespace std;

main()
{
   int i;
   int numloops;
   int success;
   GALDEC galdec("Asmall.txt",1);
   galdec.printsparseA();
   double y[] = {-0.63,-0.83,-0.73,-0.04,0.1,0.95,-0.76,.66,-0.55,.58};
   double *p = new double[sizeof(y)/sizeof(double)];
   galdec.meastoprobAWGN(y,p,2,2);
   VECDUMP(p,galdec.N);
   success = galdec.decode(p,1,10,numloops);  // decode and print
   cout << "Success=" << success << 
	  "  number of loops=" << numloops << endl;
}

/*
Local Variables:
compile-command:"g++ -o galtest -g galtest.cc  galdec.cc"
End:
*/
