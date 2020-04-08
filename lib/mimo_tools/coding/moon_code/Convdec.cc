//  Program: Convdec.cc -- Convolutional decoding
// using the Viterbi algorithm
//  Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "Convdec.h"
#include "BinConv.h"
#include "matalloc.h"
#define LARGE 1e99

Convdec::Convdec(BinConv &encoder, int in_pathmem)
{
   n = encoder.n;
   k = encoder.k;

   pathmem = in_pathmem;

   // build up the information needed to do the decoding
   numstates = (1<<encoder.nu);
   numbranches = (1<<k);
   metrics1 = new double[numstates];
   metrics2 = new double[numstates];

   buildprev(encoder);   // build the state/previous state array

   // Build up the path information
   // Fill in the blanks ...
   startstate = 0;              // default state state
   setpaths();                  // initialize the path variables
}

int
Convdec::viterbi()
// r is a array of the received values
{
   // Fill in the blanks
   // ...
}

int
Convdec::getinpnow(int adv)
// Get the best input from the current state of the trellis
// if adv is set, move forward
{
   // Fill in the blanks
   // ...
}

void 
Convdec::buildprev(BinConv &encoder)
{
   unsigned int savestate = encoder.getstate();

   CALLOCMATRIX(prevstate,unsigned int, numstates, numbranches);
   CALLOCMATRIX(inputfrom,unsigned int, numstates, numbranches);

   // first build the nextstate 
   unsigned int **nextstate;
   CALLOCMATRIX(nextstate,unsigned int,numstates,numbranches);
   unsigned char ins[k];
   unsigned int state;
   unsigned int inp;
   int i;
   unsigned int nextst;
   int *nfrom = new int[numstates];

   for(state = 0; state < numstates; state++) {
      for(inp = 0; inp < numbranches; inp++)  {
         encoder.setstate(state);
         // convert inp to array
         for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
         encoder.encode(ins);
         nextst = encoder.getstate();
         prevstate[nextst][nfrom[nextst]] = state;
         inputfrom[nextst][nfrom[nextst]] = inp;
         nfrom[nextst]++;
      }
   }
      
   for(state = 0; state < numstates; state++) { // for each state
      for(inp = 0; inp < numbranches; inp++)  {
         encoder.setstate(state);
         // convert inp to array
         for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
         encoder.encode(ins);
         nextstate[state][inp] = encoder.getstate();
      }
   }
   for(state = 0; state < numstates; state++) {
      unsigned char *outs;
      cout << "state=" << state << ": ";
      for(inp=0; inp < numbranches; inp++) {
         cout << nextstate[state][inp] << " ";
         encoder.setstate(state);
         for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
         outs = encoder.encode(ins);
         cout << "(";
         for(i = 0; i < n; i++) {
            cout << int(outs[i]);
         }
         cout << ") ";
      }
      cout << endl;
   }
   encoder.setstate(savestate);

   // print the prevstate table
   cout << "fromstates: " << endl;
   for(state = 0; state < numstates; state++) {
      cout << "state=" << state << ": ";
      for(inp=0; inp < (1<<k); inp++) {
         cout << prevstate[state][inp] << " (";
         cout << inputfrom[state][inp] << ") ";
      }
      cout << endl;
   }
   delete[] nfrom;
   encoder.setstate(savestate);
   FREEMATRIX(nextstate);
}



void Convdec::setpaths()
{
   int i;

   fpath = 0;
   bpath = 0;
   numbranchdec = 1;
   for(i = 0; i < numstates; i++) {
      metrics1[i] = LARGE;
   }
   metrics1[startstate] = 0;
   metrics = metrics1;
   othermetrics = metrics2;
}


void Convdec::showpaths(void)
// Used to print out the paths in the trellis 
// (for debugging purposes)
{
   // Fill in the blanks ...
}

/*
Local Variables:
compile-command: "g++ -c -g Convdec.cc"
End:
*/


