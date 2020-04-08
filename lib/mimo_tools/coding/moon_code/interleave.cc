// interleave.cc -- A random interleaver
// Todd K. Moon

// This is meant to be a very elementary interleaver, 
// operating simply on the basis of random selection.
// Interleavers requirless less memory (and less setup) 
// could be made (e.g., using a 2-d array).

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only


#include "interleave.h"
#include <iostream>
using namespace std;
extern "C" {
#include <stdlib.h>
}

interleave::interleave(int in_size, unsigned int seed)
{
   int i,j,r,numleft;
   size = in_size;
   pi = new int[size];
   piinv = new int[size];
   double rm;
   // use a brute-force approach to selection without replacement
   int *notused = new int[size];
   for(i = 0; i < size; i++) {
	  notused[i] = i;
   }
   rm = double(RAND_MAX);
   rm = rm+1;
   numleft = size;
   if(seed) srand(seed); else srand(1);
   for(i = 0; i < size; i++) {
	  r = int(numleft*(double(rand())/rm));
	  pi[i] = notused[r];
	  piinv[pi[i]] = i;
	  // get rid of the entry used
	  for(j = r; j < numleft-1; j++) {
		 notused[j] = notused[j+1];
	  }
	  numleft--;
   }   
   // check it
//    for(i = 0; i < size; i++) {
// 	  cout << pi[i] << " " << piinv[i] << " " << pi[piinv[i]] << " " <<
// 		 piinv[pi[i]] << endl;
//    }
   delete[] notused;
}

void
interleave::Pi(double *in, double *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[pi[i]];
   }
}

void
interleave::Pi(double **in, double **out, int nrow)
{
   int i,j;
   for(i = 0; i < size; i++) {
	  for(j = 0; j < nrow; j++) {
		 out[i][j] = in[pi[i]][j];
	  }
   }
}

void
interleave::Pi(unsigned char *in, unsigned char *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[pi[i]];
   }
}

void
interleave::Piinv(double *in, double *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[piinv[i]];
   }
}

void
interleave::Piinv(double **in, double **out, int nrow)
{
   int i,j;
   for(i = 0; i < size; i++) {
	  for(j = 0; j < nrow; j++) {
		 out[i][j] = in[piinv[i]][j];
	  }
   }
}

void
interleave::PiinvTimesoverlay(double **in, double **out, int nrow)
{
   int i,j;
   for(i = 0; i < size; i++) {
	  for(j = 0; j < nrow; j++) {
		 out[i][j] *= in[piinv[i]][j];
	  }
   }
}

void
interleave::Piinv(unsigned char *in, unsigned char *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[piinv[i]];
   }
}


/*
Local Variables:
compile-command: "g++ -c interleave.cc"
End:
*/
