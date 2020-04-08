//
//
//  Program:  fht.cc
//  Compute the fast Hadamard transform on the integer sequence F
//  where F has 2^m points in it.
//
//  Todd K. Moon
//  Utah State University
//
//  Date:  March 24, 2004
//
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only


void fht(int *F, int m)
{
   int i;
   int n;
   int skip,skip2;
   int j, j1, j1skip, k1;
   int tmp;

   n = (1<<m);					// number of points
   skip = 1;
   skip2 = 2;
   for(i = 0; i < m; i++) {		// over each iteration
	  j = 0;					// start at the first line
	  while(j < n) {
		 j1 = j;
		 j1skip = j1+skip;
		 for(k1 = 0, j1=j,j1skip=j1+skip; k1 < skip; k1++,j1++,j1skip++) {
			tmp = F[j1];
			F[j1] += F[j1skip];
			F[j1skip] = tmp - F[j1skip];
		 }
		 j += skip2;
	  }
	  skip = skip2;
	  skip2 <<= 1;
   }
}


/*
Local Variables:
compile-command: "g++ -o testfht testfht.cc fht.cc"
End:
*/


