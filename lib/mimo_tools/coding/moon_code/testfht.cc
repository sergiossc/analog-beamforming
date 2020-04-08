//
//
//  Program: testfht.cc --- test the fast Hadamard transform routine
//
//
//  Todd K. Moon
//
//  Date:  March 24, 2004
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <ostream.h>
void fht(int *F, int m);
void fhtinv(int *F, int m);

main()
{
   int m = 3;
   int n = (1<<m);
   int F[] = {-1,1,1,-1,1,1,-1,1};
   int i;
   cout << "Original:\n";
   for(i = 0; i < n; i++) {
	  cout << F[i] << " ";
   }
   cout << endl;
   

   fht(F,3);
   cout << "transform:\n";
   for(i = 0; i < n; i++) {
	  cout << F[i] << " ";
   }
   cout << endl;
   cout << endl;
}



/*
Local Variables:
compile-command: "g++ -o testfht testfht.cc fht.cc"
End:
*/


